#include "lotina.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void usleep(unsigned int);

#define LINUX

#define DOMAIN_WIDTH 100
#define DOMAIN_HEIGHT 50
#define DOMAIN_SIZE (DOMAIN_WIDTH*DOMAIN_HEIGHT)
#define F_IX(x, y) ((x) + (y)*DOMAIN_WIDTH)
char g_screenbuf[(DOMAIN_WIDTH + 1)*DOMAIN_HEIGHT];

void draw_char(int x, int y, char ch)
{
	if (x < 0 || x >= DOMAIN_WIDTH)
		return;
	if (y < 0 || y >= DOMAIN_HEIGHT)
		return;

	g_screenbuf[x + y*(DOMAIN_WIDTH + 1)] = ch;
}

void clear_screen(void)
{
	int i;
	for (i = 0; i < DOMAIN_SIZE; ++i)
		g_screenbuf[i] = ' ';
}

void blit_screen(void)
{
	int i;
	for (i = 0; i < DOMAIN_HEIGHT; ++i)
		g_screenbuf[i*(DOMAIN_WIDTH + 1) + DOMAIN_WIDTH] = '\n';
	printf("\033[0;0H");
	printf("%.*s", (DOMAIN_WIDTH + 1)*DOMAIN_HEIGHT, g_screenbuf);
}

/* d3q9 LBM */
const struct Ltn_V2f v_e[9] = {
	{0, 0},
	{1, 0}, {0, 1}, {-1, 0}, {0, -1},
	{1, 1}, {-1, 1}, {-1, -1}, {1, -1}
};

const float w[9] = {
	4.0f/9,
	1.0f/9, 1.0f/9, 1.0f/9, 1.0f/9,
	1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36,
};

const int index_shift[9] = {
	0,
	1, DOMAIN_WIDTH, -1, -DOMAIN_WIDTH,
	1 + DOMAIN_WIDTH, -1 + DOMAIN_WIDTH, -1 - DOMAIN_WIDTH, 1 - DOMAIN_WIDTH
};

const int mirrored_dir_index[9] = {
	0,
	3, 4, 1, 2,
	7, 8, 5, 6
};

int main(void)
{
	float f[DOMAIN_SIZE][9] = {{0}};
	float f_tmp[DOMAIN_SIZE][9] = {{0}};
	float rho[DOMAIN_SIZE] = {0};
	struct Ltn_V2f u[DOMAIN_SIZE] = {{0, 0}};
	char boundary[DOMAIN_SIZE] = {0};
	int i, k, x, y;
	const float c = 1.0f;
	const float tau = 2.5f;
	const float dx = 1.0f;
	float rho_0 = 1.0f;

	/* Initial conditions */
	for (y = 0; y < DOMAIN_HEIGHT; ++y) {
		for (x = 0; x < DOMAIN_WIDTH; ++x) {
			if (	x > 70 && x + 1 < DOMAIN_WIDTH &&
					y > 2 && y < 10 && y + 1 < DOMAIN_HEIGHT) {
				for (k = 0; k < 9; ++k)
					f[F_IX(x,y)][k] = w[k]*1.0;
			}

			if (x == 0 || y == 0 || x + 1 == DOMAIN_WIDTH || y + 1 == DOMAIN_HEIGHT)
				boundary[F_IX(x,y)] = 1;

			if (x == 15 && y > 10)
				boundary[F_IX(x,y)] = 1;

			if (x == 70  && y > 25)
				boundary[F_IX(x,y)] = 1;

			if (y == 20 && x > 25)
				boundary[F_IX(x,y)] = 1;
		}
	}

	while (1) {
		float mass_check = 0.0f, mass_check_2 = 0.0f;

		/* Update density and velocity fields */
		float total_mass = 0.0;
		for (i = 0; i < DOMAIN_SIZE; ++i) {
			rho[i] = 0.0;
			for (k = 0; k < 9; ++k)
				rho[i] += f[i][k];
			total_mass += rho[i]*dx*dx;

			u[i].x = u[i].y = 0.0;
			for (k = 0; k < 9; ++k) {
				u[i].x += c*v_e[k].x*f[i][k];
				u[i].y += c*v_e[k].y*f[i][k];
			}
			u[i].x /= (rho[i] + 0.0001);
			u[i].y /= (rho[i] + 0.0001);
		}

		/* Collision step */
		for (y = 0; y < DOMAIN_HEIGHT; ++y) {
			for (x = 0; x < DOMAIN_WIDTH; ++x) {
				float new_rho = 0.0f;
				int i = F_IX(x,y);
				for (k = 0; k < 9; ++k) {
					struct Ltn_V2f u_xy = u[i];
					float f_k = f[i][k];
					float f_eq, f_force;
					{ /* BGK */
						float s_k = w[k]*(	3.0f*ltn_dot_v2f(v_e[k], u_xy)/c +
											9.0f/2*ltn_dot_v2f(v_e[k], u_xy)*ltn_dot_v2f(v_e[k], u_xy)/(c*c) +
											-3.0f/2*ltn_dot_v2f(u_xy, u_xy)/(c*c));
						f_eq = rho[i]*w[k] + rho_0*s_k;
					}

					{ /* Apply external forces */
						const struct Ltn_V2f g = {
							0.0, 1
						};
						float m = rho[i]*dx*dx;
						struct Ltn_V2f f;
						f.x = g.x*m;
						f.y = g.y*m;

						/* @todo */
						f_force = w[k]*ltn_dot_v2f(v_e[k], f);
						f_force = 0.0;
					}

					{
						float new_f_k = f_k - (f_k - f_eq)/tau + f_force;
						f_tmp[i][k] = LTN_GREATEST(new_f_k, 0.0);
						new_rho += f_tmp[i][k];
					}
				}

				/* HACK -- preserve mass */
				for (k = 0; k < 9; ++k)
					f_tmp[i][k] *= (rho[i] + 0.0000001)/(new_rho + 0.00000001);
			}
		}

		for (i = 0; i < DOMAIN_SIZE; ++i) {
			float rho = 0.0f;
			for (k = 0; k < 9; ++k)
				rho += f_tmp[i][k];
			mass_check += rho*dx*dx;
		}

		/* Streaming step */
		for (y = 1; y < DOMAIN_HEIGHT - 1; ++y) {
			for (x = 1; x < DOMAIN_WIDTH - 1; ++x) {
				int i = F_IX(x,y);
				for (k = 0; k < 9; ++k)
					f[i][k] = 0.0f;

				for (k = 0; k < 9; ++k) {
					int from_i = i - index_shift[k];
#if 1
					if (boundary[i]) {
						f[i][k] = 0.0f;
						continue;
					}

					if (boundary[from_i] == 0) /* Free flow */
						f[i][k] = f_tmp[from_i][k];
					else { /* Mirror flow from boundary */
						f[i][k] = f_tmp[i][mirrored_dir_index[k]];
					}
#else
					f[i][k] = f_tmp[from_i][k];
#endif
				}
			}
		}

		for (i = 0; i < DOMAIN_SIZE; ++i) {
			float rho = 0.0f;
			for (k = 0; k < 9; ++k)
				rho += f[i][k];
			mass_check_2 += rho*dx*dx;
		}

		clear_screen();

		/* Draw domain */
		for (y = 0; y < DOMAIN_HEIGHT; ++y) {
			for (x = 0; x < DOMAIN_WIDTH; ++x) {
				float rho_xy = rho[F_IX(x,y)];

				if (boundary[F_IX(x,y)] == 1)
					draw_char(x, y, '#');

				if (rho_xy > 0.01) {
					char ch;
					if (rho_xy < 0.3)
						ch = '.';
					else if (rho_xy < 0.6)
						ch = '~';
					else
						ch = '*';

					draw_char(x, y, ch);
				}
			}
		}

		blit_screen();
		printf("Fluid mass: %f\n", total_mass);
		printf("dif: %f\n", mass_check - total_mass);
		printf("dif2: %f\n", mass_check_2 - mass_check);

#if defined(LINUX)
		usleep(1000000/60);
#else
		if (getchar() == 'q')
			break;
#endif
	}
	return 0;
}
