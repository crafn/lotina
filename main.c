#include "lotina.h"

#include <assert.h>
#include <math.h>
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

const int point_mirrored_dir_index[9] = {
	0,
	3, 4, 1, 2,
	7, 8, 5, 6
};

/* Prefers vertical mirroring */
const int side_mirrored_dir_index[9] = {
	0,
	3, 4, 1, 2,
	8, 7, 6, 5
};

char is_domain_edge(int x, int y)
{
	return x == 0 || y == 0 || x + 1 == DOMAIN_WIDTH || y + 1 == DOMAIN_HEIGHT;
}

int main(void)
{
	float f[DOMAIN_SIZE][9] = {{0}};
	float f_tmp[DOMAIN_SIZE][9] = {{0}};
	float rho[DOMAIN_SIZE] = {0};
	float psi[DOMAIN_SIZE] = {0}; /* Cohesion */
	float psi2[DOMAIN_SIZE] = {0};
	struct Ltn_V2f u[DOMAIN_SIZE] = {{0, 0}};
	char boundary[DOMAIN_SIZE] = {0};
	int i, k, x, y;
	int step = 0;
	const float c = 1.01f;
	const float tau = 0.9f;
	const float dx = 1.0f;
	char * map = NULL;

	{ /* Load map (leaks) */
		long length;
		FILE * f = fopen("map", "rb");
		assert(f && "Where's map?");
		fseek(f, 0, SEEK_END);
		length = ftell(f);
		assert(length == (DOMAIN_WIDTH + 1)*DOMAIN_HEIGHT);
		fseek(f, 0, SEEK_SET);
		map = malloc(length);
		fread(map, 1, length, f);
		fclose(f);
	}

	system("clear");


	/* Initial conditions */
	for (y = 0; y < DOMAIN_HEIGHT; ++y) {
		for (x = 0; x < DOMAIN_WIDTH; ++x) {
			char ch = map[x + (DOMAIN_WIDTH + 1)*y];

			if (ch == 'w')
				for (k = 0; k < 9; ++k)
					f[F_IX(x,y)][k] =  w[k]*20;
			else if (ch == '[')
				boundary[F_IX(x,y)] = 1;

			if (is_domain_edge(x, y))
				boundary[F_IX(x,y)] = 1;

		}
	}

	while (1) {
		float mass_check = 0.0f, mass_check_2 = 0.0f;

		/* Update some fields */
		float total_mass = 0.0;
		float max_u = 0.0;
		float max_rho = 0.0;
		for (i = 0; i < DOMAIN_SIZE; ++i) {
			rho[i] = 0.0;
			for (k = 0; k < 9; ++k)
				rho[i] += f[i][k];
			total_mass += rho[i]*dx*dx;

			u[i].x = u[i].y = 0.0;
			for (k = 0; k < 9; ++k) {
				u[i].x += v_e[k].x*f[i][k];
				u[i].y += v_e[k].y*f[i][k];
			}
			u[i].x /= (rho[i] + 0.00001);
			u[i].y /= (rho[i] + 0.00001);

			if (ltn_dot_v2f(u[i], u[i]) > max_u*max_u)
				max_u = sqrt(ltn_dot_v2f(u[i], u[i]));

			if (rho[i] > max_rho)
				max_rho = rho[i];

			/* Multiplier determines cohesion strength */
			psi[i] = 2.5*exp(-0.7/rho[i]);

			psi2[i] = rho[i]*0.0;
		}

		/* Collision step */
		for (y = 0; y < DOMAIN_HEIGHT; ++y) {
			for (x = 0; x < DOMAIN_WIDTH; ++x) {
				float new_rho = 0.0f;
				int i = F_IX(x,y);
				struct Ltn_V2f force = {0, 0};
				struct Ltn_V2f u_xy = u[i];
				float total_feq = 0.0;

				{ /* Gravity */
					const struct Ltn_V2f g = {
						0.0, 0.1
					};
					float m = dx*dx*rho[i];
					force.x += g.x*m;
					force.y += g.y*m;
				}

				if (!is_domain_edge(x,y)) { /* Cohesion */
					for (k = 1; k < 9; ++k) {
						force.x += psi[i]*w[k]*psi[i + index_shift[k]]*v_e[k].x;
						force.y += psi[i]*w[k]*psi[i + index_shift[k]]*v_e[k].y;

						force.x += w[k]*(psi2[i] - psi2[i + index_shift[k]])*v_e[k].x;
						force.y += w[k]*(psi2[i] - psi2[i + index_shift[k]])*v_e[k].y;
					}
				}

				if (1){ /* Apply external forces */
					u_xy.x += tau*force.x/(rho[i] + 0.000001);
					u_xy.y += tau*force.y/(rho[i] + 0.000001);
				}

				for (k = 0; k < 9; ++k) {
					float f_k = f[i][k];
					float f_eq;

					{ /* BGK */
						/* This produces negative values sometimes. Maybe too large velocities. */
						float a1, a2, a3;
						a1 = 3.0f*ltn_dot_v2f(v_e[k], u_xy)/c;
						a2 = 9.0f/2*ltn_dot_v2f(v_e[k], u_xy)*ltn_dot_v2f(v_e[k], u_xy)/(c*c);
						a3 = -3.0f/2*ltn_dot_v2f(u_xy, u_xy)/(c*c);
						f_eq = rho[i]*w[k]*(1 + a1 + a2 + a3);
						total_feq += f_eq;
					}

					{
						float new_f_k = LTN_GREATEST(f_k - (f_k - f_eq)/tau, 0);
						f_tmp[i][k] = new_f_k;
						new_rho += f_tmp[i][k];
					}
				}

#if 1
				/* HACK -- preserve mass */
				for (k = 0; k < 9; ++k)
					f_tmp[i][k] *= (rho[i] + 0.0000001)/(new_rho + 0.00000001);
#endif
			}
		}

		for (i = 0; i < DOMAIN_SIZE; ++i) {
			float rho = 0.0f;
			for (k = 0; k < 9; ++k)
				rho += f_tmp[i][k];
			mass_check += rho*dx*dx;
		}

		/* Streaming step */
		for (i = 0; i < DOMAIN_SIZE; ++i)
			for (k = 0; k < 9; ++k)
				f[i][k] = 0; /* This shouldn't be necessary, but otherwise mass grows o_O */
		for (y = 1; y < DOMAIN_HEIGHT - 1; ++y) {
			for (x = 1; x < DOMAIN_WIDTH - 1; ++x) {
				int i = F_IX(x,y);
				for (k = 0; k < 9; ++k) {
					int to_i = i + index_shift[k];
					f[to_i][k] = f_tmp[i][k];
					f_tmp[i][k] = 0.0; /* Shouldn't be necessary */
				}
			}
		}

		/* Mirror from boundaries */
		for (y = 0; y < DOMAIN_HEIGHT; ++y) {
			for (x = 0; x < DOMAIN_WIDTH; ++x) {
				int i = F_IX(x,y);

				if (boundary[i] == 0)
					continue;

				for (k = 0; k < 9; ++k) {
					int from_i;
					if (f[i][k] == 0.0f)
						continue; /* This is scary, as it protects from out-of-bounds access */

					from_i = i - index_shift[k];
					assert(boundary[from_i] == 0);
					f[from_i][side_mirrored_dir_index[k]] += f[i][k];
					f[i][k] = 0.0f;
				}
			}
		}

		for (i = 0; i < DOMAIN_SIZE; ++i) {
			float rho = 0.0f;
			for (k = 0; k < 9; ++k)
				rho += f[i][k];
			mass_check_2 += rho*dx*dx;
		}

		if ((step++) % 2 == 0) {
			/* Draw domain */
			clear_screen();
			for (y = 0; y < DOMAIN_HEIGHT; ++y) {
				for (x = 0; x < DOMAIN_WIDTH; ++x) {
					float rho_xy = rho[F_IX(x,y)];

					if (boundary[F_IX(x,y)] == 1)
						draw_char(x, y, '[');

					if (rho_xy > 0.001) {
						char ch;
						if (rho_xy < 0.3)
							ch = '.';
						else if (rho_xy < 0.6)
							ch = '~';
						else if (rho_xy < 10.0)
							ch = '*';
						else
							ch = '#';

						draw_char(x, y, ch);
					}
				}
			}
			blit_screen();

			printf("Fluid mass: %f\n", total_mass);
			printf("dif: %f\n", mass_check - total_mass);
			printf("dif2: %f\n", mass_check_2 - mass_check);
			printf("max rho: %f\n", max_rho);
			printf("max u: %f\n", max_u);
#if 1
			usleep(1000000/60);
#else
			if (getchar() == 'q')
				break;
#endif
		}
	}
	return 0;
}
