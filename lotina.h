#ifndef LOTINA_H
#define LOTINA_H

#define LTN_GREATEST(a, b) ((a) > (b) ? (a) : (b))
#define LTN_SMALLEST(a, b) ((a) < (b) ? (a) : (b))

struct Ltn_V2i {
	int x, y;
};

struct Ltn_V2f {
	float x, y;
};

static float ltn_dot_v2f(struct Ltn_V2f a, struct Ltn_V2f b)
{ return a.x*b.x + a.y*b.y; }

struct Ltn_Fluid_Domain {
	int i;
};

#endif
