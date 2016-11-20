#ifndef PERLIN_H
#define PERLIN_H

#include "Vector.h"

#include <cstdint>

class Perlin
{
public:
							Perlin		(uint64_t seed = 0);
	virtual					~Perlin		(void);

	double					Noise		(const Vector3& point);
	double					Noise		(const Vector2& point);

	inline static double	Fade		(double x) { return x * x * x * (x * (x * 6.0 - 15.0) + 10.0); }
	inline static double	Lerp		(double t, double a, double b) { return a + t * (b - a); }
	static double			Gradient	(int hash, const Vector3& point);
	static double			Gradient	(int hash, const Vector2& point);

private:

	int m_permutations[512];

	static const int s_permutations[256];
};

#endif // PERLIN_H
