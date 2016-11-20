#include "Perlin.h"

#include <cmath>
#include <random>
#include <algorithm>

/*static*/ const int Perlin::s_permutations[256] =
{
	151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
	140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
	247,120,234, 75,  0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
	 57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
	 74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
	 60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
	 65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
	200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
	 52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
	207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
	119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
	129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
	218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
	 81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
	184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
	222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
};

Perlin::Perlin(uint64_t seed)
{
	std::mt19937_64 generator(seed);
	for(uint64_t i = 0; i < 256; ++i)
	{
		m_permutations[i] = s_permutations[i];
	}
	std::shuffle(&(m_permutations[0]), &(m_permutations[255]), generator);
	for(uint64_t i = 0; i < 256; ++i)
	{
		m_permutations[i + 256] = m_permutations[i];
	}
}

/*virtual*/ Perlin::~Perlin(void)
{

}

double Perlin::Noise(const Vector3& point)
{
	double x = point.X();
	double y = point.Y();
	double z = point.Z();
	int X = int(floor(x)) & 0xff;
	int Y = int(floor(y)) & 0xff;
	int Z = int(floor(z)) & 0xff;

	x -= floor(x);
	y -= floor(y);
	z -= floor(z);

	double u = Fade(x);
	double v = Fade(y);
	double w = Fade(z);

	int _A = m_permutations[X     ] + Y;
	int AA = m_permutations[_A    ] + Z;
	int AB = m_permutations[_A + 1] + Z;
	int _B = m_permutations[X + 1 ] + Y;
	int BA = m_permutations[_B    ] + Z;
	int BB = m_permutations[_B + 1] + Z;

	return Lerp(w,	Lerp(v,	Lerp(u,	Gradient(m_permutations[AA  ], Vector3(x  , y  , z  )),
									Gradient(m_permutations[BB  ], Vector3(x-1, y  , z  ))),
							Lerp(u, Gradient(m_permutations[AB  ], Vector3(x  , y-1, z  )),
									Gradient(m_permutations[BB  ], Vector3(x-1, y-1, z  )))),
					Lerp(v,	Lerp(u, Gradient(m_permutations[AA+1], Vector3(x  , y  , z-1)),
									Gradient(m_permutations[BA+1], Vector3(x-1, y  , z-1))),
							Lerp(u, Gradient(m_permutations[AB+1], Vector3(x  , y-1, z-1)),
									Gradient(m_permutations[BB+1], Vector3(x-1, y-1, z-1)))));
}

double Perlin::Noise(const Vector2& point)
{
	double x = point.X();
	double y = point.Y();
	int X = int(floor(x)) & 0xff;
	int Y = int(floor(y)) & 0xff;

	x -= floor(x);
	y -= floor(y);

	double u = Fade(x);
	double v = Fade(y);

	int A = m_permutations[X    ] + Y;
	int B = m_permutations[X + 1] + Y;

	return Lerp(v,	Lerp(u,	Gradient(m_permutations[A  ], Vector2(x  , y  )),
							Gradient(m_permutations[B  ], Vector2(x-1, y  ))),
					Lerp(u,	Gradient(m_permutations[A+1], Vector2(x  , y-1)),
							Gradient(m_permutations[B+1], Vector2(x-1, y-1))));
}

/*static*/ double Perlin::Gradient(int hash, const Vector3& point)
{
	int h = hash & 0xf;
	double u = (h < 8) ? point.X() : point.Y();
	double v = (h < 4) ? point.Y() : ((h == 12 || h == 14) ? point.X() : point.Z());
	return ((h & 0x1) == 0 ? u : -u) + ((h & 0x2) == 0 ? v : -v);
}

/*static*/ double Perlin::Gradient(int hash, const Vector2& point)
{
	int h = hash & 0x3;
	return ((h & 0x1) ? -point.X() : point.X()) + ((h & 0x2) ? -point.Y() : point.Y());
}
