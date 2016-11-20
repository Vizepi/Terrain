#include "Matrix.h"

#include <cmath>

inline Matrix2x2::Matrix2x2(void)
{
	m_matrix[0] = m_matrix[1] = m_matrix[2] = m_matrix[3] = 0;
}

inline Matrix2x2::Matrix2x2(double m0, double m1, double m2, double m3)
{
	m_matrix[0] = m0;
	m_matrix[1] = m1;
	m_matrix[2] = m2;
	m_matrix[3] = m3;
}

inline Matrix2x2::Matrix2x2(double angle)
{
	m_matrix[0] = m_matrix[3] = cos(angle);
	m_matrix[1] = sin(angle);
	m_matrix[2] = -m_matrix[1];
}

inline Matrix2x2::Matrix2x2(double m[4])
{
	m_matrix[0] = m[0];
	m_matrix[1] = m[1];
	m_matrix[2] = m[2];
	m_matrix[3] = m[3];
}

inline Matrix2x2::~Matrix2x2(void)
{

}

inline const Vector2 Matrix2x2::operator*(const Vector2& vector) const
{
	return Vector2(
				m_matrix[0] * vector.X() + m_matrix[2] * vector.Y(),
				m_matrix[1] * vector.X() + m_matrix[3] * vector.Y());
}

inline Matrix3x3::Matrix3x3(void)
{
	m_matrix[0] = m_matrix[1] = m_matrix[2] =
	m_matrix[3] = m_matrix[4] = m_matrix[5] =
	m_matrix[6] = m_matrix[7] = m_matrix[8] = 0;
}

inline Matrix3x3::Matrix3x3(double m0, double m1, double m2, double m3, double m4, double m5, double m6, double m7, double m8)
{
	m_matrix[0] = m0;
	m_matrix[1] = m1;
	m_matrix[2] = m2;
	m_matrix[3] = m3;
	m_matrix[4] = m4;
	m_matrix[5] = m5;
	m_matrix[6] = m6;
	m_matrix[7] = m7;
	m_matrix[8] = m8;
}

inline Matrix3x3::Matrix3x3(const Vector3& direction, double angle)
{
	double c = cos(angle);
	double c1 = 1.0 - c;
	double s = sin(angle);

	m_matrix[0] = direction.X() * direction.X() * c1 + c;
	m_matrix[1] = direction.X() * direction.Y() * c1 + direction.Z() * s;
	m_matrix[2] = direction.X() * direction.Z() * c1 - direction.Y() * s;
	m_matrix[3] = direction.X() * direction.Y() * c1 - direction.Z() * s;
	m_matrix[4] = direction.Y() * direction.Y() * c1 + c;
	m_matrix[5] = direction.Y() * direction.Z() * c1 + direction.X() * s;
	m_matrix[6] = direction.X() * direction.Z() * c1 + direction.Y() * s;
	m_matrix[7] = direction.Y() * direction.Z() * c1 - direction.X() * s;
	m_matrix[8] = direction.Z() * direction.Z() * c1 + c;
}

inline Matrix3x3::Matrix3x3(double m[9])
{
	m_matrix[0] = m[0];
	m_matrix[1] = m[1];
	m_matrix[2] = m[2];
	m_matrix[3] = m[3];
	m_matrix[4] = m[4];
	m_matrix[5] = m[5];
	m_matrix[6] = m[6];
	m_matrix[7] = m[7];
	m_matrix[8] = m[8];
}

inline Matrix3x3::~Matrix3x3(void)
{

}

inline const Vector3 Matrix3x3::operator*(const Vector3& vector) const
{
	return Vector3(
			m_matrix[0] * vector.X() + m_matrix[3] * vector.Y() + m_matrix[6] * vector.Z(),
			m_matrix[1] * vector.X() + m_matrix[4] * vector.Y() + m_matrix[7] * vector.Z(),
			m_matrix[2] * vector.X() + m_matrix[5] * vector.Y() + m_matrix[8] * vector.Z());
}
