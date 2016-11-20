#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.h"

class Matrix2x2
{

public:

	inline					Matrix2x2	(void);
	inline					Matrix2x2	(double m0, double m1, double m2, double m3);
	inline					Matrix2x2	(double angle);
	inline					Matrix2x2	(double m[4]);
	virtual inline			~Matrix2x2	(void);

	inline const Vector2	operator*	(const Vector2& vector) const;

private:

	double m_matrix[4];

};

class Matrix3x3
{

public:

	inline					Matrix3x3	(void);
	inline					Matrix3x3	(double m0, double m1, double m2, double m3, double m4, double m5, double m6, double m7, double m8);
	inline					Matrix3x3	(const Vector3& direction, double angle);
	inline					Matrix3x3	(double m[9]);
	virtual inline			~Matrix3x3	(void);

	inline const Vector3	operator*	(const Vector3& vector) const;

private:

	double m_matrix[9];

};

#include "Matrix.inl"

#endif // MATRIX_HPP
