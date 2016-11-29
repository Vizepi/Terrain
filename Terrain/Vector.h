#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

class Vector2
{

public:

	inline			Vector2		(void) :				m_x(0.0), m_y(0.0) {}
	inline			Vector2		(double x, double y) :	m_x(x), m_y(y) {}
	inline			Vector2		(const Vector2& rhs) :	m_x(rhs.m_x), m_y(rhs.m_y) {}
	inline virtual	~Vector2	(void)					{}

	inline void		operator=	(const Vector2& rhs)	{ m_x = rhs.m_x; m_y = rhs.m_y; }

	inline Vector2	Orthogonal	(void) const			{ return Vector2(-m_y, m_x); }
	inline double	SquareLength(void) const			{ return m_x * m_x + m_y * m_y; }
	inline double	Length		(void) const			{ return sqrt(SquareLength()); }
	inline void		Normalize	(void)					{ double length = Length(); if(length == 0.0) return; m_x /= length; m_y /= length; }

	inline double	X			(void) const			{ return m_x; }
	inline double	Y			(void) const			{ return m_y; }

	inline void		SetX		(double x)				{ m_x = x; }
	inline void		SetY		(double y)				{ m_y = y; }

protected:

	double m_x, m_y;

};

inline double	Product		(const Vector2& lhs, const Vector2& rhs) { return lhs.X() * rhs.Y() - lhs.Y() * rhs.X(); }
inline double	DotProduct	(const Vector2& lhs, const Vector2& rhs) { return lhs.X() * rhs.X() + lhs.Y() * rhs.Y(); }
inline Vector2	operator+	(const Vector2& lhs, const Vector2& rhs) { return Vector2(lhs.X() + rhs.X(), lhs.Y() + rhs.Y()); }
inline Vector2	operator-	(const Vector2& lhs, const Vector2& rhs) { return Vector2(lhs.X() - rhs.X(), lhs.Y() - rhs.Y()); }
inline Vector2	operator*	(const Vector2& lhs, const Vector2& rhs) { return Vector2(lhs.X() * rhs.X(), lhs.Y() * rhs.Y()); }
inline Vector2	operator/	(const Vector2& lhs, const Vector2& rhs) { return Vector2(lhs.X() / rhs.X(), lhs.Y() / rhs.Y()); }
inline Vector2	operator+	(const Vector2& lhs, double rhs) { return Vector2(lhs.X() + rhs, lhs.Y() + rhs); }
inline Vector2	operator-	(const Vector2& lhs, double rhs) { return Vector2(lhs.X() - rhs, lhs.Y() - rhs); }
inline Vector2	operator*	(const Vector2& lhs, double rhs) { return Vector2(lhs.X() * rhs, lhs.Y() * rhs); }
inline Vector2	operator/	(const Vector2& lhs, double rhs) { return Vector2(lhs.X() / rhs, lhs.Y() / rhs); }
inline bool		operator<	(const Vector2& lhs, const Vector2& rhs) { return lhs.X() < rhs.X() ? true : (lhs.X() == rhs.X() ? lhs.Y() < rhs.Y() : false); }

class Vector3
{

public:

	inline			Vector3		(void) :							m_x(0.0), m_y(0.0), m_z(0.0) {}
	inline			Vector3		(double x, double y, double z) :	m_x(x), m_y(y), m_z(z) {}
	inline			Vector3		(const Vector3& rhs) :				m_x(rhs.m_x), m_y(rhs.m_y), m_z(rhs.m_z) {}
	inline virtual	~Vector3	(void)								{}

	inline void		operator=	(const Vector3& rhs)				{ m_x = rhs.m_x; m_y = rhs.m_y; m_z = rhs.m_z; }

	inline double	SquareLength(void) const						{ return m_x * m_x + m_y * m_y + m_z * m_z; }
	inline double	Length		(void) const						{ return sqrt(SquareLength()); }
	inline void		Normalize	(void)								{ double length = Length(); if (length == 0.0) return; m_x /= length; m_y /= length; m_z /= length; }

	inline double	X			(void) const						{ return m_x; }
	inline double	Y			(void) const						{ return m_y; }
	inline double	Z			(void) const						{ return m_z; }

	inline void		SetX		(double x)							{ m_x = x; }
	inline void		SetY		(double y)							{ m_y = y; }
	inline void		SetZ		(double z)							{ m_z = z; }

protected:

	double m_x, m_y, m_z;

};

inline Vector3	Product		(const Vector3& lhs, const Vector3& rhs) { return Vector3(lhs.Y() * rhs.Z() - lhs.Z() * rhs.Y(), lhs.Z() * rhs.X() - lhs.X() * rhs.Z(), lhs.X() * rhs.Y() - lhs.Y() * rhs.X()); }
inline double	DotProduct	(const Vector3& lhs, const Vector3& rhs) { return lhs.X() * rhs.X() + lhs.Y() * rhs.Y() + lhs.Z() * rhs.Z(); }
inline Vector3	operator+	(const Vector3& lhs, const Vector3& rhs) { return Vector3(lhs.X() + rhs.X(), lhs.Y() + rhs.Y(), lhs.Z() + rhs.Z()); }
inline Vector3	operator-	(const Vector3& lhs, const Vector3& rhs) { return Vector3(lhs.X() - rhs.X(), lhs.Y() - rhs.Y(), lhs.Z() - rhs.Z()); }
inline Vector3	operator*	(const Vector3& lhs, const Vector3& rhs) { return Vector3(lhs.X() * rhs.X(), lhs.Y() * rhs.Y(), lhs.Z() * rhs.Z()); }
inline Vector3	operator/	(const Vector3& lhs, const Vector3& rhs) { return Vector3(lhs.X() / rhs.X(), lhs.Y() / rhs.Y(), lhs.Z() / rhs.Z()); }
inline Vector3	operator+	(const Vector3& lhs, double rhs) { return Vector3(lhs.X() + rhs, lhs.Y() + rhs, lhs.Z() + rhs); }
inline Vector3	operator-	(const Vector3& lhs, double rhs) { return Vector3(lhs.X() - rhs, lhs.Y() - rhs, lhs.Z() - rhs); }
inline Vector3	operator*	(const Vector3& lhs, double rhs) { return Vector3(lhs.X() * rhs, lhs.Y() * rhs, lhs.Z() * rhs); }
inline Vector3	operator/	(const Vector3& lhs, double rhs) { return Vector3(lhs.X() / rhs, lhs.Y() / rhs, lhs.Z() / rhs); }
inline bool		operator<	(const Vector3& lhs, const Vector3& rhs) { return lhs.X() < rhs.X() ? true : (lhs.X() == rhs.X() ? (lhs.Y() < rhs.Y() ? true : (lhs.Y() == rhs.Y() ? lhs.Z() < rhs.Z() : false)) : false); }

#endif // VECTOR_H
