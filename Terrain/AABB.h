#ifndef AABB_H
#define AABB_H

#include "Vector.h"

class AABB2
{

public:

	inline					AABB2		(void) : m_a(), m_b() {}
	inline					AABB2		(const Vector2& a, const Vector2& b) : m_a(a), m_b(b) { SolveOrder(); }
	inline					~AABB2		(void) {}

	inline const Vector2&	A			(void) const { return m_a; }
	inline const Vector2&	B			(void) const { return m_b; }
	inline Vector2			Size		(void) const { return Vector2(m_b.X() - m_a.X(), m_b.Y() - m_a.Y()); }

	inline void				SetA		(const Vector2& a) { m_a = a; SolveOrder(); }
	inline void				SetB		(const Vector2& b) { m_b = b; SolveOrder(); }

	inline double			Area		(void) const { return (m_b.X() - m_a.X()) * (m_b.Y() - m_a.Y()); }
	inline double			Perimeter	(void) const { return 2.0 * ((m_b.X() - m_a.X()) + (m_b.Y() - m_a.Y())); }
	inline bool				Contains	(const Vector2& point) { return !(point.X() < m_a.X() || point.X() > m_b.X() || point.Y() < m_a.Y() || point.Y() > m_b.Y()); }

private:

	inline void SolveOrder(void)
	{
		if(m_a.X() > m_b.X())
		{
			double tmp = m_a.X(); m_a.SetX(m_b.X()); m_b.SetX(tmp);
		}
		if(m_a.Y() > m_b.Y())
		{
			double tmp = m_a.Y(); m_a.SetY(m_b.Y()); m_b.SetY(tmp);
		}
	}

	Vector2 m_a, m_b;

};

class AABB3
{

public:

	inline					AABB3		(void) : m_a(), m_b() {}
	inline					AABB3		(const Vector3& a, const Vector3& b) : m_a(a), m_b(b) { SolveOrder(); }
	inline					~AABB3		(void) {}

	inline const Vector3&	A			(void) const { return m_a; }
	inline const Vector3&	B			(void) const { return m_b; }
	inline Vector3			Size		(void) const { return Vector3(m_b.X() - m_a.X(), m_b.Y() - m_a.Y(), m_b.Z() - m_a.Z()); }

	inline void				SetA		(const Vector3& a) { m_a = a; SolveOrder(); }
	inline void				SetB		(const Vector3& b) { m_b = b; SolveOrder(); }

	inline double			Volume		(void) const { return (m_b.X() - m_a.X()) * (m_b.Y() - m_a.Y()) * (m_b.Z() - m_a.Z()); }
	inline double			Area		(void) const { return 2.0 * ((m_b.X() - m_a.X()) * (m_b.Y() - m_a.Y()) + (m_b.X() - m_a.X()) * (m_b.Z() - m_a.Z()) + (m_b.Y() - m_a.Y()) * (m_b.Y() - m_a.Y())); }
	inline double			Perimeter	(void) const { return 4.0 * ((m_b.X() - m_a.X()) + (m_b.Y() - m_a.Y()) + (m_b.Z() - m_a.Z())); }
	inline bool				Contains	(const Vector3& point) { return !(point.X() < m_a.X() || point.X() > m_b.X() || point.Y() < m_a.Y() || point.Y() > m_b.Y() || point.Z() < m_a.Z() || point.Z() > m_b.Z()); }

private:

	inline void SolveOrder(void)
	{
		if(m_a.X() > m_b.X())
		{
			double tmp = m_a.X(); m_a.SetX(m_b.X()); m_b.SetX(tmp);
		}
		if(m_a.Y() > m_b.Y())
		{
			double tmp = m_a.Y(); m_a.SetY(m_b.Y()); m_b.SetY(tmp);
		}
		if(m_a.Z() > m_b.Z())
		{
			double tmp = m_a.Z(); m_a.SetZ(m_b.Z()); m_b.SetZ(tmp);
		}
	}

	Vector3 m_a, m_b;

};

#endif // AABB_H
