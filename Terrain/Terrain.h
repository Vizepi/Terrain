#ifndef TERRAIN_H
#define TERRAIN_H

#include "TerrainBuilder.h"
#include "AABB.h"

#include <cstdint>
#include <string>

class Terrain
{
public:
					Terrain	(const TerrainBuilder& builder, uint64_t resolution, const AABB3& aabb, uint64_t seed = 0);
	virtual			~Terrain(void);

	bool			ExportOBJ	(const std::string& filename, bool exportNormals = false);
	bool			ExportIMG	(const std::string& filename, bool doublePrecision);

	void			Carve	(double* height, double* mask, uint64_t resolution, const Vector2& position);
	void			Erode	(uint64_t passCount);
	void			Ridge	(const TerrainBuilder& builder, const Vector2& altitude, uint64_t seed = 1);

	inline uint64_t	Index	(uint64_t x, uint64_t y) { return x * m_resolution + y; }
	double			Height	(const Vector2& position);
	Vector2			Point2	(uint64_t x, uint64_t y);
	Vector3			Point3	(uint64_t x, uint64_t y);

private:

	double			Bilinear(double* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ);

	uint64_t	m_resolution;
	double*		m_bufferRock;
	double*		m_bufferDirt;
	Vector2*	m_gradient;
	AABB3		m_aabb;

};

#endif // TERRAIN_H
