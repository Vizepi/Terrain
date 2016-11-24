#ifndef TERRAIN_H
#define TERRAIN_H

#include "TerrainBuilder.h"
#include "AABB.h"

#include <cstdint>
#include <string>
#include <random>

class Terrain
{
public:
					Terrain	(const TerrainBuilder& builder, uint64_t resolution, const AABB3& aabb, uint64_t seed = 0, bool verbose = false);
	virtual			~Terrain(void);

	inline void		SetVerbose	(bool verbose) { m_verbose = verbose; }

	bool			ExportOBJ	(const std::string& filename, bool exportNormals = false);
	bool			ExportIMG	(const std::string& filename, bool doublePrecision);

	void			Carve       (bool carveRock, bool carveDirt, double* rock, double* dirt, double* mask, uint64_t resolution, const Vector2& position, double scale = 1.0, double rotation = 0.0);
	double          GetMaxSlope (const Vector2& crtPos, Vector2* nextPos);
    void			Erode       (uint64_t passCount, double maxSlopeForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingSlope);
	void			Ridge       (const TerrainBuilder& builder, const Vector2& altitude);
	void			Gradient	(void);
	void			Influence	(const TerrainBuilder& builder);

	inline uint64_t	Index	(uint64_t x, uint64_t y) { return x * m_resolution + y; }
	double			Height	(const Vector2& position);
	Vector2			Point2	(uint64_t x, uint64_t y);
	Vector3			Point3	(uint64_t x, uint64_t y);

private:

	double			Bilinear(double* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ);

	uint64_t		m_resolution;
	double*			m_bufferRock;
	double*			m_bufferDirt;
	Vector2*		m_gradient;
	AABB3			m_aabb;
	bool			m_verbose;
	std::mt19937_64 m_generator;

};

#endif // TERRAIN_H
