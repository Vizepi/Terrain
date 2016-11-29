#ifndef TERRAIN_H
#define TERRAIN_H

#include "TerrainBuilder.h"
#include "AABB.h"
#include "Tree.h"

#include <cstdint>
#include <string>
#include <random>
#include <vector>

class Terrain
{
public:
					Terrain				(const TerrainBuilder& builder, uint64_t resolution, const AABB3& aabb, uint64_t seed = 0, bool verbose = false);
	virtual			~Terrain			(void);

	inline void		SetVerbose			(bool verbose) { m_verbose = verbose; }

	bool			ExportOBJ			(const std::string& filename, bool exportNormals = false, bool exportUV = false);
	bool			ExportIMG			(const std::string& filename, bool doublePrecision);

	void			Carve				(bool carveRock, bool carveDirt, double* rock, double* dirt, double* mask, uint64_t resolution, const Vector2& position, double scale = 1.0, double rotation = 0.0);
	double          GetMaxSlope			(const Vector2& crtPos, Vector2* nextPos);
	double          GetMaxSlopeWithDirt	(const Vector2& crtPos, Vector2* nextPos);
	void			Erode				(uint64_t passCount, uint64_t passWaterCount, double maxAngleForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingAngle, double maxSedimentTransported);
	void			Ridge				(const TerrainBuilder& builder, const Vector2& altitude);
	void			Gradient			(void);
	void			Influence			(const TerrainBuilder& builder);
	void			AddVegetation		(const Tree::Builder& builder, uint64_t passCount, uint64_t lockBreak = 0xffffffffffffffff);

	inline uint64_t	Index				(uint64_t x, uint64_t y) { return x * m_resolution + y; }
	double			Height				(const Vector2& position);
	double			HeightRock			(const Vector2& position);
	double			HeightDirt			(const Vector2& position);
	double			Slope				(const Vector2& position);
	Vector2			Point2				(uint64_t x, uint64_t y);
	Vector3			Point3				(uint64_t x, uint64_t y);

private:

	double			Bilinear			(double* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ);
	Vector2			Bilinear			(Vector2* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ);
	inline void		ReSeed				(void) { m_generator.seed(m_seed); }

	uint64_t					m_resolution;
	double*						m_bufferRock;
	double*						m_bufferDirt;
	Vector2*					m_gradient;
	AABB3						m_aabb;
	bool						m_verbose;
	uint64_t					m_seed;
	std::mt19937_64				m_generator;
	std::vector<Tree::Instance> m_vegetation;

};

#endif // TERRAIN_H
