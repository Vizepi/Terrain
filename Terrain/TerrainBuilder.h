#ifndef TERRAINBUILDER_H
#define TERRAINBUILDER_H

#include <cstdint>

#include <Vector.h>

class TerrainBuilder
{
public:
	inline			TerrainBuilder			(uint64_t octaveCount = 0, uint64_t influencePointCount = 0, double minPercentHeightCoef = 0.1);
	inline virtual	~TerrainBuilder			(void);

	inline bool		SetOctave				(const uint64_t octave, const double frequency, const double amplitude, const Vector2& offset, const double rotation);
	inline bool		GetOctave				(const uint64_t octave, double& frequency, double& amplitude, Vector2& offset, double& rotation) const;
	inline bool		SetInfluencePoint		(const uint64_t point, const Vector2& position, const double radius);
	inline bool		GetInfluencePoint		(const uint64_t point, Vector2& position, double& radius) const;
	inline uint64_t GetOctaveCount			(void) const;
	inline uint64_t GetInfluencePointCount	(void) const;
	inline double	GetMinPercentHeightCoef	(void) const;

private:

	uint64_t	m_octaveCount;
	uint64_t	m_influencePointCount;

	double*		m_frequencies;
	double*		m_amplitudes;
	Vector2*	m_offsets;
	double*		m_rotations;

	Vector2*	m_positions;
	double*		m_radius;
	double		m_minPercent;

};

#include "TerrainBuilder.inl"

#endif // TERRAINBUILDER_H
