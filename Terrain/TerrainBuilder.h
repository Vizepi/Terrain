#ifndef TERRAINBUILDER_H
#define TERRAINBUILDER_H

#include <cstdint>

#include <Vector.h>

class TerrainBuilder
{
public:
	inline			TerrainBuilder	(uint32_t octaveCount = 0);
	inline virtual	~TerrainBuilder	(void);

	inline bool		SetOctave		(const uint32_t octave, const double frequency, const double amplitude, const Vector2& offset, const double rotation);
	inline bool		GetOctave		(const uint32_t octave, double& frequency, double& amplitude, Vector2& offset, double& rotation) const;
	inline uint32_t GetOctaveCount	(void) const;

private:

	uint32_t	m_octaveCount;

	double*		m_frequencies;
	double*		m_amplitudes;
	Vector2*	m_offsets;
	double*		m_rotations;

};

#include "TerrainBuilder.inl"

#endif // TERRAINBUILDER_H
