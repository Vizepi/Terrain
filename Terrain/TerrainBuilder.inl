#include "TerrainBuilder.h"

/**
 * @brief Create a new terrain builder.
 * @param octaveCount Number of octaves to define the terrain.
 */
inline TerrainBuilder::TerrainBuilder(uint32_t octaveCount)
	: m_octaveCount(octaveCount)
	, m_frequencies(octaveCount ? new double[octaveCount] : nullptr)
	, m_amplitudes(octaveCount ? new double[octaveCount] : nullptr)
	, m_offsets(octaveCount ? new Vector2[octaveCount] : nullptr)
	, m_rotations(octaveCount ? new double[octaveCount] : nullptr)
{

}

/**
 * @brief Destroy the terrain builder.
 */
inline /*virtual*/ TerrainBuilder::~TerrainBuilder(void)
{
	if(m_octaveCount)
	{
		delete [] m_frequencies;
		delete [] m_amplitudes;
		delete [] m_offsets;
		delete [] m_rotations;
	}
}

/**
 * @brief Set the values for one octave.
 * @param octave Index octave for which to set values.
 * @param frequency Frequency to set.
 * @param amplitude Amplitude to set.
 * @param offset Offset to set.
 * @param rotation Rotation to set.
 * @return On success, returns true. If octave index is invalid, function returns false.
 */
inline bool TerrainBuilder::SetOctave(const uint32_t octave, const double frequency, const double amplitude, const Vector2& offset, const double rotation)
{
	if(octave < m_octaveCount)
	{
		m_frequencies[octave] = frequency;
		m_amplitudes[octave] = amplitude;
		m_offsets[octave] = offset;
		m_rotations[octave] = rotation;
		return true;
	}
	return false;
}

/**
 * @brief Get the values of one octave.
 * @param octave Index of the octave for which to get values.
 * @param frequency Frequency of the octave.
 * @param amplitude Amplitude of the octave.
 * @param offset Offset of the octave.
 * @param rotation Rotation of the octave.
 * @return On success, returns true. If octave index is invalid, function returns false.
 */
inline bool TerrainBuilder::GetOctave(const uint32_t octave, double& frequency, double& amplitude, Vector2& offset, double& rotation) const
{
	if(octave < m_octaveCount)
	{
		frequency = m_frequencies[octave];
		amplitude = m_amplitudes[octave];
		offset = m_offsets[octave];
		rotation = m_rotations[octave];
		return true;
	}
	return false;
}

/**
 * @brief Get the number of octaves in the builder.
 * @return The octave count.
 */
inline uint32_t TerrainBuilder::GetOctaveCount(void) const
{
	return m_octaveCount;
}
