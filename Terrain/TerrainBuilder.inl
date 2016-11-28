#include "TerrainBuilder.h"

/**
 * @brief Create a new terrain builder.
 * @param octaveCount Number of octaves to define the terrain.
 * @param influencePointCount Number of influence point.
 * @param minPercentHeightCoef The minimal height  multiplier when a point is to far from an influence point.
 */
inline TerrainBuilder::TerrainBuilder(uint64_t octaveCount, uint64_t influencePointCount, double minPercentHeightCoef)
	: m_octaveCount(octaveCount)
	, m_influencePointCount(influencePointCount)
	, m_frequencies(octaveCount ? new double[octaveCount] : nullptr)
	, m_amplitudes(octaveCount ? new double[octaveCount] : nullptr)
	, m_offsets(octaveCount ? new Vector2[octaveCount] : nullptr)
	, m_rotations(octaveCount ? new double[octaveCount] : nullptr)
	, m_positions(influencePointCount ? new Vector2[influencePointCount] : nullptr)
	, m_radius(influencePointCount ? new double[influencePointCount] : nullptr)
	, m_minPercent(minPercentHeightCoef)
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
	if(m_influencePointCount)
	{
		delete [] m_positions;
		delete [] m_radius;
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
inline bool TerrainBuilder::SetOctave(const uint64_t octave, const double frequency, const double amplitude, const Vector2& offset, const double rotation)
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
inline bool TerrainBuilder::GetOctave(const uint64_t octave, double& frequency, double& amplitude, Vector2& offset, double& rotation) const
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
 * @brief Set values for one influence point.
 * @param point The index of the point to set.
 * @param position Position of the point in.
 * @param radius Radius of the point.
 * @return
 */
inline bool TerrainBuilder::SetInfluencePoint(const uint64_t point, const Vector2& position, const double radius)
{
	if(point < m_influencePointCount)
	{
		m_positions[point] = position;
		m_radius[point] = radius;
		return true;
	}
	return false;
}

/**
 * @brief Get values of one influence point.
 * @param point Index of the inluence point to get.
 * @param position Reference filled with the position of the influence point.
 * @param radius Reference filled with radius of the influence point.
 * @return
 */
inline bool TerrainBuilder::GetInfluencePoint(const uint64_t point, Vector2& position, double& radius) const
{
	if(point < m_influencePointCount)
	{
		position = m_positions[point];
		radius = m_radius[point];
		return true;
	}
	return false;
}

/**
 * @brief Get the number of octaves in the builder.
 * @return The octave count.
 */
inline uint64_t TerrainBuilder::GetOctaveCount(void) const
{
	return m_octaveCount;
}

/**
 * @brief Get the number of influence points in the builder.
 * @return  The influence point count.
 */
inline uint64_t TerrainBuilder::GetInfluencePointCount(void) const
{
	return m_influencePointCount;
}

/**
 * @brief Get the min height coeficient for point far from influence points.
 * @return The coeficient.
 */
inline double TerrainBuilder::GetMinPercentHeightCoef(void) const
{
	return m_minPercent;
}
