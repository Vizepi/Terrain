#include "Terrain.h"

#include "Perlin.h"
#include "Matrix.h"

#include <cmath>
#include "PI.h"
#include <fstream>
#include <iostream>

#include <QImage>
#include <QVector>

#define VERBOSE(msg) if(m_verbose) { std::cout << msg << std::endl; }

/**
 * @brief Create a terrain given parameters.
 * @param builder A set of octave parameters for perlin generation.
 * @param resolution Dimension of the terrain height matrix.
 * @param aabb Bounding box that contains the terrain.
 * @param seed A number to seed the perlin permutation shuffle.
 */
Terrain::Terrain(const TerrainBuilder& builder, uint64_t resolution, const AABB3& aabb, uint64_t seed, bool verbose)
	: m_resolution(resolution)
	, m_bufferRock(nullptr)
	, m_bufferDirt(nullptr)
	, m_gradient(nullptr)
	, m_aabb(aabb)
	, m_verbose(verbose)
	, m_generator(seed)
{
	VERBOSE("Generating terrain")
	if(resolution > 0)
	{
		uint64_t size = resolution * resolution;
		m_bufferRock = new double[size];
		m_bufferDirt = new double[size];
		m_gradient = new Vector2[size];

		Perlin perlin(std::uniform_int_distribution<uint64_t>()(m_generator));
		double frequency = 0.0, amplitude = 0.0, rotation = 0.0;
		Vector2 offset;

		double height = 0.0;
		for(uint64_t octave = 0; octave < builder.GetOctaveCount(); ++octave)
		{
			builder.GetOctave(octave, frequency, amplitude, offset, rotation);
			height += amplitude;
		}
		for(uint64_t j = 0; j < m_resolution; ++j)
		{
			for(uint64_t i = 0; i < m_resolution; ++i)
			{
				double& h = m_bufferRock[Index(i, j)];
				m_bufferDirt[Index(i, j)] = 0.0;
				h = 0.0;
				for(uint64_t octave = 0; octave < builder.GetOctaveCount(); ++octave)
				{
					builder.GetOctave(octave, frequency, amplitude, offset, rotation);

					h += (1.0 + perlin.Noise(Matrix2x2(M_PI  * rotation / 180.0) * Point2(i, j) * frequency + offset)) * amplitude / 2.0;
				}
				h /= height;
			}
		}
	}
}

/*virtual*/ Terrain::~Terrain(void)
{

}

/**
 * @brief Export the terrain in a .obj file.
 * @param filename Name of the obj file in which to export terrain.
 * @param exportNormals Set to true to export normals in the .obj.
 * @return On success, returns true. If the file cannot be open, returns false.
 */
bool Terrain::ExportOBJ(const std::string& filename, bool exportNormals)
{
	VERBOSE("Exporting OBJ")
	std::ofstream file(filename, std::ios::out);
	file << "o terrain\ng height_map\n";
	if(file.is_open())
	{
		for(uint64_t j = 0; j < m_resolution; ++j)
		{
			for(uint64_t i = 0; i < m_resolution; ++i)
			{
				Vector2 p = Point2(i, j);
				uint64_t index = Index(i, j);
				file << "v " << p.X() << " " << p.Y() << " " <<
						(m_bufferRock[index] + m_bufferDirt[index]) * m_aabb.Size().Z() + m_aabb.A().Z() << "\n";
			}
		}
		for(uint64_t j = 1; j < m_resolution; ++j)
		{
			for(uint64_t i = 1; i < m_resolution; ++i)
			{
				uint64_t a = i + j * m_resolution - m_resolution;
				uint64_t b = i + j * m_resolution + 1 - m_resolution;
				uint64_t c = i + 1 + (j + 1) * m_resolution - m_resolution;
				uint64_t d = i + (j + 1) * m_resolution - m_resolution;
				if(exportNormals)
				{
					file << "f " << a << "//" << a
						 << " " << b << "//" << b
						 << " " << c << "//" << c
						 << " " << d << "//" << d << "\n";
				}
				else
				{
				file << "f " << a
					 << " " << b
					 << " " << c
					 << " " << d << "\n";
				}
			}
		}
		file.close();
		return true;
	}
	return false;
}

/**
 * @brief Export an image of the terrain.
 * Function generate images ***_dirt.*** and ***_rock.***. Each image contains
 * one of the two layers of the terrain. For example, if filename is "Terrain.png",
 * function generate images Terrain_dirt.png and Terrain_rock.png.
 * @param filename Name of the image.
 * @param doublePrecision If set, generate a 16-bits per pixels image, otherwise a 8-bits per pixels image.
 * @return If images can be saved,returns true, otherwise returns false.
 */
bool Terrain::ExportIMG(const std::string& filename, bool doublePrecision)
{
	VERBOSE("Exporting image")
	QImage* rock = nullptr;
	QImage* dirt = nullptr;
	if(doublePrecision)
	{
		rock = new QImage(m_resolution, m_resolution, QImage::Format_RGB16);
		dirt = new QImage(m_resolution, m_resolution, QImage::Format_RGB16);
	}
	else
	{
		rock = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
		dirt = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
		QVector<QRgb> colorTable(256);
		for(uint64_t i = 0; i < 256; ++i)
		{
			colorTable[i] = qRgb(i, i, i);
		}
		rock->setColorTable(colorTable);
		rock->setColorCount(256);
		dirt->setColorTable(colorTable);
		dirt->setColorCount(256);
	}
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			rock->setPixel(i, j, int(fmin(1.0, fmax(0.0, m_bufferRock[Index(i, j)])) * (doublePrecision ? 0x10000 : 0x100)));
			dirt->setPixel(i, j, int(fmin(1.0, fmax(0.0, m_bufferDirt[Index(i, j)])) * (doublePrecision ? 0x10000 : 0x100)));
		}
	}
	size_t dotPosition = filename.find_last_of('.');
	std::string rockname = filename.substr(0, dotPosition) + "_rock" + filename.substr(dotPosition);
	std::string dirtname = filename.substr(0, dotPosition) + "_dirt" + filename.substr(dotPosition);

	bool success = true;
	success = success && rock->save(rockname.c_str());
	success = success && dirt->save(dirtname.c_str());
	delete rock;
	delete dirt;
	return success;
}

/**
 * @brief Get the height of the terrain at a given position.
 * @param position Position where to get the height.
 * @return Height of the terrain at the given position or 0.0 if position is out of the terrain space.
 */
double Terrain::Height(const Vector2& position)
{
	double deltaX = m_aabb.Size().X();
	double deltaY = m_aabb.Size().Y();
	double u = (position.X() - m_aabb.A().X()) / deltaX;
	double v = (position.Y() - m_aabb.A().Y()) / deltaY;

	if(0.0 > u || 1.0 <= u || 0.0 > v || 1.0 <= v)
	{
		return 0.0;
	}

	uint64_t i = u * m_resolution;
	uint64_t j = v * m_resolution;
	double cu = ((position.X() - m_aabb.A().X()) - ((i * deltaX) / (m_resolution - 1))) / (deltaX / (m_resolution - 1));
	double cv = ((position.Y() - m_aabb.A().Y()) - ((j * deltaY) / (m_resolution - 1))) / (deltaY / (m_resolution - 1));
	return (Bilinear(m_bufferRock, cu, cv, i, j) + Bilinear(m_bufferDirt, cu, cv, i, j)) * m_aabb.Size().Z() + m_aabb.A().Z();
}

/**
 * @brief Get the position at a given index of the matrix.
 * @param x The column of the matrix.
 * @param y The line of the matrix.
 * @return A vector containing position at the given index.
 */
Vector2 Terrain::Point2(uint64_t x, uint64_t y)
{
	return Vector2((double(x) * m_aabb.Size().X() / double(m_resolution)) + m_aabb.A().X(), (double(y) * m_aabb.Size().Y()/ double(m_resolution)) + m_aabb.A().Y());
}

/**
 * @brief Get the position of the terrain at a given index of the matrix.
 * @param x The column of the matrix.
 * @param y The line of the matrix.
 * @return A vector containing position and altitude at the given index.
 */
Vector3 Terrain::Point3(uint64_t x, uint64_t y)
{
	Vector2 p = Point2(x, y);
	return Vector3(p.X(), p.Y(), Height(p));
}

/**
 * @brief Compute the bilinear interpolation of a point in a quad.
 * @param buffer The buffer in which to take the quad.
 * @param squarePositionX The position X in the quad.
 * @param squarePositionY The position Y in the quad.
 * @param squareIndexI The index I of the quad.
 * @param squareIndexJ The index J of the quad.
 * @return Bilinear interpolation of (squarePositionX,squarePositionY) in quad (squareIndexI, squareIndexJ).
 */
double Terrain::Bilinear(double* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ)
{
	return	(1 - squarePositionX) * (1 - squarePositionY) * buffer[Index(squareIndexI, squareIndexJ)] +
			squarePositionX *		(1 - squarePositionY) * buffer[Index(squareIndexI + 1, squareIndexJ)] +
			squarePositionX *		squarePositionY *		buffer[Index(squareIndexI + 1, squareIndexJ + 1)] +
			(1 - squarePositionX) * squarePositionY *		buffer[Index(squareIndexI, squareIndexJ + 1)];
}

double Terrain::GetMaxSlope(const Vector2& crtPos, Vector2* nextPos)
{
    //
    // For the eight direction, compute the slope of the position
    double maxSlope = -1000;
    int newX = -1;
    int newY = -1;
    for(uint x = -1; x < 1; ++x)
    {
        for(uint y = -1; y < 1; ++y)
        {
            double slope = abs(Height(Vector2(crtPos.X()+x, crtPos.Y()+y) - crtPos));
            if(slope > maxSlope)
            {
                newX = x + crtPos.X();
                newY = y + crtPos.Y();
                maxSlope = slope;
            }        }
    }
    nextPos->SetX(newX);
    nextPos->SetY(newY);
    return maxSlope;
}

//
// Simulator
//
// TODO Compute with rock speed
//
// TODO Use a better slope system
/**
 * @brief Terrain::Erode
 * @param passCount
 * @param maxSlopeForDirt
 * @param maxDirtLevel
 * @param minDrop
 * @param maxDrop
 * @param stoppingSpeed
 */
void Terrain::Erode (uint64_t passCount, double maxSlopeForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingSlope)
{
	VERBOSE("Eroding terrain")
	//
	// Create random
	std::uniform_int_distribution<uint64_t> rand_u64(0, m_resolution - 1);
	std::uniform_real_distribution<double> rand_dbl(0.0, maxDrop - minDrop);

    //
    // First generation
    Vector2* tmpVec2 = new Vector2();
    for(uint i = 0; i < m_resolution; ++i)
    {
        for(uint j = 0; j < m_resolution; ++j)
        {
            //
            // For the eight direction, compute the slope of the position
			double maxSlope = GetMaxSlope(Vector2(i, j), tmpVec2);

            //
            // Compute the dirt level on this point
            double dirtLevel = maxDirtLevel - ((maxSlope/maxSlopeForDirt) * maxDirtLevel);
            dirtLevel = std::max(0.0, dirtLevel);
            m_bufferDirt[Index(i, j)] = dirtLevel;
        }
    }

    //
    // Simulation loop drop passCount rock
    for(uint64_t nPass = 0; nPass < passCount; nPass++)
    {
        //
        // Choose a random position
		int x = rand_u64(m_generator);
		int y = rand_u64(m_generator);
        //
        // Compute the level of rock falling
        // TODO change ?
		double fallingRock = rand_dbl(m_generator) + minDrop;
        double crtDirt = m_bufferDirt[Index(x, y)];
        fallingRock = std::max(fallingRock, crtDirt);

        //
        // Tear off the rock
        m_bufferDirt[Index(x, y)] = crtDirt - fallingRock;

        //
        // Make the rock fall
        bool stoped = false;
        while(!stoped)
        {
            //
            // Compute the new position
			double slope = GetMaxSlope(Vector2(x, y), tmpVec2);

            //
            //Check the stopping state
            if(slope <= stoppingSlope)
            {
                stoped = true;
            }
            else
            {
                x = tmpVec2->X();
                y = tmpVec2->Y();
            }
        }
    }
}

/**
 * @brief Add ridges on the terrain.
 * @param builder A terrain builder to generate ridge heightfield.
 * @param altitude A vector containing maximal altitude in X and minimal atlitude in Y.
 */
void Terrain::Ridge(const TerrainBuilder& builder, const Vector2& altitude)
{
	VERBOSE("Adding ridge")
	Perlin perlin(std::uniform_int_distribution<uint64_t>()(m_generator));
	double frequency = 0.0, amplitude = 0.0, rotation = 0.0;
	Vector2 offset;

	double height = 0.0;
	for(uint64_t octave = 0; octave < builder.GetOctaveCount(); ++octave)
	{
		builder.GetOctave(octave, frequency, amplitude, offset, rotation);
		height += amplitude;
	}
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			double& hRock = m_bufferRock[Index(i, j)];
			double h = 0.0;
			for(uint64_t octave = 0; octave < builder.GetOctaveCount(); ++octave)
			{
				builder.GetOctave(octave, frequency, amplitude, offset, rotation);

				h += (1.0 + perlin.Noise(Matrix2x2(M_PI  * rotation / 180.0) * Point2(i, j) * frequency + offset)) * amplitude / 2.0;
			}
			h /= height;

			h = (h * (altitude.X() - altitude.Y()) + altitude.Y());

			double hRock2 = hRock * m_aabb.Size().Z() + m_aabb.A().Z();
			double hDirt2 = m_bufferDirt[Index(i, j)] * m_aabb.Size().Z();
			h -= hDirt2;

			if(hRock2 > h)
			{
				hRock = (h + h - hRock2 - m_aabb.A().Z()) / m_aabb.Size().Z();
			}
		}
	}
}
