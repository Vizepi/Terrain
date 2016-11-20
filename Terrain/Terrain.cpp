#include "Terrain.h"

#include "Perlin.h"
#include "Matrix.h"

#include <cmath>
#include "PI.h"
#include <fstream>

#include <QImage>
#include <QVector>

Terrain::Terrain(const TerrainBuilder& builder, uint64_t resolution, const AABB3& aabb, uint64_t seed)
	: m_resolution(resolution)
	, m_bufferRock(nullptr)
	, m_bufferDirt(nullptr)
	, m_gradient(nullptr)
	, m_aabb(aabb)
{
	if(resolution > 0)
	{
		uint64_t size = resolution * resolution;
		m_bufferRock = new double[size];
		m_bufferDirt = new double[size];
		m_gradient = new Vector2[size];

		Perlin perlin(seed);
		double frequency = 0.0, amplitude = 0.0, rotation = 0.0;
		Vector2 offset;

		double height = 0.0;
		for(uint64_t octave = 0; octave < builder.GetOctaveCount(); ++octave)
		{
			builder.GetOctave(octave, frequency, amplitude, offset, rotation);
			height += amplitude;
		}

		for(uint64_t j = 0; j < resolution; ++j)
		{
			for(uint64_t i = 0; i < resolution; ++i)
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

bool Terrain::ExportOBJ(const std::string& filename, bool exportNormals)
{
	std::ofstream file(filename, std::ios::out);
	if(file.is_open())
	{
		for(uint64_t j = 0; j < m_resolution; ++j)
		{
			for(uint64_t i = 0; i < m_resolution; ++i)
			{
				Vector3 p = Point3(i, j);
				file << "v " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
			}
		}
		for(uint64_t j = 0; j < m_resolution - 1; ++j)
		{
			for(uint64_t i = 0; i < m_resolution - 1; ++i)
			{
				uint64_t a = i + j * m_resolution + 1;
				uint64_t b = a + m_resolution;
				uint64_t c = b + 1;
				uint64_t d = a + 1;
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

bool Terrain::ExportIMG(const std::string& filename, bool doublePrecision)
{
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
			rock->setPixel(i, j, int(m_bufferRock[Index(i, j)] * (doublePrecision ? 0x10000 : 0x100)));
			dirt->setPixel(i, j, int(m_bufferDirt[Index(i, j)] * (doublePrecision ? 0x10000 : 0x100)));
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

#include <iostream>
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
	double cv = ((position.Y() - m_aabb.A().Y()) - ((i * deltaY) / (m_resolution - 1))) / (deltaY / (m_resolution - 1));
	return (Bilinear(m_bufferRock, cu, cv, i, j) + Bilinear(m_bufferDirt, cu, cv, i, j)) * m_aabb.Size().Z() + m_aabb.A().Z();
}

Vector2 Terrain::Point2(uint64_t x, uint64_t y)
{
	return Vector2((double(x) * m_aabb.Size().X() / double(m_resolution)) + m_aabb.A().X(), (double(y) * m_aabb.Size().Y()/ double(m_resolution)) + m_aabb.A().Y());
}

Vector3 Terrain::Point3(uint64_t x, uint64_t y)
{
	Vector2 p = Point2(x, y);
	return Vector3(p.X(), p.Y(), Height(p));
}

double Terrain::Bilinear(double* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ)
{
	return	(1 - squarePositionX) * (1 - squarePositionY) * buffer[Index(squareIndexI, squareIndexJ)] +
			squarePositionX *		(1 - squarePositionY) * buffer[Index(squareIndexI + 1, squareIndexJ)] +
			squarePositionX *		squarePositionY *		buffer[Index(squareIndexI + 1, squareIndexJ + 1)] +
			(1 - squarePositionX) * squarePositionY *		buffer[Index(squareIndexI, squareIndexJ + 1)];
}

//
// Simulator
void Terrain::Erode (uint64_t passCount, double maxSlopeForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingSpeed)
{
    //
    // First generation
    for(uint i = 0; i < m_resolution; ++i)
    {
        for(uint j = 0; j < m_resolution; ++j)
        {
            //
            // For the eight direction, compute the slope of the position
            double maxSlope = -1000;
            for(uint x = -1; x < 1; ++x)
            {
                for(uint y = -1; y < 1; ++y)
                {
                    double slope = abs(Height(Vector2(i+x, j+y) - Vector2(i, j)));
                    maxSlope = std::max(slope, maxSlope);
                }
            }

            //
            // Compute the dirt level on this point
            double dirtLevel = maxDirtLevel - ((maxSlope/maxSlopeForDirt) * maxDirtLevel);
            dirtLevel = std::max(0.0, dirtLevel);
            m_bufferDirt[Index(i, j)] = dirtLevel;
        }
    }

    //
    // Simulation loop
    for(int nPass = 0; nPass < passCount; nPass++)
    {
        //
        // Drop a new rock

        //
        // Choose a random position
        int x = std::rand() % m_resolution;
        int y = std::rand() % m_resolution;
        //
        // Compute the level of rock falling
        // TODO change
        double fallingRock = (double)(std::rand() % ((int)maxDrop - (int)minDrop)) + minDrop;

        double speed = 1.0f;
        bool stoped = false;
        while(!stoped)
        {
            //
            // Compute the new position
            double maxSlope = -1000;
            int newX = -1;
            int newY = -1;
            for(uint crtX = -1; x < 1; ++crtX)
            {
                for(uint crtY = -1; y < 1; ++crtY)
                {
                    double slope = abs(Height(Vector2(x+crtX, y+crtY) - Vector2(x, y)));
                    if(slope > maxSlope)
                    {
                        newX = x + crtX;
                        newY = y + crtY;
                    }
                }
            }
            //
            // TODO Compute the new speed

            //
            //Check the stopping state
            if(speed <= stoppingSpeed)
            {
                stoped = true;
            }
            else
            {
                x = newX;
                y = newY;
            }
        }
    }
}
