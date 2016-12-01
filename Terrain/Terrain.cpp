#include "Terrain.h"

#include "Perlin.h"
#include "Matrix.h"

#include <cmath>
#include "PI.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <list>
#include <algorithm>

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
	, m_seed(seed)
	, m_generator(seed)
{
	VERBOSE("Generating terrain");

	//
	// Allocate buffers
	uint64_t size = resolution * resolution;
	m_bufferRock = new double[size];
	m_bufferDirt = new double[size];
	m_gradient = new Vector2[size];

	//
	// Create temporaries
	ReSeed();
	Perlin perlin(std::uniform_int_distribution<uint64_t>()(m_generator));
	uint64_t octaveCount = builder.GetOctaveCount();
	double* frequencies = new double[octaveCount];
	double* amplitudes = new double[octaveCount];
	double* rotations = new double[octaveCount];
	Vector2* offsets = new Vector2[octaveCount];

	//
	// Copy builder values to get fast access and total amplitude
	for(uint64_t octave = 0; octave < octaveCount; ++octave)
	{
		builder.GetOctave(octave, frequencies[octave], amplitudes[octave], offsets[octave], rotations[octave]);
	}

	//
	// Generate each height
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			//
			// Initialize each height to 0.0
			double& h = m_bufferRock[Index(i, j)];
			m_bufferDirt[Index(i, j)] = 0.0;
			h = 0.0;

			//
			// Generate each octave and sums them
			for(uint64_t octave = 0; octave < octaveCount; ++octave)
			{
				h += (1.0 + perlin.Noise(Matrix2x2(M_PI  * rotations[octave] / 180.0) * Point2(i, j) * frequencies[octave] + offsets[octave])) * amplitudes[octave] / 2.0;
			}

			//
			// Clamp terrain and put value in range [0;1]
			if(h < m_aabb.A().Z())
			{
				h = m_aabb.A().Z();
			}
			if(h > m_aabb.B().Z())
			{
				h = m_aabb.B().Z();
			}
			h = (h - m_aabb.A().Z()) / m_aabb.Size().Z();
		}
	}

	//
	// Clean memory
	delete[] frequencies;
	delete[] amplitudes;
	delete[] rotations;
	delete[] offsets;
}

/*virtual*/ Terrain::~Terrain(void)
{
	if(m_resolution > 0)
	{
		delete[] m_bufferDirt;
		delete[] m_bufferRock;
		delete[] m_gradient;
	}
}

/**
 * @brief Export the terrain in a .obj file.
 * @param filename Name of the obj file in which to export terrain.
 * @param exportNormals Set to true to export normals in the .obj.
 * @param exportUV Export colors in .obj.
 * @return On success, returns true. If the file cannot be open, returns false.
 */
bool Terrain::ExportOBJ(const std::string& filename, bool exportNormals, bool exportUV)
{
	VERBOSE("Exporting OBJ");

	//
	// Create file
	std::ofstream file(filename, std::ios::out);
	file << "o terrain\ng height_map\n";
	if(file.is_open())
	{
		double maxHeight = std::numeric_limits<double>().min();
		//
		// Save vertices
		VERBOSE("\tVertices...");
		for(uint64_t j = 0; j < m_resolution; ++j)
		{
			for(uint64_t i = 0; i < m_resolution; ++i)
			{
				Vector2 p = Point2(i, j);
				uint64_t index = Index(i, j);

				double height = (m_bufferRock[index] + m_bufferDirt[index]);
				file << "v " << p.X() << " " << p.Y() << " " << height * m_aabb.Size().Z() + m_aabb.A().Z() << "\n";
				if(height > maxHeight)
				{
					maxHeight = height;
				}
			}
		}
		//
		// Save normals
		if(exportNormals)
		{
			VERBOSE("\tNormals...");
			Gradient();
			for(uint64_t j = 0; j < m_resolution; ++j)
			{
				for(uint64_t i = 0; i < m_resolution; ++i)
				{
					uint64_t index = Index(i, j);
					Vector2 gradient = m_gradient[index];
					Vector3 normal(gradient.X(), gradient.Y(), sqrt(gradient.Length() - pow(gradient.X(), 2.0) - pow(gradient.Y(), 2.0)));
					normal.Normalize();
					file << "vn " << normal.X() << " " << normal.Y() << " " << normal.Z() << "\n";
				}
			}
		}
		// Save UVs
		if(exportUV)
		{
			VERBOSE("\tUVs...");
			for(uint64_t i = 0; i < 256; ++i)
			{
				file << "vt 0.5 " << (i + 0.5) / 256.0 << "\n";
			}
		}
		//
		// Save faces
		VERBOSE("\tFaces...");
		char buffer[10];
		for(uint64_t j = 1; j < m_resolution; ++j)
		{
			for(uint64_t i = 1; i < m_resolution; ++i)
			{
				uint64_t a = i + j * m_resolution - m_resolution;
				uint64_t b = i + j * m_resolution + 1 - m_resolution;
				uint64_t c = i + 1 + (j + 1) * m_resolution - m_resolution;
				uint64_t d = i + (j + 1) * m_resolution - m_resolution;
				uint64_t n = j + i * m_resolution - m_resolution;;
				if(exportNormals)
				{
					buffer[0] = 0;
					if(exportUV)
					{
						sprintf(buffer, "%d", int((m_bufferRock[n] + m_bufferDirt[n]) * 255.0 / maxHeight) + 1);
					}
					//
					// Save vertices and normals
					file << "f " << a << "/" << buffer << "/" << a
						 << " " << b << "/" << buffer << "/" << b
						 << " " << c << "/" << buffer << "/" << c
						 << " " << d << "/" << buffer << "/" << d << "\n";
				}
				else
				{
					if(exportUV)
					{
						sprintf(buffer, "%d", int((m_bufferRock[n] + m_bufferDirt[n]) * 255.0 / maxHeight) + 1);
						//
						// Save vertices and UVs
						file << "f " << a << "/" << buffer
							 << " " << b << "/" << buffer
							 << " " << c << "/" << buffer
							 << " " << d << "/" << buffer << "\n";
					}
					else
					{
						//
						// Save vertices
						file << "f " << a
							 << " " << b
							 << " " << c
							 << " " << d << "\n";
					}
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
	VERBOSE("Exporting image");

	//
	// Allocate images
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

	//
	// Copy pixels in the images
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			rock->setPixel(i, j, int(fmin(1.0, fmax(0.0, m_bufferRock[Index(i, j)])) * (doublePrecision ? 0xffff : 0xff)));
			dirt->setPixel(i, j, int(fmin(1.0, fmax(0.0, m_bufferDirt[Index(i, j)])) * (doublePrecision ? 0xffff : 0xff)));
		}
	}

	//
	// Create filenames
	size_t dotPosition = filename.find_last_of('.');
	std::string rockname = filename.substr(0, dotPosition) + "_rock" + filename.substr(dotPosition);
	std::string dirtname = filename.substr(0, dotPosition) + "_dirt" + filename.substr(dotPosition);

	//
	// Save images, clean memory and return
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
 * @brief Get the height of rock of the terrain at a given position.
 * @param position Position where to get the height.
 * @return Height of the terrain at the given position or 0.0 if position is out of the terrain space.
 */
double Terrain::HeightRock(const Vector2& position)
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
	return Bilinear(m_bufferRock, cu, cv, i, j) * m_aabb.Size().Z() + m_aabb.A().Z();
}

/**
 * @brief Get the height of dirt of the terrain at a given position.
 * @param position Position where to get the height.
 * @return Height of the terrain at the given position or 0.0 if position is out of the terrain space.
 */
double Terrain::HeightDirt(const Vector2& position)
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
	return Bilinear(m_bufferDirt, cu, cv, i, j) * m_aabb.Size().Z();
}

/**
 * @brief Get the slope of the terrain at a given position.
 * @param position Position where to get the height.
 * @return Slope from the gradient.
 */
double Terrain::Slope(const Vector2& position)
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
	Vector2 bilinear = Bilinear(m_gradient, cu, cv, i, j);
	return bilinear.Length();
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
 * @brief Import a 2D map into this terrain.
 * @param carveRock Set to true id rock buffer is used for carving.
 * @param carveDirt Set to true if dirt buffer is used for carving.
 * @param rock Rock buffer to carve on the terrain.
 * @param dirt Dirt buffer to carve on the terrain.
 * @param mask Blending mask. A value of 1.0 correspond to 100% carving buffer, and 0.0 to 100% terrain buffer.
 * @param resolution Length of the carved surface square side.
 * @param position The position where to carve the surface, corresponding to the center of the carved surface.
 * @param scale Scaling of the carved surface.
 * @param rotation Rotation angle of the carved surface in degrees.
 */
void Terrain::Carve(bool carveRock, bool carveDirt, double* rock, double* dirt, double* mask, uint64_t resolution, const Vector2& position, double scale, double rotation)
{
	Matrix2x2 rotationMatrix(rotation);
	if(carveRock)
	{
		Vector3* rockBuffer = new Vector3[resolution * resolution];
		for(uint64_t j = 0; j < resolution; ++j)
		{
			for(uint64_t i = 0; i < resolution; ++i)
			{
				Vector2 p((i - (resolution / 2.0)) * scale, (j - (resolution / 2.0)) * scale);
				p = rotationMatrix * p;
				rockBuffer[i + j * resolution] = Vector3(p.X(), p.Y(), rock[i + j * resolution]);
			}
		}

	}
	if(carveDirt)
	{

	}
}

double Terrain::GetMaxSlope(const Vector2& crtPos, Vector2* nextPos)
{
    //
    // Get the gradient at the current position
    /*Vector2 grad = m_gradient[Index(crtPos.X(), crtPos.Y())];
    grad.Normalize();
    double angle = DotProduct(Vector2(1.0, 0.0), grad);

    nextPos->SetX(crtPos.X());
    nextPos->SetY(crtPos.Y());

    double squarelength = m_aabb.Size().X()/(double)m_resolution;
    //
    // Find the next position
    if(angle > M_PI/6.0 && angle < (5*M_PI)/6.0)
    {
        nextPos->SetY(crtPos.Y()+1);
    }
    else if(angle > (7*M_PI)/6.0 && angle < (11*M_PI)/6.0)
    {
        nextPos->SetY(crtPos.Y()-1);
    }

    if(angle < M_PI/3.0 && angle < (5*M_PI)/3.0)
    {
        nextPos->SetX(crtPos.X()+1);
    }
    else if(angle > (2*M_PI)/3.0 && angle < (4*M_PI)/3.0)
    {
        nextPos->SetX(crtPos.X()-1);
    }

    double height = std::abs(Height(crtPos) - Height(*nextPos));*/

    nextPos->SetX(crtPos.X());
    nextPos->SetY(crtPos.Y());
    double squarelength = m_aabb.Size().X()/(double)m_resolution;

    double maxHeight = 0;
    double height = 0;
    for(int crtX = -1; crtX <= 1; ++crtX)
    {
        if(crtPos.X()+crtX > 0 && crtPos.X()+crtX < m_resolution)
        {
            for(int crtY = -1; crtY <= 1; ++crtY)
            {
                if(crtPos.Y()+crtY > 0 && crtPos.Y()+crtY < m_resolution && !(crtX == 0 && crtY ==0))
                {
                    //height = Height(crtPos) - Height(Vector2(crtPos.X()+crtX, crtPos.Y()+crtY));//m_bufferRock[Index(crtPos.X(), crtPos.X())] - m_bufferRock[Index(crtPos.X()+crtX, crtPos.X()+crtY)];
                    double hcrt = m_bufferRock[Index(crtPos.X(), crtPos.Y())] * m_aabb.Size().Z() + m_aabb.A().Z();
                    double hnext = m_bufferRock[Index(crtPos.X()+crtX, crtPos.Y()+crtY)] * m_aabb.Size().Z() + m_aabb.A().Z();
                    height = hcrt-hnext;
                    if(height > maxHeight)
                    {
                        nextPos->SetX(crtPos.X()+crtX);
                        nextPos->SetY(crtPos.Y()+crtY);
                        maxHeight = height;
                    }
                }
            }
        }
    }
    height = maxHeight;

    double length;
    if(nextPos->X() != crtPos.X() && nextPos->Y() != crtPos.Y())
    {
        length = sqrt(squarelength*squarelength*2);
    }
    else
    {
        length = squarelength;
    }
    double angleRes = atan(height/length)*180.0/M_PI;
    return angleRes;
}

double Terrain::GetMaxSlopeWithDirt(const Vector2& crtPos, Vector2* nextPos)
{
    nextPos->SetX(crtPos.X());
    nextPos->SetY(crtPos.Y());
    double squarelength = m_aabb.Size().X()/(double)m_resolution;

    double maxHeight = 0;
    double height = 0;
    for(int crtX = -1; crtX <= 1; ++crtX)
    {
        if(crtPos.X()+crtX > 0 && crtPos.X()+crtX < m_resolution)
        {
            for(int crtY = -1; crtY <= 1; ++crtY)
            {
                if(crtPos.Y()+crtY > 0 && crtPos.Y()+crtY < m_resolution && !(crtX == 0 && crtY ==0))
                {
                    double hcrt = (m_bufferRock[Index(crtPos.X(), crtPos.Y())] + m_bufferDirt[Index(crtPos.X(), crtPos.Y())])
                            * m_aabb.Size().Z() + m_aabb.A().Z();
                    double hnext = (m_bufferRock[Index(crtPos.X()+crtX, crtPos.Y()+crtY)] + m_bufferDirt[Index(crtPos.X()+crtX, crtPos.Y()+crtY)])
                            * m_aabb.Size().Z() + m_aabb.A().Z();
                    height = hcrt-hnext;
                    if(height > maxHeight)
                    {
                        nextPos->SetX(crtPos.X()+crtX);
                        nextPos->SetY(crtPos.Y()+crtY);
                        maxHeight = height;
                    }
                }
            }
        }
    }
    height = maxHeight;

    double length;
    if(nextPos->X() != crtPos.X() && nextPos->Y() != crtPos.Y())
    {
        length = sqrt(squarelength*squarelength*2);
    }
    else
    {
        length = squarelength;
    }
    double angleRes = atan(height/length)*180.0/M_PI;
    return angleRes;
}

//
// Simulator
//
// TODO Compute with rock speed
//
// TODO Use a better slope system
/**
 * @brief Terrain::Erode
 * @param passCount Number of iterations for the ersion process
 * @param maxAngleForDirt the maximum angle in degree where we can put dirt
 * @param maxDirtLevel The maximum level of dirt in one position at the first iteration. In %
 * @param minDrop The minimum dirt Tear off level
 * @param maxDrop The maximum dirt Tear off level
 * @param stoppingAngle The angle at wich the falling rock stop
 */
void Terrain::Erode (uint64_t passCount, uint64_t passWaterCount, double maxAngleForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingAngle, double maxSedimentTransported)
{
    //
    // Get the gradient in each point
    Gradient();
	VERBOSE("Eroding terrain");

	//
	// Create random
	std::uniform_int_distribution<uint64_t> rand_u64(0, m_resolution - 1);
	std::uniform_real_distribution<double> rand_dbl(minDrop/100.0, maxDrop/100.0);
	ReSeed();

    //
    // First generation
    Vector2* tmpVec2 = new Vector2();

    for(uint i = 0; i < m_resolution; ++i)
    {
        for(uint j = 0; j < m_resolution; ++j)
        {
            //
            // For the eight direction, compute the slope of the position
            double angleSlope = GetMaxSlope(Vector2(i, j), tmpVec2);

            //
            // Compute the dirt level on this point
            double dirtFade = 1.0 - std::min(angleSlope/maxAngleForDirt, 1.0);
            double dirtLevel = maxDirtLevel * dirtFade;
            dirtLevel = dirtLevel;
            m_bufferDirt[Index(i, j)] = dirtLevel/100.0;//dirtLevel/100.0;
        }
    }
    VERBOSE("First level of dirt added");

    DirtSmooth();
    VERBOSE("First level of dirt smoothed");

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
        double fallingRock = rand_dbl(m_generator);
        double crtDirt = m_bufferDirt[Index(x, y)];
        //fallingRock = std::max(fallingRock, crtDirt);

        //
        // Tear off the rock
        //m_bufferDirt[Index(x, y)] = std::max(crtDirt - fallingRock, 0.0);
        if(m_bufferDirt[Index(x, y)] == 0)
        {
            m_bufferRock[Index(x, y)] = m_bufferRock[Index(x, y)] - (fallingRock - crtDirt);
        }

        //
        // Make the rock fall
        bool stoped = false;
        while(!stoped)
        {
            //
            // Compute the new position
            double angleSlope = GetMaxSlopeWithDirt(Vector2(x, y), tmpVec2);

            //
            //Check the stopping state
            if(angleSlope <= stoppingAngle)
            {
                stoped = true;
            }
            else
            {
                x = tmpVec2->X();
                y = tmpVec2->Y();
            }
        }
        m_bufferDirt[Index(x, y)] = m_bufferDirt[Index(x, y)] + fallingRock;
    }
    DirtSmooth();

    //
    // Water Erosion
    for(uint64_t nPassWater = 0; nPassWater < passWaterCount; nPassWater++)
    {
        //
        // Choose a random position
        int x = rand_u64(m_generator);
        int y = rand_u64(m_generator);

        //
        // Compute the level of sediment transported
        double fallingSediment = rand_dbl(m_generator);


        /*double crtDirt = m_bufferDirt[Index(x, y)];
        double crtDirtTransported = std::max(std::max(crtDirt, fallingSediment), maxSedimentTransported);*/

        //
        // Tear off the first level of sediment
        m_bufferRock[Index(x, y)] -= fallingSediment;
        m_bufferDirt[Index(x, y)] = fallingSediment;

        //m_bufferRock[Index(x, y)] - crtDirtTransported
        // Make the sediment fall
        bool stoped = false;
        while(!stoped)
        {
            //
            // Compute the new position
            double angleSlope = GetMaxSlopeWithDirt(Vector2(x, y), tmpVec2);

            //
            //Check the stopping state
            if(angleSlope <= stoppingAngle)
            {
                stoped = true;
            }
            else
            {
                x = tmpVec2->X();
                y = tmpVec2->Y();
                double crtDirt = m_bufferDirt[Index(x, y)];

                //
                // Tear off sediment

                //
                // Drop off sediment

                //double NcrtDirtTransported = std::max(crtDirtTransported + std::max(crtDirt, fallingSediment), maxSedimentTransported);
                //m_bufferDirt[Index(x, y)] = m_bufferDirt[Index(x, y)] - (NcrtDirtTransported - crtDirtTransported);
                //crtDirtTransported = NcrtDirtTransported;
                //
                // get sediment depending of the slop.
                //m_bufferDirt[Index(x, y)] = 1.0 - (angleSlope/stoppingAngle) * crtDirtTransported;
            }
        }
        m_bufferDirt[Index(x, y)] = m_bufferDirt[Index(x, y)] + fallingSediment;
    }
    DirtSmooth();
}

void Terrain::DirtSmooth()
{
    int conv[9] = {1, 1, 1,
                   1, 1, 1,
                   1, 1, 1};

    for(uint x = 0; x < m_resolution; ++x)
    {
        for(uint y = 0; y < m_resolution; ++y)
        {
            //
            // Apply smooth on dirt
            double convValue = 0.0;
            int nConv = 0;
            double tConv = 0.0;
            for(int crtX = -1; crtX <= 1; ++crtX)
            {
                for(int crtY = -1; crtY <= 1; ++crtY)
                {
                    if(y+crtY > 0 && y+crtY < m_resolution && x+crtX > 0 && x+crtX < m_resolution)
                    {
                        convValue += conv[nConv] * m_bufferDirt[Index(x+crtX, y+crtY)];
                        tConv++;
                    }
                    nConv++;
                }
            }

            m_bufferDirt[Index(x, y)] = convValue/tConv;
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
	VERBOSE("Adding ridge");

	//
	// Create temporaries
	ReSeed();
	Perlin perlin(std::uniform_int_distribution<uint64_t>()(m_generator));
	uint64_t octaveCount = builder.GetOctaveCount();
	double* frequencies = new double[octaveCount];
	double* amplitudes = new double[octaveCount];
	double* rotations = new double[octaveCount];
	Vector2* offsets = new Vector2[octaveCount];

	//
	// Copy builder values to get fast access and total amplitude
	for(uint64_t octave = 0; octave < octaveCount; ++octave)
	{
		builder.GetOctave(octave, frequencies[octave], amplitudes[octave], offsets[octave], rotations[octave]);
	}

	//
	// Compute ridge
	double minRidge = std::numeric_limits<double>().max();
	double maxRidge = std::numeric_limits<double>().min();
	double* ridge = new double[m_resolution * m_resolution];
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			double h = 0.0;

			//
			// Generate octaves and sums them
			for(uint64_t octave = 0; octave < octaveCount; ++octave)
			{
				h += (1.0 + perlin.Noise(Matrix2x2(M_PI  * rotations[octave] / 180.0) * Point2(i, j) * frequencies[octave] + offsets[octave])) * amplitudes[octave] / 2.0;
			}
			if(h < minRidge)
			{
				minRidge = h;
			}
			if(h > maxRidge)
			{
				maxRidge = h;
			}
			ridge[Index(i, j)] = h;
		}
	}
	double height = maxRidge - minRidge;

	//
	// Apply ridge
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			double& hRock = m_bufferRock[Index(i, j)];
			double h = ridge[Index(i, j)];

			h = (((h - minRidge) / height) * (altitude.X() - altitude.Y()) + altitude.Y());
			//
			// Compare ridge with terrain
			double hRock2 = hRock * m_aabb.Size().Z() + m_aabb.A().Z();
			double hDirt2 = m_bufferDirt[Index(i, j)] * m_aabb.Size().Z();
			h -= hDirt2;

			if(hRock2 > h)
			{
				hRock = (h + h - hRock2 - m_aabb.A().Z()) / m_aabb.Size().Z();
			}
		}
	}

	//
	// Clean memory
	delete[] ridge;
	delete[] frequencies;
	delete[] amplitudes;
	delete[] rotations;
	delete[] offsets;
}

/**
 * @brief Compute gradient on the entire terrain
 */
void Terrain::Gradient(void)
{
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		int64_t jp = int64_t(j)+1, jm = int64_t(j)-1;
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
            double point = m_bufferRock[Index(i, j)] + m_bufferDirt[Index(i, j)];
			int64_t ip = int64_t(i)+1, im = int64_t(i)-1;
            double s = (jm >= 0 ?			m_bufferRock[Index(i, jm)] + m_bufferDirt[Index(i, jm)] : point + point - (m_bufferRock[Index(i, jp)] + m_bufferDirt[Index(i, jp)]));
            double n = (jp < m_resolution ? m_bufferRock[Index(i, jp)] + m_bufferDirt[Index(i, jp)] : point + point - (m_bufferRock[Index(i, jm)] + m_bufferDirt[Index(i, jm)]));
            double w = (im >= 0 ?			m_bufferRock[Index(im, j)] + m_bufferDirt[Index(im, j)] : point + point - (m_bufferRock[Index(ip, j)] + m_bufferDirt[Index(ip, j)]));
            double e = (ip < m_resolution ? m_bufferRock[Index(ip, j)] + m_bufferDirt[Index(ip, j)] : point + point - (m_bufferRock[Index(im, j)] + m_bufferDirt[Index(im, j)]));

			m_gradient[Index(i, j)] = Vector2(e - w, n - s);

		}
	}
}

/**
 * @brief Add influence point on terrain.
 * @param builder The terrain builder that contains influence points definition.
 */
void Terrain::Influence(const TerrainBuilder& builder)
{
	Vector2* ipPositions;
	double* ipRadius;
	uint64_t ipCount = builder.GetInfluencePointCount();
	double minCoef = builder.GetMinPercentHeightCoef();
	if(ipCount)
	{
		ipPositions = new Vector2[ipCount];
		ipRadius = new double[ipCount];
		for(uint64_t nIp = 0; nIp < ipCount; ++nIp)
		{
			builder.GetInfluencePoint(nIp, ipPositions[nIp], ipRadius[nIp]);
		}

		//
		// Rescale height
		for(uint64_t j = 0; j < m_resolution; ++j)
		{
			for(uint64_t i = 0; i < m_resolution; ++i)
			{
				double& h = m_bufferRock[Index(i, j)];

				//
				// Use influence points
				double multiplier = 0.0;
				for(uint64_t nIp = 0; nIp < ipCount; ++nIp)
				{
					double distance = (Point2(i, j) - ipPositions[nIp]).Length();
					if(distance < ipRadius[nIp])
					{
						multiplier += (sin(M_PI * (distance / ipRadius[nIp] + 0.5)) / 2.0) + 0.5;
						multiplier = fmin(multiplier, 1.0);
					}
				}
				if(multiplier  != 0.0)
				{
					h *= multiplier * (1.0 - minCoef) + minCoef;
				}
				else
				{
					h *= minCoef;
				}
			}
		}
	}
}

void Terrain::AddVegetation(const Tree::Builder& builder, uint64_t passCount, uint64_t lockBreak)
{
	VERBOSE("Adding vegetation");
	std::vector<Tree::Instance> fullVegetation;
	std::map<std::string, std::vector<Tree::Instance>> trees;
	ReSeed();
	std::uniform_real_distribution<double> randX(m_aabb.A().X(), m_aabb.B().X() - m_aabb.Size().X() / m_resolution);
	std::uniform_real_distribution<double> randY(m_aabb.A().Y(), m_aabb.B().Y() - m_aabb.Size().Y() / m_resolution);
	std::uniform_real_distribution<double> randA(0.0, 360.0);
	std::uniform_real_distribution<double> randS(builder.survivalMinValue, builder.survivalMaxValue);

	//
	int count;
	double accH = 0.0, accD = 0.0, accS = 0.0;
	double minH = 10000000.0, minD = 100000000.0, minS = 1000000000.0;
	double maxH = -1000000000000.0, maxD = -10000000000000000.0, maxS = -10000000000.0;
	QImage* veget1 = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
	QImage* veget2 = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
	QImage* veget3 = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
	QVector<QRgb> colorTable(256);
	for(uint64_t i = 0; i < 256; ++i)
	{
		colorTable[i] = qRgb(i, 255-i, 0);
	}
	colorTable[0] = qRgb(0, 0, 0);
	veget1->setColorTable(colorTable);
	veget1->setColorCount(256);
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			veget1->setPixel(i, j, 0);
		}
	}
	veget2->setColorTable(colorTable);
	veget2->setColorCount(256);
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			veget2->setPixel(i, j, 0);
		}
	}
	veget3->setColorTable(colorTable);
	veget3->setColorCount(256);
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			veget3->setPixel(i, j, 0);
		}
	}
	//
	VERBOSE("\tGenerating layers...");
	for(uint64_t tree = 0; tree < builder.trees.size(); ++tree)
	{
		trees[builder.trees[tree].GetName()] = std::vector<Tree::Instance>();
		std::vector<Tree::Instance>& treeInstance = trees[builder.trees[tree].GetName()];
		std::uniform_real_distribution<double> randR(builder.trees[tree].radiusMin, builder.trees[tree].radiusMax);

		//
		// Generate tree map
		for(uint64_t pass = 0; pass < passCount; ++pass)
		{
			bool notPlaced = true;
			uint64_t lockBreaker = 0;
			Tree::Instance instance;
			instance.rotation = randA(m_generator);
			instance.name = builder.trees[tree].GetName();

			//
			// Try to place a tree
			while(notPlaced && lockBreaker < lockBreak)
			{
				lockBreaker++;
				notPlaced = false;
				instance.position = Vector2(randX(m_generator), randY(m_generator));
				instance.scale = randR(m_generator);
				double squareRadius = instance.scale * instance.scale;

				//
				// Check if tree collides with other trees
				for(uint64_t nTree = 0; nTree < treeInstance.size() && !notPlaced; ++nTree)
				{
					if((treeInstance[nTree].position - instance.position).SquareLength() <= squareRadius)
					{
						notPlaced = false;
					}
				}
				if(!notPlaced)
				{
					treeInstance.push_back(instance);
				}
			}
		}

		//
		// Check tree survival
		for(int64_t pass = treeInstance.size()-1; pass >= 0; --pass)
		{
			double height = Height(treeInstance[pass].position);
			double dirt = HeightDirt(treeInstance[pass].position);
			double slope = Slope(treeInstance[pass].position);

			double valueHeight = builder.trees[tree].heightCurve.Value(height);
			double valueDirt = builder.trees[tree].dirtCurve.Value(dirt);
			double valueSlope = builder.trees[tree].slopeCurve.Value(slope);
			accH += valueHeight;
			accD += valueDirt;
			accS += valueSlope;
			minH = fmin(minH, valueHeight);
			minD = fmin(minD, valueDirt);
			minS = fmin(minS, valueSlope);
			maxH = fmax(maxH, valueHeight);
			maxD = fmax(maxD, valueDirt);
			maxS = fmax(maxS, valueSlope);
			veget1->setPixel(m_resolution * (treeInstance[pass].position.X() - m_aabb.A().X()) / m_aabb.Size().X(), m_resolution * (treeInstance[pass].position.Y() - m_aabb.A().Y()) / m_aabb.Size().Y(), valueHeight * 255);
			veget2->setPixel(m_resolution * (treeInstance[pass].position.X() - m_aabb.A().X()) / m_aabb.Size().X(), m_resolution * (treeInstance[pass].position.Y() - m_aabb.A().Y()) / m_aabb.Size().Y(), valueDirt * 255);
			veget3->setPixel(m_resolution * (treeInstance[pass].position.X() - m_aabb.A().X()) / m_aabb.Size().X(), m_resolution * (treeInstance[pass].position.Y() - m_aabb.A().Y()) / m_aabb.Size().Y(), valueSlope * 255);
			count++;

			//
			// Check if tree can survive regarding to its environment
			bool canSurvive = true;
			if(builder.useMinAsSurvivalRule)
			{
				canSurvive = 3.0 * randS(m_generator) < (valueHeight + valueDirt + valueSlope);
			}
			else
			{
				canSurvive = randS(m_generator) < (valueHeight * valueDirt * valueSlope);
			}
			if(canSurvive)
			{
				//
				// Insert tree in vegetation
				if(fullVegetation.empty())
				{
					fullVegetation.push_back(treeInstance[pass]);
				}
				else
				{
					fullVegetation.insert(std::upper_bound(fullVegetation.begin(), fullVegetation.end(), treeInstance[pass]), treeInstance[pass]);
				}
			}
			else
			{
				treeInstance.erase(treeInstance.begin() + pass);
			}
		}
	}
	std::cout << accH / count << " " << accD / count << " " << accS / count << "\n";
	std::cout << minH << " " << minD << " " << minS << "\n";
	std::cout << maxH << " " << maxD << " " << maxS << "\n";
	veget1->save("Output/veget1.png");
	delete veget1;
	veget2->save("Output/veget2.png");
	delete veget2;
	veget3->save("Output/veget3.png");
	delete veget3;

	VERBOSE("Merging layers...");

	//
	// While trees are colliding, perform selection
	bool colliding = true;
	while(colliding)
	{
		colliding = false;
		std::list<Tree::Collision> collidingTrees;
		uint64_t id = 0;

		//
		// For ecah tree, compute chances to survive against other trees
		for(std::vector<Tree::Instance>::iterator it = fullVegetation.begin(); it != fullVegetation.end(); ++it)
		{
			Tree::Collision c;
			c.instanceId = id;
			c.count = 0;
			c.probability = 0.0;

			//
			// Compute chance to survive against each tree colliding
			for(std::vector<Tree::Instance>::iterator it2 = it+1; it2 != fullVegetation.end(); ++it2)
			{
				if((it->position - it2->position).Length() < (it->scale + it2->scale))
				{
					colliding = true;
					c.count++;
					double minP = 0.0, maxP = 0.0;
					builder.GetRule(it->name, it2->name, minP, maxP);
					c.probability += std::uniform_real_distribution<double>(minP, maxP)(m_generator);
				}

				//
				// Not necessary to continue if trees are to far
				if(fabs(it->position.X() - it2->position.X()) > (it->scale + it2->scale))
				{
					break;
				}
			}
			if(c.count != 0)
			{
				c.probability /= double(c.count);
				collidingTrees.insert(std::upper_bound(collidingTrees.begin(), collidingTrees.end(), c), c);
			}
			id++;
		}

		uint64_t treesToDeleteCount = builder.perPassTreeSelectionMinimal;
		if(collidingTrees.size() > treesToDeleteCount)
		{
			treesToDeleteCount = collidingTrees.size() * builder.perPassTreeSelectionPercent;
		}
		std::list<uint64_t> ids;

		//
		// Remove last trees
		for(std::list<Tree::Collision>::iterator it = collidingTrees.end(); it != collidingTrees.begin();)
		{
			--it;
			ids.insert(std::upper_bound(ids.begin(), ids.end(), it->instanceId), it->instanceId);
			treesToDeleteCount--;
			if(treesToDeleteCount <= 0)
			{
				break;
			}
		}
		ids.reverse();
		for(std::list<uint64_t>::iterator it = ids.begin(); it != ids.end(); ++it)
		{
			fullVegetation.erase(fullVegetation.begin() + *it);
		}
	}

	m_vegetation.clear();
	for(std::vector<Tree::Instance>::iterator it = fullVegetation.begin(); it != fullVegetation.end(); ++it)
	{
		m_vegetation.push_back(*it);
	}

	//
	// DEBUG
	//
	//
	// Allocate images
	QImage* veget = new QImage(m_resolution, m_resolution, QImage::Format_Indexed8);
	/*QVector<QRgb> colorTable(256);
	for(uint64_t i = 0; i < 256; ++i)
	{
		colorTable[i] = qRgb(i, 255-i, 0);
	}
	colorTable[0] = qRgb(0, 0, 0);*/
	veget->setColorTable(colorTable);
	veget->setColorCount(256);

	//
	// Copy pixels in the images
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			veget->setPixel(i, j, 0);
		}
	}
	//
	// Create output file
	std::ofstream bytes("Output/veget.bytes", std::ios::out | std::ios::binary);
	uint64_t countTree = 2;
	bytes.write((char*)(&countTree), sizeof(uint64_t));
	countTree = m_vegetation.size();
	bytes.write((char*)(&countTree), sizeof(uint64_t));

	for(uint64_t i = 0; i < m_vegetation.size(); ++i)
	{
		uint64_t oak = m_vegetation[i].name == "Oak";
		bytes.write((char*)(&oak), sizeof(uint64_t));
		double tmp = m_vegetation[i].position.X();
		bytes.write((const char*)&tmp, sizeof(double));
		tmp = m_vegetation[i].position.Y();
		bytes.write((char*)&tmp, sizeof(double));
		tmp = Height(m_vegetation[i].position);
		bytes.write((char*)&tmp, sizeof(double));
		bytes.write((char*)&m_vegetation[i].rotation, sizeof(double));
		bytes.write((char*)&m_vegetation[i].scale, sizeof(double));
		veget->setPixel(m_resolution * (m_vegetation[i].position.X() - m_aabb.A().X()) / m_aabb.Size().X(), m_resolution * (m_vegetation[i].position.Y() - m_aabb.A().Y()) / m_aabb.Size().Y(), m_vegetation[i].name == "Pine" ? 255 : 1);
	}
	bytes.close();

	//
	// Save images, clean memory and return
	veget->save("Output/veget.png");
	delete veget;
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

/**
 * @brief Compute the bilinear interpolation of a point in a quad.
 * @param buffer The buffer in which to take the quad.
 * @param squarePositionX The position X in the quad.
 * @param squarePositionY The position Y in the quad.
 * @param squareIndexI The index I of the quad.
 * @param squareIndexJ The index J of the quad.
 * @return Bilinear interpolation of (squarePositionX,squarePositionY) in quad (squareIndexI, squareIndexJ).
 */
Vector2 Terrain::Bilinear(Vector2* buffer, double squarePositionX, double squarePositionY, uint64_t squareIndexI, uint64_t squareIndexJ)
{
	return	buffer[Index(squareIndexI, squareIndexJ)] *			(1 - squarePositionX) * (1 - squarePositionY) +
			buffer[Index(squareIndexI + 1, squareIndexJ)] *		squarePositionX *		(1 - squarePositionY) +
			buffer[Index(squareIndexI + 1, squareIndexJ + 1)] * squarePositionX *		squarePositionY +
			buffer[Index(squareIndexI, squareIndexJ + 1)] *		(1 - squarePositionX) * squarePositionY;
}
