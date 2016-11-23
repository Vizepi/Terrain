#include "Terrain.h"

#include "Perlin.h"
#include "Matrix.h"

#include <cmath>
#include "PI.h"
#include <fstream>
#include <iostream>
#include <limits>

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
	VERBOSE("Generating terrain");

	//
	// Allocate buffers
	uint64_t size = resolution * resolution;
	m_bufferRock = new double[size];
	m_bufferDirt = new double[size];
	m_gradient = new Vector2[size];

	//
	// Create temporaries
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
	double minHeight = std::numeric_limits<double>().max();
	double maxHeight = std::numeric_limits<double>().min();
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
			if(h < minHeight)
			{
				minHeight = h;
			}
			if(h > maxHeight)
			{
				maxHeight = h;
			}
		}
	}

	double height = maxHeight - minHeight;

	//
	// Rescale height
	for(uint64_t j = 0; j < m_resolution; ++j)
	{
		for(uint64_t i = 0; i < m_resolution; ++i)
		{
			double& h = m_bufferRock[Index(i, j)];
			h = (h - minHeight) / height;
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
 * @return On success, returns true. If the file cannot be open, returns false.
 */
bool Terrain::ExportOBJ(const std::string& filename, bool exportNormals)
{
	VERBOSE("Exporting OBJ");

	//
	// Create file
	std::ofstream file(filename, std::ios::out);
	file << "o terrain\ng height_map\n";
	if(file.is_open())
	{
		//
		// Save vertices
		VERBOSE("\tVertices...");
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
		//
		// Save faces
		VERBOSE("\tFaces...");
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
					//
					// Save vertices and normals
					file << "f " << a << "//" << a
						 << " " << b << "//" << b
						 << " " << c << "//" << c
						 << " " << d << "//" << d << "\n";
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
    Vector2 grad = m_gradient[Index(crtPos.X(), crtPos.Y())];
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

    if(nextPos->Y() < 0)
    {
        nextPos->SetY(0);
    }
    else if(nextPos->Y() >= m_resolution)
    {
        nextPos->SetY(m_resolution-1);
    }
    if(nextPos->X() < 0)
    {
        nextPos->SetX(0);
    }
    else if(nextPos->X() >= m_resolution)
    {
        nextPos->SetX(m_resolution-1);
    }

    //
    // Compute the slope angle
    double height = std::max(Height(*nextPos) - Height(crtPos), 0.0);
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
    /*VERBOSE("ANGLE = " << angleRes);
    VERBOSE("height = " << height);
    VERBOSE("length = " << length);*/
    /*VERBOSE("Height(*nextPos) = " << Height(*nextPos));
    VERBOSE("Height(crtPos) = " << Height(crtPos));*/
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
 * @param passCount
 * @param maxSlopeForDirt
 * @param maxDirtLevel
 * @param minDrop
 * @param maxDrop
 * @param stoppingSpeed
 */
void Terrain::Erode (uint64_t passCount, double maxAngleForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingAngle)
{
    //
    // Get the gradient in each point
    Gradient();
	VERBOSE("Eroding terrain");

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
            double angleSlope = GetMaxSlope(Vector2(i, j), tmpVec2);

            //
            // Compute the dirt level on this point
            double dirtFade = 1.0 - std::min(angleSlope/maxAngleForDirt, 1.0);
            double dirtLevel = maxDirtLevel * dirtFade;
            dirtLevel = dirtLevel;
            VERBOSE("dirtLevel = " << dirtLevel)
            m_bufferDirt[Index(i, j)] = dirtLevel;
        }
    }
    VERBOSE("First level of dirt added");

    //
    // Simulation loop drop passCount rock
    /*for(uint64_t nPass = 0; nPass < passCount; nPass++)
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
            double angleSlope = GetMaxSlope(Vector2(x, y), tmpVec2);

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
    }*/
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
			double& h = ridge[Index(i, j)];
			h = 0.0;

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
				//std::cout << h << std::endl;
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
