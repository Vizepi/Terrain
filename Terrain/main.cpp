#include "Terrain.h"

int main(int argc, char *argv[])
{
	TerrainBuilder builder(8);
	builder.SetOctave(0, 1.0/10000.0,	500,	Vector2(0.0, 0.0), 0.0);
	builder.SetOctave(1, 1.0/5000.0,	150,	Vector2(0.0, 0.0), 10.0);
	builder.SetOctave(2, 1.0/1000.0,	75,		Vector2(0.0, 0.0), 25.0);
	builder.SetOctave(3, 1.0/500.0,		25,		Vector2(0.0, 0.0), 20.0);
	builder.SetOctave(4, 1.0/100.0,		5,		Vector2(0.0, 0.0), 30.0);
	builder.SetOctave(5, 1.0/50.0,		3,		Vector2(0.0, 0.0), 50.0);
	builder.SetOctave(6, 1.0/20.0,		1,		Vector2(0.0, 0.0), 40.0);
	builder.SetOctave(7, 1.0/1.0,		0.01,	Vector2(0.0, 0.0), 40.0);

	TerrainBuilder ridge(4);
	ridge.SetOctave(0, 1.0/50000.0,		300,		Vector2(0.0, 0.0), 60.0);
	ridge.SetOctave(1, 1.0/10000.0,		100,		Vector2(0.0, 0.0), 30.0);
	ridge.SetOctave(2, 1.0/1000.0,		50,		Vector2(0.0, 0.0), 10.0);
	ridge.SetOctave(3, 1.0/100.0,		5,		Vector2(0.0, 0.0), 0.0);

	Terrain t(builder, 1000, AABB3(Vector3(-5000, -5000, 0), Vector3(5000, 5000, 1500)), 1, true);
	t.Ridge(ridge, Vector2(1250, 750));
	t.ExportOBJ("Output/Terrain.obj", false);
	t.ExportIMG("Output/Terrain.png", false);

	return 0;
}
