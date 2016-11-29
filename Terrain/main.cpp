#include "Terrain.h"

int main(int argc, char *argv[])
{
	TerrainBuilder builder(8, 4, 0.4);
	builder.SetOctave(0, 1.0/5000.0,	1500,	Vector2(0.0, 0.0), 0.0);
	builder.SetOctave(1, 1.0/1000.0,	750,	Vector2(0.0, 0.0), 10.0);
	builder.SetOctave(2, 1.0/500.0,		50,		Vector2(0.0, 0.0), 25.0);
	builder.SetOctave(3, 1.0/100.0,		10,		Vector2(0.0, 0.0), 20.0);
	builder.SetOctave(4, 1.0/50.0,		5,		Vector2(0.0, 0.0), 30.0);
	builder.SetOctave(5, 1.0/20.0,		2,		Vector2(0.0, 0.0), 50.0);
	builder.SetOctave(6, 1.0/10.0,		1,		Vector2(0.0, 0.0), 40.0);
	builder.SetOctave(7, 1.0/1.0,		0.05,	Vector2(0.0, 0.0), 40.0);

	builder.SetInfluencePoint(0, Vector2(-3000, -2000), 3000);
	builder.SetInfluencePoint(1, Vector2(-2500, 2000), 2000);
	builder.SetInfluencePoint(2, Vector2(2000, 4000), 2500);
	builder.SetInfluencePoint(3, Vector2(1000, -1000), 1500);

	TerrainBuilder ridge(4);
	ridge.SetOctave(0, 1.0/5000.0,		2000,	Vector2(0.0, 0.0), 0.0);
	ridge.SetOctave(1, 1.0/1000.0,		300,	Vector2(0.0, 0.0), 10.0);
	ridge.SetOctave(2, 1.0/100.0,		50,		Vector2(0.0, 0.0), 20.0);
	ridge.SetOctave(3, 1.0/10.0,		10,		Vector2(0.0, 0.0), 40.0);


	Terrain t(builder, 250, AABB3(Vector3(-5000, -5000, 1000), Vector3(5000, 5000, 3000)), 666, true);
    t.Ridge(ridge, Vector2(2500, 2250));
    t.Gradient();

    //Erode(uint64_t passCOunt, double maxSlopeForDirt, double maxDirtLevelInPourcentage, double minDrop, double maxDrop, double stoppingSlope)
    t.Erode(0, 0, 35, 1.0, 0.0, 1.0, 10, 1.0);

	//t.Ridge(ridge, Vector2(1800, 1400));
	t.Influence(builder);
	t.Gradient();

	//
	// Adding vegetation
	//
	std::vector<Vector2> pineHeightCurve;
		pineHeightCurve.push_back(Vector2(250.0, 0.0));
		pineHeightCurve.push_back(Vector2(500.0, 1.0));
		pineHeightCurve.push_back(Vector2(2000.0, 1.0));
		pineHeightCurve.push_back(Vector2(3000.0, 0.0));
	std::vector<Vector2> pineDirtCurve;
		pineDirtCurve.push_back(Vector2(5.0, 0.0));
		pineDirtCurve.push_back(Vector2(15.0, 1.0));
	std::vector<Vector2> pineSlopeCurve;
		pineSlopeCurve.push_back(Vector2(0.0, 1.0));
		pineSlopeCurve.push_back(Vector2(0.1, 1.0));
		pineSlopeCurve.push_back(Vector2(0.2, 0.0));

	std::vector<Vector2> oakHeightCurve;
		oakHeightCurve.push_back(Vector2(500.0, 0.0));
		oakHeightCurve.push_back(Vector2(1250.0, 1.0));
		oakHeightCurve.push_back(Vector2(2500.0, 0.0));
	std::vector<Vector2> oakDirtCurve;
		oakDirtCurve.push_back(Vector2(10.0, 0.0));
		oakDirtCurve.push_back(Vector2(50.0, 1.0));
	std::vector<Vector2> oakSlopeCurve;
		oakSlopeCurve.push_back(Vector2(0.0, 1.0));
		oakSlopeCurve.push_back(Vector2(0.05, 0.75));
		oakSlopeCurve.push_back(Vector2(0.1, 0.0));

	Tree::Builder treeBuilder;
	treeBuilder.AddTree("Pine", pineHeightCurve, pineDirtCurve, pineSlopeCurve, 1.0, 3.0);
	treeBuilder.AddTree("Oak", oakHeightCurve, oakDirtCurve, oakSlopeCurve, 1.0, 5.0);
	treeBuilder.AddRule("Pine", "Oak", 0.0, 0.5);
	//treeBuilder.useMinAsSurvivalRule = false;

	t.AddVegetation(treeBuilder, 10000, 1000);

    //Erode(uint64_t assCOunt, double maxSlopeForDirt, double maxDirtLevel, double minDrop, double maxDrop, double stoppingSlope)
	//t.Erode(10000, 35, 10, 5, 10, 35);

	//t.ExportOBJ("Output/Terrain.obj", false);
	t.ExportOBJ("Output/TerrainN.obj", true);
	t.ExportIMG("Output/Terrain.png", false);

	return 0;
}
