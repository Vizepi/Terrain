#ifndef TREE_H
#define TREE_H

#include <Vector.h>
#include <vector>
#include <string>

struct Tree
{
	inline Tree(const std::string& name) : m_name(name) {}

	std::vector<Vector2> heightCurve;
	std::vector<Vector2> humidityCurve;
	std::vector<Vector2> slopeCurve;
	std::vector<Vector2> lightningCurve;

	struct Instance
	{
		Vector2		position;
		double		rotation;
		double		scale;
		std::string	name;
	};

private:

	std::string m_name;

};

#endif // TREE_H
