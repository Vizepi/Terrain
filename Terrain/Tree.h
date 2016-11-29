#ifndef TREE_H
#define TREE_H

#include <Vector.h>
#include <vector>
#include <map>
#include <string>
#include <cstdint>

struct Tree
{
	struct Curve
	{
		inline void AddPoints(const std::vector<Vector2>& points)
		{
			for(uint64_t i = 0; i < points.size(); ++i)
			{
				AddPoint(points[i]);
			}
		}

		inline void AddPoint(const Vector2& point)
		{
			for(std::vector<Vector2>::iterator it = m_points.begin(); it != m_points.end(); ++it)
			{
				if(it->X() == point.X())
				{
					it->SetY(point.Y());
					return;
				}
				if(it->X() < point.X() && (it+1)->X() > point.X())
				{
					m_points.insert(it+1, point);
					return;
				}
			}
			m_points.push_back(point);
		}

		inline double Value(double index) const
		{
			if(m_points.empty())
			{
				return 0.0;
			}
			if(m_points.size() == 0 || m_points.front().X() > index)
			{
				return m_points.front().Y();
			}
			for(uint64_t i = 0; i < m_points.size() - 1; ++i)
			{
				if(m_points[i].X() <= index && m_points[i+1].X() > index)
				{
					double t = (index - m_points[i].X()) / (m_points[i+1].X() - m_points[i].X());
					return m_points[i].Y() * (1.0 - t) + m_points[i+1].Y() * t;
				}
			}
			return m_points.back().Y();
		}
		std::vector<Vector2> m_points;
	};

	inline Tree(const std::string& name) : m_name(name) {}
	inline virtual ~Tree(void) {}
	inline const std::string& GetName(void) const { return m_name; }

	Curve	heightCurve;
	Curve	dirtCurve;
	Curve	slopeCurve;
	double	radiusMin;
	double	radiusMax;

	struct Instance
	{
		Instance(void) {}
		Instance(const Instance& rhs)
			: position(rhs.position)
			, rotation(rhs.rotation)
			, scale(rhs.scale)
			, name(rhs.name)
		{}
		virtual ~Instance(void) {}
		Vector2		position;
		double		rotation;
		double		scale;
		std::string	name;
		inline bool	operator<(const Instance& rhs) const
		{
			return (position - scale) < (rhs.position - scale);
		}
	};

	struct Collision
	{
		uint64_t	instanceId;
		double		probability;
		uint64_t	count;
		inline bool operator<(const Collision& rhs) const
		{
			return probability < rhs.probability;
		}
	};

	struct Superior
	{
		std::string name;
		double		minProbability;
		double		maxProbability;
	};

	struct Rule
	{
		std::string inferior;
		std::map<std::string, Superior> superiors;
	};

	struct Builder
	{
		Builder(void)
			: useMinAsSurvivalRule(true)
			, survivalMinValue(0.0)
			, survivalMaxValue(0.5)
			, perPassTreeSelectionPercent(0.2)
			, perPassTreeSelectionMinimal(1000)
		{

		}
		inline void AddTree(const std::string& name, const std::vector<Vector2>& height, const std::vector<Vector2>& dirt, const std::vector<Vector2>& slope, double radiusMin, double radiusMax)
		{
			trees.push_back(name);
			trees.back().heightCurve.AddPoints(height);
			trees.back().dirtCurve.AddPoints(dirt);
			trees.back().slopeCurve.AddPoints(slope);
			trees.back().radiusMin = radiusMin;
			trees.back().radiusMax = radiusMax;
		}
		inline void AddRule(const std::string& inferior, const std::string& superior, double minProbability, double maxProbability)
		{
			if(minProbability > maxProbability)
			{
				minProbability += maxProbability;
				maxProbability = minProbability - maxProbability;
				minProbability -= maxProbability;
			}
			rules[inferior] = Rule();
			rules[inferior].inferior = inferior;
			rules[inferior].superiors[superior].name = superior;
			rules[inferior].superiors[superior].minProbability = minProbability < 0.0 ? 0.0 : minProbability > 1.0 ? 1.0 : minProbability;
			rules[inferior].superiors[superior].maxProbability = maxProbability < 0.0 ? 0.0 : maxProbability > 1.0 ? 1.0 : maxProbability;
		}
		inline void GetRule(const std::string& inferior, const std::string& superior, double& minProbability, double& maxProbability) const
		{
			std::map<std::string, Rule>::const_iterator i1 = rules.find(inferior);
			std::map<std::string, Rule>::const_iterator i2 = rules.find(superior);
			if(i1 != rules.end())
			{
				std::map<std::string, Superior>::const_iterator it = i1->second.superiors.find(superior);
				if(it != i1->second.superiors.end())
				{
					minProbability = it->second.minProbability;
					maxProbability = it->second.maxProbability;
					return;
				}
			}
			if(i2 != rules.end())
			{
				std::map<std::string, Superior>::const_iterator it = i2->second.superiors.find(inferior);
				if(it != i2->second.superiors.end())
				{
					minProbability = 1.0 - it->second.maxProbability;
					maxProbability = 1.0 - it->second.minProbability;
					return;
				}
			}
			minProbability = 0.0;
			maxProbability = 1.0;
		}

		std::vector<Tree> trees;
		std::map<std::string, Rule> rules;
		bool useMinAsSurvivalRule;
		double survivalMinValue;
		double survivalMaxValue;
		double perPassTreeSelectionPercent;
		uint64_t perPassTreeSelectionMinimal;
	};

private:

	std::string m_name;

};

#endif // TREE_H
