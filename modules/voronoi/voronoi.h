#ifndef VORONOI_H
#define VORONOI_H

#include "core/reference.h"
#include <string>
#include <map>
#include <vector>


class Voronoi : public Reference
{
	GDCLASS(Voronoi, Reference);

	Dictionary vertexes = Dictionary();
	Array faces = Array();
	Array fragments = Array();

	protected:
		static void _bind_methods();
		void parse_output2d(const std::string &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y);
		void parse_output3d(const std::string &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y, const double &min_z, const double &max_z);
		std::vector<double> generate_random_points2d(int n_points, double min_x, double max_x, double min_y, double max_y);
		std::vector<double> generate_random_points3d(int n_points, const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z);
		bool manual_seed = false;
		double seed = seed;

	public:
		
		void voronoi2d_in_box(PoolVector2Array bounding_points, int n_points);
		void voronoi2d(int n_points);
		void voronoi3d_in_box(PoolVector3Array bounding_points, int n_points);
		void voronoi3d(int n_points);
		Dictionary get_vertexes();
		Array get_faces();
		Array get_fragments();
		void set_seed(double seed);
		double get_seed();

	Voronoi();
};

#endif

