#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/RboxPoints.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "voronoi.h"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

std::vector<double> Voronoi::generate_random_points2d(int n_points, double min_x, double max_x, double min_y, double max_y) {
	std::vector<double> random_points(n_points * 2);
	std::mt19937 eng;
	if (manual_seed)
	{
		eng = std::mt19937(seed);
	}
	else
	{
		std::random_device rd;
		eng = std::mt19937(rd);
	}
	
	std::uniform_real_distribution<> distr(0.0, 1.0);
	double rand_value = 0;
	for (int i = 0; i < n_points * 2; i += 2) {
		rand_value = distr(eng);
		random_points[i] = rand_value * (max_x - min_x) + min_x;
		rand_value = distr(eng);
		random_points[i + 1] = rand_value * (max_y - min_y) + min_y;
	}
	return random_points;
}

std::vector<double> Voronoi::generate_random_points3d(int n_points, const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z) {
	std::vector<double> random_points(n_points * 3);
	std::mt19937 eng;
	if (manual_seed) {
		eng = std::mt19937(seed);
	} else {
		std::random_device rd;
		eng = std::mt19937(rd);
	}
	std::uniform_real_distribution<> distr(0.0, 1.0);
	double rand_value = 0;
	for (int i = 0; i < n_points * 3; i += 3) {
		rand_value = distr(eng);
		random_points[i] = rand_value * (max_x - min_x) + min_x;
		rand_value = distr(eng);
		random_points[i + 1] = rand_value * (max_y - min_y) + min_y;
		rand_value = distr(eng);
		random_points[i + 2] = rand_value * (max_z - min_z) + min_z;
	}

	return random_points;
}

void Voronoi::voronoi2d(int n_points) {
	PoolVector2Array bounding_points = PoolVector2Array();
	bounding_points.append(Vector2(0, 0));
	bounding_points.append(Vector2(1, 0));
	bounding_points.append(Vector2(1, 1));
	bounding_points.append(Vector2(0, 1));
	Voronoi::voronoi2d_in_box(bounding_points, n_points);
}

void Voronoi::voronoi2d_in_box(PoolVector2Array bounding_points, int n_points) {
	int size = bounding_points.size();

	//get min/max values of bounding box
	double max_x = bounding_points[0].x;
	double max_y = bounding_points[0].y;
	double min_x = bounding_points[0].x;
	double min_y = bounding_points[0].y;
	for (int i = 0; i < size; ++i) {
		if (bounding_points[i].x > max_x) {
			max_x = bounding_points[i].x;
		} else if (bounding_points[i].x < min_x) {
			min_x = bounding_points[i].x;
		}
		if (bounding_points[i].y > max_y) {
			max_y = bounding_points[i].y;
		} else if (bounding_points[i].y < min_y) {
			min_y = bounding_points[i].y;
		}
	}

	//create some random points inside the bounding box
	std::vector<double> points = generate_random_points2d(n_points, min_x, max_x, min_y, max_y);
	//assign the random points to the vector
	std::vector<double> all_points(n_points * 2 * 5);
	for (int i = 0; i < n_points; ++i) {
		all_points[i * 2] = points[i];
		all_points[i * 2 + 1] = points[i + 1];
	}

	//adjust/mirror the points
	//TODO put it in an elegant loop
	for (int i = 0; i < n_points; ++i) {
		//left mirror
		all_points[2 * n_points * 1 + i * 2] = 2 * min_x - all_points[i * 2];
		all_points[2 * n_points * 1 + i * 2 + 1] = all_points[i * 2 + 1];
		//right mirror
		all_points[2 * n_points * 2 + i * 2] = 2 * max_x - all_points[i * 2];
		all_points[2 * n_points * 2 + i * 2 + 1] = all_points[i * 2 + 1];
		//top mirror
		all_points[2 * n_points * 3 + i * 2] = all_points[i * 2];
		all_points[2 * n_points * 3 + i * 2 + 1] = 2 * min_y - all_points[i * 2 + 1];
		//bottom mirror
		all_points[2 * n_points * 4 + i * 2] = all_points[i * 2];
		all_points[2 * n_points * 4 + i * 2 + 1] = 2 * max_y - all_points[i * 2 + 1];
	}

	orgQhull::RboxPoints rbox;

	std::stringstream stringStream;

	stringStream << "2 " << 5 * n_points;
	for (size_t i = 0; i < all_points.size(); ++i) {
		stringStream << ' ' << all_points[i];
	}
	std::istringstream istringStream(stringStream.str());

	rbox.appendPoints(istringStream);
	orgQhull::Qhull qhull;
	std::stringstream output;
	qhull.setOutputStream(&output);
	qhull.runQhull(rbox, "v Qbb");
	qhull.outputQhull("p");
	qhull.outputQhull("FN");
	std::string qhull_output = output.str();

	Voronoi::parse_output2d(qhull_output, min_x, max_x, min_y, max_y);
}

void Voronoi::voronoi3d(int n_points) {
	PoolVector3Array bounding_points = PoolVector3Array();
	bounding_points.append(Vector3(0, 0, 0));
	bounding_points.append(Vector3(1, 0, 0));
	bounding_points.append(Vector3(1, 1, 0));
	bounding_points.append(Vector3(0, 1, 0));
	bounding_points.append(Vector3(0, 1, 1));
	bounding_points.append(Vector3(1, 1, 1));
	bounding_points.append(Vector3(1, 0, 1));
	bounding_points.append(Vector3(0, 0, 1));
	Voronoi::voronoi3d_in_box(bounding_points, n_points);
}

void Voronoi::voronoi3d_in_box(PoolVector3Array bounding_points, int n_points) {

	int dim = 3;
	int size = bounding_points.size();

	//get min/max values of bounding box
	double max_x = bounding_points[0].x;
	double max_y = bounding_points[0].y;
	double max_z = bounding_points[0].z;
	double min_x = bounding_points[0].x;
	double min_y = bounding_points[0].y;
	double min_z = bounding_points[0].z;

	for (int i = 0; i < size; ++i) {
		if (bounding_points[i].x > max_x) {
			max_x = bounding_points[i].x;
		} else if (bounding_points[i].x < min_x) {
			min_x = bounding_points[i].x;
		}
		if (bounding_points[i].y > max_y) {
			max_y = bounding_points[i].y;
		} else if (bounding_points[i].y < min_y) {
			min_y = bounding_points[i].y;
		}
		if (bounding_points[i].z > max_z) {
			max_z = bounding_points[i].z;
		} else if (bounding_points[i].z < min_z) {
			min_z = bounding_points[i].z;
		}
	}

	//create some random points inside the bounding box
	std::vector<double> points = generate_random_points3d(n_points, min_x, max_x, min_y, max_y, min_z, max_z);
	//assign the random points to the vector
	std::vector<double> all_points(n_points * 3 * 7);
	for (int i = 0; i < n_points; ++i) {
		for (int d = 0; d < dim; ++d) {
			all_points[i * dim + d] = points[i + d];
		}
	}


	//adjust/mirror the points
	//TODO put it in an elegant loop
	for (int i = 0; i < n_points; ++i) {
		//left mirror
		all_points[3 * n_points * 1 + i * 3] = 2 * min_x - all_points[i * 3];
		all_points[3 * n_points * 1 + i * 3 + 1] = all_points[i * 3 + 1];
		all_points[3 * n_points * 1 + i * 3 + 2] = all_points[i * 3 + 2];
		//right mirror
		all_points[3 * n_points * 2 + i * 3] = 2 * max_x - all_points[i * 3];
		all_points[3 * n_points * 2 + i * 3 + 1] = all_points[i * 3 + 1];
		all_points[3 * n_points * 2 + i * 3 + 2] = all_points[i * 3 + 2];
		//top mirror
		all_points[3 * n_points * 3 + i * 3] = all_points[i * 3];
		all_points[3 * n_points * 3 + i * 3 + 1] = 2 * min_y - all_points[i * 3 + 1];
		all_points[3 * n_points * 3 + i * 3 + 2] = all_points[i * 3 + 2];
		//bottom mirror
		all_points[3 * n_points * 4 + i * 3] = all_points[i * 3];
		all_points[3 * n_points * 4 + i * 3 + 1] = 2 * max_y - all_points[i * 3 + 1];
		all_points[3 * n_points * 4 + i * 3 + 2] = all_points[i * 3 + 2];
		//up mirror
		all_points[3 * n_points * 5 + i * 3] = all_points[i * 3];
		all_points[3 * n_points * 5 + i * 3 + 1] = all_points[i * 3 + 1];
		all_points[3 * n_points * 5 + i * 3 + 2] = 2 * min_z - all_points[i * 3 + 2];
		//down mirror
		all_points[3 * n_points * 6 + i * 3] = all_points[i * 3];
		all_points[3 * n_points * 6 + i * 3 + 1] = all_points[i * 3 + 1];
		all_points[3 * n_points * 6 + i * 3 + 2] = 2 * max_z - all_points[i * 3 + 2];
	}

	orgQhull::RboxPoints rbox;

	std::stringstream stringStream;

	stringStream << "3 " << 7 * n_points;
	for (size_t i = 0; i < all_points.size(); ++i) {
		stringStream << ' ' << all_points[i];
	}
	std::istringstream istringStream(stringStream.str());

	rbox.appendPoints(istringStream);
	orgQhull::Qhull qhull;
	std::stringstream output;
	qhull.setOutputStream(&output);
	qhull.runQhull(rbox, "v Qbb");
	qhull.outputQhull("o");
	std::string qhull_output = output.str();

	Voronoi::parse_output3d(qhull_output, min_x, max_x, min_y, max_y, min_z, max_z);
}

void Voronoi::parse_output2d(const std::string &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y) {

	std::vector<bool> good_vertexes;
	std::vector<std::vector<int> > vectfaces;

	double epsilon = 0.000001;
	int n_points = 0;
	int n_faces = 0;

	std::string buff{ "" };

	int line_counter = 0;

	bool good_vertex;
	for (auto n : output) {
		if (n != '\n') {
			buff += n;
		} else {
			if (line_counter == 1) {
				n_points = atoi(buff.c_str());
				good_vertexes.reserve(n_points);
			} else if (line_counter == n_points + 2) {
				n_faces = atoi(buff.c_str());
				vectfaces = std::vector<std::vector<int> >(n_faces);
			} else if (line_counter > 1 && line_counter < n_points + 2) {
				std::string::size_type sz;
				double x = std::stod(buff, &sz);
				double y = std::stod(buff.substr(sz));
				good_vertex = (x >= min_x - epsilon && x <= max_x + epsilon && y >= min_y - epsilon && y <= max_y + epsilon);
				good_vertexes.push_back(good_vertex);
				if (good_vertex) {
					vertexes[line_counter - 2] = Vector2(x, y);
				}
			} else if (line_counter > (n_points + 1)) {
				std::string::size_type sz;
				std::string::size_type pos = 0;
				int vect_index = line_counter - 3 - n_points;
				while (pos < buff.size() - 1) {
					int value = std::stoi(buff.substr(pos), &sz);

					if (pos == 0) {
						vectfaces[vect_index] = std::vector<int>();
					} else {
						vectfaces[vect_index].push_back(value);
					}
					pos += sz;
				}
			}
			line_counter++;
			buff.clear();
		}
	}

	faces = Array();

	for (int i = 0; i < vectfaces.size(); ++i) {
		bool good_face = True;
		for (int j = 0; j < vectfaces[i].size(); ++j) {
			if (!good_vertexes[std::abs(vectfaces[i][j])]) {
				good_face = False;
				break;
			}
		}
		if (good_face) {
			Array cur_face = Array();
			for (int j = 0; j < vectfaces[i].size(); ++j) {
				cur_face.append(vectfaces[i][j]);
			}
			faces.append(cur_face);
		}
	}
}

void Voronoi::parse_output3d(const std::string &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y, const double &min_z, const double &max_z) {

	std::vector<bool> good_vertexes;
	std::vector<std::vector<int> > vectfaces;

	double epsilon = 0.000001;
	int n_points = 0;
	int n_faces = 0;

	std::string buff{ "" };

	int line_counter = 0;

	bool good_vertex;

	std::vector<std::vector<double> > voro_vertexes = std::vector<std::vector<double> >();

	for (auto n : output) {
		if (n != '\n') {
			buff += n;
		} else {
			std::string::size_type sz;
			std::string::size_type pos = 0;
			if (line_counter == 1) {
				n_points = std::stoi(buff, &sz);
				pos += sz;
				n_faces = std::stoi(buff.substr(pos));
				vectfaces = std::vector<std::vector<int> >(n_faces);
				good_vertexes.reserve(n_points);
			} else if (line_counter > 1 && line_counter < n_points + 2) {
				double x = std::stod(buff, &sz);
				pos += sz;
				double y = std::stod(buff.substr(pos), &sz);
				pos += sz;
				double z = std::stod(buff.substr(pos), &sz);
				good_vertex = (x >= min_x - epsilon && x <= max_x + epsilon && y >= min_y - epsilon && y <= max_y + epsilon && z >= min_z - epsilon && z <= max_z + epsilon);
				good_vertexes.push_back(good_vertex);
				if (good_vertex) {
					vertexes[line_counter - 2] = Vector3(x, y, z);
				}
				voro_vertexes.push_back(std::vector<double>{ x, y, z });

			} else if (line_counter >= (n_points + 2)) {
				int counter = -1;
				int vect_index = line_counter - 2 - n_points;
				while (pos < buff.size() - 1) {
					int value = std::stoi(buff.substr(pos), &sz);
					if (pos == 0) {
						vectfaces[vect_index] = std::vector<int>();
					} else {
						vectfaces[vect_index].push_back(value);
					}
					pos += sz;
				}
				counter++;
			}
			line_counter++;
			buff.clear();
		}
	}

	faces = Array();
	fragments = Array();
	int fragment_counter = -1;

	for (int i = 0; i < vectfaces.size(); ++i) {
		bool good_face = True;
		for (int j = 0; j < vectfaces[i].size(); ++j) {
			if (!good_vertexes[std::abs(vectfaces[i][j])]) {
				good_face = False;
				break;
			}
		}
		if (good_face && vectfaces[i].size() > 0) {
			fragment_counter++;
			Array cur_face = Array();

			std::stringstream output_hull;
			orgQhull::Qhull qhull;
			qhull.setOutputStream(&output_hull);
			orgQhull::RboxPoints rbox;

			std::stringstream stringStream;
			stringStream << "3 " << vectfaces[i].size();

			for (int j = 0; j < vectfaces[i].size(); ++j) {
				cur_face.append(vectfaces[i][j]);
				for (int n = 0; n < 3; ++n)
				{
					stringStream << ' ' << voro_vertexes[vectfaces[i][j]][n];
				}
			}

			faces.append(cur_face);
			std::istringstream istringStream(stringStream.str());
			rbox.appendPoints(istringStream);
			qhull.runQhull(rbox, "i Qt");
			qhull.outputQhull();

			buff.clear();
			Array godot_triangles = Array();
			bool first_line = true;
			
			for (auto n : output_hull.str()) {
				if (n != '\n') {
					buff += n;
				} else if (first_line) {
					first_line = false;
					buff.clear();
				} else {
					std::string::size_type sz;
					std::string::size_type pos = 0;
					int a = std::stoi(buff, &sz);
					pos += sz;
					int b = std::stoi(buff.substr(pos), &sz);
					pos += sz;
					int c = std::stoi(buff.substr(pos), &sz);

					Array p = Array();
					p.append(Vector3(voro_vertexes[vectfaces[i][a]][0], voro_vertexes[vectfaces[i][a]][1], voro_vertexes[vectfaces[i][a]][2]));
					p.append(Vector3(voro_vertexes[vectfaces[i][b]][0], voro_vertexes[vectfaces[i][b]][1], voro_vertexes[vectfaces[i][b]][2]));
					p.append(Vector3(voro_vertexes[vectfaces[i][c]][0], voro_vertexes[vectfaces[i][c]][1], voro_vertexes[vectfaces[i][c]][2]));
					godot_triangles.append(p);
					buff.clear();
				}
			}
			fragments.append(godot_triangles);
		}
	}
}

Dictionary Voronoi::get_vertexes() {
	return vertexes;
}

Array Voronoi::get_faces() {
	return faces;
}

Array Voronoi::get_fragments() {
	return fragments;
}

void Voronoi::set_seed(double seedInput)
{
	manual_seed = true;
	seed = seedInput;
}

double Voronoi::get_seed() {
	return seed;
}

void Voronoi::_bind_methods() {
	ClassDB::bind_method(D_METHOD("voronoi2d_in_box", "bounding_points", "n_points"), &Voronoi::voronoi2d_in_box);
	ClassDB::bind_method(D_METHOD("voronoi2d"), &Voronoi::voronoi2d);
	ClassDB::bind_method(D_METHOD("voronoi3d_in_box", "bounding_points", "n_points"), &Voronoi::voronoi3d_in_box);
	ClassDB::bind_method(D_METHOD("voronoi3d"), &Voronoi::voronoi3d);
	ClassDB::bind_method(D_METHOD("get_vertexes"), &Voronoi::get_vertexes);
	ClassDB::bind_method(D_METHOD("get_faces"), &Voronoi::get_faces);
	ClassDB::bind_method(D_METHOD("get_fragments"), &Voronoi::get_fragments);
	ClassDB::bind_method(D_METHOD("set_seed"), &Voronoi::set_seed);
	ClassDB::bind_method(D_METHOD("get_seed"), &Voronoi::get_seed);
}

Voronoi::Voronoi() {
}
