#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/RboxPoints.h"

#include "voronoi.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

std::vector<double> Voronoi::generate_random_points2d(int n_points, double min_x, double max_x, double min_y, double max_y) {
	return generate_random_points(n_points, 2, std::vector<double>{ min_x, min_y }, std::vector<double>{ max_x, max_y });
}

std::vector<double> Voronoi::generate_random_points3d(int n_points, const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z) {
	return generate_random_points(n_points, 3, std::vector<double>{ min_x, min_y, min_z }, std::vector<double>{ max_x, max_y, max_z });
}

std::vector<double> Voronoi::generate_random_points(const int n_points, const int n_dim, std::vector<double> val_min, std::vector<double> val_max) {
	std::vector<double> random_points(n_points * n_dim);
	std::mt19937 eng;
	if (manual_seed) {
		eng.seed(seed);
	} else {
		eng.seed(std::random_device()());
	}
	std::uniform_real_distribution<> distr(0.0, 1.0);
	double rand_value = 0;
	for (int i = 0; i < n_points * n_dim; i += n_dim) {
		for (int n = 0; n < n_dim; ++n) {
			rand_value = distr(eng);
			random_points[i + n] = rand_value * (val_max[n] - val_min[n]) + val_min[n];
		}
	}

	return random_points;
}

void Voronoi::voronoi2d(int n_points) {
	Array bounding_points = Array();
	bounding_points.append(Vector2(0, 0));
	bounding_points.append(Vector2(1, 0));
	bounding_points.append(Vector2(1, 1));
	bounding_points.append(Vector2(0, 1));
	Voronoi::voronoi2d_in_box(bounding_points, n_points);
}

void Voronoi::get_boundaries(Array const &bounding_points, int const &n_dim, std::vector<double> &val_min, std::vector<double> &val_max) {

	int size = bounding_points.size();
	val_min = std::vector<double>();
	val_max = std::vector<double>();

	for (int n = 0; n < n_dim; ++n) {
		val_min.push_back(bounding_points.get(0).get(n));
		val_max.push_back(bounding_points.get(0).get(n));
	}

	for (int i = 0; i < size; ++i) {
		for (int n = 0; n < n_dim; ++n) {
			double val = static_cast<double>(bounding_points.get(i).get(n));
			if (val > val_max[n]) {
				val_max[n] = val;
			} else if (val < val_min[n]) {
				val_min[n] = val;
			}
		}
	}
}

void Voronoi::voronoi2d_in_box(Array bounding_points, int n_points) {
	Voronoi::voronoi_in_box(bounding_points, n_points, 2);
}

void Voronoi::voronoi_in_box(Array bounding_points, int n_points, int n_dim) {
	int precision = std::numeric_limits<double>::digits10 + 1;

	std::vector<double> val_min;
	std::vector<double> val_max;
	Voronoi::get_boundaries(bounding_points, n_dim, val_min, val_max);

	//create some random points inside the bounding box
	std::vector<double> r_points = generate_random_points(n_points, n_dim, val_min, val_max);
	//assign the random points to the vector
	std::vector<double> all_points(n_points * n_dim * (2 * n_dim + 1));
	std::copy_n(r_points.begin(), n_points * n_dim, all_points.begin());
	
	//adjust/mirror the points
	std::vector<double> offset = std::vector<double>();
	std::vector<int> prefix = std::vector<int>();

	//TODO move to lower loop
	for (int i = 0; i < n_dim; ++i) {
		for (int j = 0; j < 2 * n_dim; ++j) {
			offset.push_back(2 * (val_min[j % n_dim] * ((j - i) == 0) + val_max[j % n_dim] * ((j - i) == n_dim)));
			prefix.push_back((j % n_dim != i) * 2 - 1);
		}
	}

	//fancy loop which just mirrors the values in all dimensions
	for (int n = 0; n < n_dim * n_dim * 2; ++n) {
		for (int i = 0; i < n_points; ++i) {
			all_points[n_dim * n_points * (n / n_dim + 1) + i * n_dim + n % n_dim] = offset[n] + prefix[n] * all_points[i * n_dim + n % n_dim];
		}
	}

	std::stringstream stringStream;

	stringStream << n_dim << " " << (2 * n_dim + 1) * n_points;
	for (size_t i = 0; i < all_points.size(); i += n_dim) {
		stringStream << "\n";
		for (size_t d = 0; d < n_dim; ++d)
		{
			stringStream << std::setprecision(precision) << all_points[i + d] << " ";
		}
	}
	orgQhull::RboxPoints rbox;
	std::istringstream iStringStream(stringStream.str());
	rbox.appendPoints(std::istringstream(stringStream.str()));
	
	orgQhull::Qhull qhull;
	std::stringstream output;
	qhull.setOutputStream(&output);

	qhull.runQhull(rbox, "v Qbb");
	if (n_dim == 2) {
		qhull.outputQhull("p");
		qhull.outputQhull("FN");
		Voronoi::parse_output2d(output, val_min[0], val_max[0], val_min[1], val_max[1]);
	} else if (n_dim == 3) {
		qhull.outputQhull("o");
		Voronoi::parse_output3d(output, val_min[0], val_max[0], val_min[1], val_max[1], val_min[2], val_max[2]);
	}
}

void Voronoi::voronoi3d(int n_points) {
	Array bounding_points = Array();
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

void Voronoi::voronoi3d_in_box(Array bounding_points, int n_points) {
	Voronoi::voronoi_in_box(bounding_points, n_points, 3);
}

void Voronoi::parse_output2d(std::stringstream &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y) {
	std::vector<bool> good_vertexes;
	std::vector<std::vector<int> > vectfaces;

	double epsilon = 0.000001;
	int n_points = 0;
	int n_faces = 0;
	int line_counter = 0;

	bool good_vertex;
	
	output.seekg(0);
	// output has the following format: Dimension\nNumberOfFaces\nFaces\nNumberOfVertexes\nVertexes
	for (std::string line; std::getline(output, line);) {

		if (line_counter == 1) {
			// 2nd line/2nd part of header states number of points/vertexes
			n_points = atoi(line.c_str());
			good_vertexes.reserve(n_points);
		} else if (line_counter == n_points + 2) {
			// 2nd header states number of faces
			n_faces = atoi(line.c_str());
			vectfaces = std::vector<std::vector<int> >(n_faces);
		} else if (line_counter > 1 && line_counter < n_points + 2) {
			// parse points
			int space_pos = line.find(" ");

			// sometimes there are trailing spaces in the QHull output which we need to skip
			if (space_pos == 0) {
				space_pos = 1;
				while (line.substr(space_pos, 1) == " ") {
					space_pos += 1;
				}
				line = line.substr(space_pos);
				space_pos = line.find(" ");
			}

			double x = atof(line.substr(0, space_pos).c_str());
			double y = atof(line.substr(space_pos).c_str());

			good_vertex = (x >= min_x - epsilon && x <= max_x + epsilon && y >= min_y - epsilon && y <= max_y + epsilon);
			good_vertexes.push_back(good_vertex);

			if (good_vertex) {
				vertexes[line_counter - 2] = Vector2(x, y);
			}
		} else if (line_counter > (n_points + 1)) {
			//parse faces
			int pos = 0;
			int offset = 0;
			int vect_index = line_counter - 3 - n_points;
			while (pos > -1) {
				pos = line.substr(offset).find(" ");
				int value = atoi(line.substr(offset, pos).c_str());

				if (offset == 0) {
					vectfaces[vect_index] = std::vector<int>();
					vectfaces[vect_index].reserve(value);
				} else {
					vectfaces[vect_index].push_back(value);
				}
				offset += pos + 1;
			}
		}
		line_counter++;
	}

	faces = Array();

	bool good_face;
	for (int i = 0; i < n_faces; ++i) {
		good_face = True;
		Array cur_face = Array();
		for (int j = 0; j < vectfaces[i].size(); ++j) {
			if (!good_vertexes[std::abs(vectfaces[i][j])]) {
				good_face = false;
				break;
			}
			cur_face.append(std::abs(vectfaces[i][j]));
		}

		if (good_face) {
			faces.append(cur_face);
		}
	}
}

void Voronoi::parse_output3d(std::stringstream &output, const double &min_x, const double &max_x, const double &min_y, const double &max_y, const double &min_z, const double &max_z) {

	std::vector<bool> good_vertexes;
	std::vector<std::vector<int> > vectfaces;

	double epsilon = 0.000001;
	int n_points = 0;
	int n_faces = 0;
	int line_counter = 0;

	bool good_vertex;

	std::vector<std::vector<double> > voro_vertexes = std::vector<std::vector<double> >();
	output.seekg(0);

	for (std::string line; std::getline(output, line);) {
		std::string::size_type pos = 0;
		if (line_counter == 1) {
			pos = line.find(" ");
			n_points = atoi(line.substr(0, pos).c_str());
			n_faces = atoi(line.substr(pos + 1).c_str());
			vectfaces = std::vector<std::vector<int> >(n_faces);

		} else if (line_counter > 1 && line_counter < n_points + 2) {
			pos = line.find(" ");
			int offset = 0;
			double x = atof(line.substr(0, pos).c_str());
			offset += pos + 1;
			pos = line.substr(offset).find(" ");
			double y = atof(line.substr(offset, pos).c_str());
			offset += pos + 1;
			pos = line.substr(offset).find(" ");
			double z = atof(line.substr(offset, pos).c_str());
			good_vertex = (x >= min_x - epsilon && x <= max_x + epsilon && y >= min_y - epsilon && y <= max_y + epsilon && z >= min_z - epsilon && z <= max_z + epsilon);
			good_vertexes.push_back(good_vertex);
			if (good_vertex) {
				vertexes[line_counter - 2] = Vector3(x, y, z);
			}
			voro_vertexes.push_back(std::vector<double>{ x, y, z });

		} else if (line_counter >= (n_points + 2)) {
			int counter = -1;
			int vect_index = line_counter - 2 - n_points;
			int offset = 0;
			int pos = 0;
			while (pos > -1) {
				pos = line.substr(offset).find(" ");
				int value = atoi(line.substr(offset, pos).c_str());

				if (offset == 0) {
					vectfaces[vect_index] = std::vector<int>();
				} else {
					vectfaces[vect_index].push_back(value);
				}
				offset += pos + 1;
			}
			counter++;
		}
		line_counter++;
	}

	faces = Array();
	fragments = Array();

	for (int i = 0; i < n_faces; ++i) {
		bool good_face = True;
		for (int j = 0; j < vectfaces[i].size(); ++j) {
			if (!good_vertexes[std::abs(vectfaces[i][j])]) {
				good_face = False;
				break;
			}
		}
		if (good_face && vectfaces[i].size() > 0) {
			std::stringstream output_hull;
			std::stringstream stringStream;
			orgQhull::RboxPoints rbox;
			orgQhull::Qhull qhull;
			qhull.setOutputStream(&output_hull);

			Array cur_face = Array();
			stringStream << "3 " << vectfaces[i].size();
			for (int j = 0; j < vectfaces[i].size(); ++j) {
				cur_face.append(vectfaces[i][j]);
				stringStream << "\n" << voro_vertexes[vectfaces[i][j]][0];
				stringStream << " " << voro_vertexes[vectfaces[i][j]][1];
				stringStream << " " << voro_vertexes[vectfaces[i][j]][2];
			}
			stringStream << "\n";
			faces.append(cur_face);

			std::istringstream iStringStream = std::istringstream(stringStream.str());
			rbox.appendPoints(iStringStream);

			qhull.runQhull(rbox, "i Qt QJ Pp"); // we need QJ to make sure that planar points can be calculated
			qhull.outputQhull();

			Array godot_triangles = Array();
			bool first_line = true;
			output_hull.seekg(0);
			for (std::string line; std::getline(output_hull, line);) {
				if (first_line) {
					first_line = false;
				} else {
					std::string::size_type offset = 0;
					std::string::size_type pos = 0;
					pos = line.find(" ");
					int a = atoi(line.substr(offset, pos).c_str());
					offset += pos + 1;
					pos = line.substr(offset).find(" ");
					int b = atoi(line.substr(offset, pos).c_str());
					offset += pos + 1;
					pos = line.substr(offset).find(" ");
					int c = atoi(line.substr(offset, pos).c_str());

					Array p = Array();
					p.append(Vector3(voro_vertexes[vectfaces[i][a]][0], voro_vertexes[vectfaces[i][a]][1], voro_vertexes[vectfaces[i][a]][2]));
					p.append(Vector3(voro_vertexes[vectfaces[i][b]][0], voro_vertexes[vectfaces[i][b]][1], voro_vertexes[vectfaces[i][b]][2]));
					p.append(Vector3(voro_vertexes[vectfaces[i][c]][0], voro_vertexes[vectfaces[i][c]][1], voro_vertexes[vectfaces[i][c]][2]));
					godot_triangles.append(p);
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

void Voronoi::set_seed(double seedInput) {
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
