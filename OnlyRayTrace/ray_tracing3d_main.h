#pragma once
#include "ray_tracing3d_vtk_resourse.h"

typedef double Type;

//#define DRAW_RAY
#define CHECK_FULL_TIME

template<typename Type>
size_t ReadStartSettings(std::string name_file_settings, size_t& class_file_vtk, std::string& name_file_vtk, std::string& out_file_grid_vtk, std::string& out_file_ray_vtk,
	Type* center_sphere, Type& radius_sphere, Type& internal_radius_disk, Type& external_radius_disk,
	Type& width_plane, Type& height_plane, int& pixels_wide, int& pixels_high, size_t& number_of_sort_subdomains, size_t& number_of_global_steps) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings is not open !\n";
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	getline(ifile, str);
	getline(ifile, name_file_vtk);
	getline(ifile, out_file_grid_vtk);
	getline(ifile, out_file_ray_vtk);

	for (size_t i = 0; i < 3; ++i) {
		ifile >> center_sphere[i];
	}
	getline(ifile, str);

	ifile >> radius_sphere;
	getline(ifile, str);
	ifile >> internal_radius_disk;
	getline(ifile, str);
	ifile >> external_radius_disk;
	getline(ifile, str);

	ifile >> width_plane;
	getline(ifile, str);
	ifile >> height_plane;
	getline(ifile, str);
	ifile >> pixels_wide;
	getline(ifile, str);
	ifile >> pixels_high;
	getline(ifile, str);
	ifile >> number_of_sort_subdomains;
	getline(ifile, str);
	ifile >> number_of_global_steps;
	getline(ifile, str);


	ifile.close();
	return 0;
}