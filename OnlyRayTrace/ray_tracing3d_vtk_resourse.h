#pragma once

#include <algorithm>
#include <ctime>
#include <execution>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>


#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkDataSet.h>
#include <vtk-9.0\vtkDataSetAttributes.h>
#include <vtk-9.0\vtkDataObject.h>
#include <vtk-9.0\vtkDoubleArray.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkLine.h>
#include <vtk-9.0\vtkLineSource.h>
#include <vtk-9.0\vtkMath.h>
#include <vtk-9.0\vtkNamedColors.h>
#include <vtk-9.0\vtkPointData.h>
#include <vtk-9.0\vtkPoints.h>
#include <vtk-9.0\vtkQuad.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkTetra.h>
#include <vtk-9.0\vtkTriangle.h>
#include <vtk-9.0\vtkUnsignedCharArray.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>

#include <eigen3/Eigen/Dense>

#include "ray_tracing3d_simple_calculation.h"
#include "memory_pack.h"

#define PI 3.14159265358979323846
typedef double Type;

size_t Build2dPlane(const Type* start_ray, Type**& local_basis, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkSmartPointer<vtkPoints>& points_of_plane,
	vtkSmartPointer<vtkUnstructuredGrid>& plane_grid);

size_t CheckIntersectionWithGeometryObjects(Type* normal_to_picture_plane, const Type* start_ray, const Type* center_sphere, const Type radius_sphere,
	const Type internal_radius_disk, const Type external_radius_disk);

size_t FindNumberOfRoots(std::vector<Type>& roots, const size_t n, const  Type a, const Type b, const Type* start, Type*& direction, Type(*f)(const Type*, Type*&, Type));

size_t FindIdAndPoints(const std::vector<vtkIdType>& set_subdomain_indexes, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid, const vtkSmartPointer<vtkPoints>& points_of_plane,
	const Type* ray_position_on_plane, Type* normal_to_picture_plane, const Type* start_ray,
	std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id, size_t& number_of_intersections);

Type MakeIllum(vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id,
	size_t& index, bool is_intersection_with_rosh);

size_t OneFullTraceStep(const size_t global_step, const size_t number_of_global_steps,
	const Type* center_sphere, const Type radius_sphere, const Type internal_radius_disk, const Type external_radius_disk,
	const Type width_plane, const Type height_plane, const int pixels_wide, const int pixels_high, const size_t number_of_sort_subdomains,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate,
	std::string  out_file_grid_vtk, std::string  out_file_ray_vtk);

size_t RayTracingAndBuildImage(const vtkSmartPointer<vtkPoints>& points_of_plane, const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid,
	vtkSmartPointer<vtkUnstructuredGrid>& image_plane_grid, vtkSmartPointer<vtkUnstructuredGrid>& ray_grid,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate,
	const size_t number_of_sort_subdomains, const size_t max_subdomain, const  std::vector<std::vector<vtkIdType>>& set_subdomain_indexes,
	const Type* start_ray, Type**& local_basis, Type*& normal_to_picture_plane, const Type width_plane, const Type height_plane, const int pixels_wide, const int pixels_high,
	const Type* center_sphere, const Type radius_sphere, const Type internal_radius_disk, const Type external_radius_disk);

size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print = false);

size_t Sort2dGrid(const size_t N/*число подобластей по каждому направлению*/, const Type width_plane, const Type height_plane, const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid,
	const vtkSmartPointer<vtkPoints>& points_of_plane, std::vector<std::vector<vtkIdType>>& set_subdomain_indexes);

size_t SortIntersectionPoints(const Type* start_ray, size_t& number_of_intersections, std::vector<std::pair<std::vector<Type>, vtkIdType >>& points_and_id);

size_t WriteFileVTK(const std::string out_file_grid_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);