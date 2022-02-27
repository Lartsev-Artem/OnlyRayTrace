#pragma once
//#ifndef RAY_TRACING3D_SIMPLE_CALCULATION 
//#define RAY_TRACING3D_SIMPLE_CALCULATION
//#endif 

#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

typedef double Type;

Type CheckLength(const std::vector<std::vector<Type>>& point);
Type Distance(const Type* point1, const Type* point2);
Type Distance(const Type point1_x, const Type point1_y, const Type point1_z,
	const Type point2_x, const Type point2_y, const Type point2_z);

size_t IntersectionWithPlane(const std::vector<Type*>& face, const Type* start_point, const Type* direction, std::vector<Type>& result);
bool InTriangle(const std::vector<Type*>& face, const Type* X);
size_t Make2dPoint(const Type* start, Type**& local_basis, const Type* point, Type* new_point);
Type MakeLength(Type* point1, Type* point2);
Type Norm(const Type* vec);
size_t Normalize(Type* vec);
Type Rosh(const Type* S, Type*& a, Type t);
size_t SetBasis(const Type* start_point, const Type* normal, Type* vec_1, Type* vec_2);
size_t SetDirect(const Type* start, const Type* end, Type* direct);