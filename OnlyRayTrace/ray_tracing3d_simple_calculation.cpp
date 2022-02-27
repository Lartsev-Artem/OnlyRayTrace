//#ifndef RAY_TRACING3D_SIMPLE_CALCULATION
#include "ray_tracing3d_simple_calculation.h"
//#endif

/*отладочная функция. Проверяет геометрическую длину пути луча внутри расчетной области*/
Type CheckLength(const std::vector<std::vector<Type>>& point) {
	if (point.size() == 0)
		return 0;
	Type sum = 0;
	Type lenght = 0;
	std::vector<Type> prev_point = point.front();
	
	for (auto el : point) {
		sum = 0;
		for (size_t i = 0; i < 3; ++i) {
			sum += pow(el[i] - prev_point[i], 2);
		}
		lenght += sqrt(sum);
		prev_point = el;
	//	std::cout << el[0] << ' ' << el[1] << ' ' << el[2] << '\n';
	}

	return lenght;
}

Type Distance(const Type* point1, const Type* point2) {
	Type s = 0;
	for (size_t i = 0; i < 3; ++i)
		s += pow((point2[i] - point1[i]), 2);
	return sqrt(s);
}
Type Distance(const Type point1_x, const Type point1_y, const Type point1_z,
	const Type point2_x, const Type point2_y, const Type point2_z) {

	Type s = pow(point1_x - point2_x, 2) + pow(point1_y - point2_y, 2) + pow(point1_z - point2_z, 2);
	return sqrt(s);
}

size_t IntersectionWithPlane(const std::vector<Type*>& face, const Type* start_point, const Type* direction, std::vector<Type>& result) {

	//вершины треугольника
	Type* A = face[0];
	Type* B = face[1];
	Type* C = face[2];

	Type a, b, c, d;  // параметры уравнения плоскости
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}

bool InTriangle(const std::vector<Type*>& face, const Type* X) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type* A = face[0];
	Type* B = face[1];
	Type* C = face[2];

	// линейная алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

size_t Make2dPoint(const Type* start, Type**& local_basis, const Type* point, Type* new_point) {

	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis[0][k];
		new_point[1] += (point[k] - start[k]) * local_basis[1][k];
	}
	return 0;
}

Type MakeLength(Type* point1, Type* point2) {
	
	Type s = 0;  // безразмерная длина

	for (size_t i = 0; i < 3; ++i)
		s += (point2[i] - point1[i]) * (point2[i] - point1[i]);

	return sqrt(s);  // *388189 * pow(10, 5);  // числа --- переход к размерным параметрам
}

Type Norm(const Type* vec) {
	Type sum = 0;
	for (size_t i = 0; i < 3; ++i)
		sum += pow(vec[i], 2);
	return sqrt(sum);
}
size_t Normalize(Type* vec) {

	Type norm = Norm(vec);
	for (int i = 0; i < 3; i++)
		vec[i] /= norm;
	return 0;
}

Type Rosh(const Type* S, Type*& a, Type t) {

	// wolfram расчет
	Type max_potential = -1.82547;  // потенциал роша
	Type m = 0.8795180722891566;

	// математические расчёты
	Type x1 = S[0] + a[0] * t;
	Type x2 = pow(S[1] + a[1] * t, 2);
	Type x3 = pow(S[2] + a[2] * t, 2);

	return  -max_potential - (pow(x1 - m, 2) + x2 + x3) / 2 -
		m / sqrt(pow(x1 - 1, 2) + x2 + x3) -
		(1 - m) / sqrt(pow(x1, 2) + x2 + x3);
}

size_t SetBasis(const Type* start_point, const Type* normal, Type* vec_1, Type* vec_2) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/

	const Type end_point[3] = { 0.879518, 0, start_point[2] };  // параметры подобраны из геометрии задачи!!

	Type N[3];
	SetDirect(start_point, end_point, N);

	if (abs(N[1]) < pow(10, -20)) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(N[0] * vec_1[0] + N[2] * vec_1[2]) / N[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
	Eigen::Vector3d a(normal[0], normal[1], normal[2]);
	Eigen::Vector3d b(vec_1[0], vec_1[1], vec_1[2]);
	Eigen::Vector3d c = a.cross(b);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	Normalize(vec_1);
	Normalize(vec_2);

	return 0;
}

size_t SetDirect(const Type* start, const Type* end, Type* direct) {

	for (size_t i = 0; i < 3; ++i)
		direct[i] = (end[i] - start[i]);
	Normalize(direct);
	return 0;
}


