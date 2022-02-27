#include "ray_tracing3d_main.h"

int Fullmain(int argc, char* argv[])
{
	std::string name_file_settings;

	if (argc > 1)
		name_file_settings = argv[1];
	else
		name_file_settings = "ray_tracing3d_start_settings_file.txt";

	size_t class_file_vtk;  // задает тип файла --- имена скалярных данных на сетке
	std::string name_file_vtk;
	std::string out_file_grid_vtk;
	std::string out_file_ray_vtk;

	// параметры аккретора и аккреционного диска вне расчётной области
	Type center_sphere[3];
	Type radius_sphere;
	Type internal_radius_disk;
	Type external_radius_disk;

	// параметры картинной плоскости
	Type width_plane;  // безразмерная ширина картинной плоскости
	Type height_plane;  // безразмерная высота картинной плоскости
	int pixels_wide;
	int pixels_high;

	size_t number_of_sort_subdomains;
	size_t number_of_global_steps;  // количесвто конечных кадров 

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, out_file_grid_vtk, out_file_ray_vtk,
		center_sphere, radius_sphere, internal_radius_disk, external_radius_disk,
		width_plane, height_plane, pixels_wide, pixels_high, number_of_sort_subdomains, number_of_global_steps)) {

		std::cout << "Error reading the start settings\n";
		return 1;
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// скалярные данные сетки (unstructured_grid)
	vtkDataArray* density;
	vtkDataArray* absorp_coef;
	vtkDataArray* rad_en_loose_rate;

	clock_t start_clock = clock();
	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	clock_t end_clock = clock();
	std::cout << "\n Reading time of the vtk file: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";

	// один шаг --- одно изображение на картинной плоскости
	for (size_t global_step = 0; global_step < number_of_global_steps; ++global_step) {
		std::cout << "Global step number: " << global_step << '\n';
		
		start_clock = clock();
		OneFullTraceStep(global_step, number_of_global_steps, center_sphere, radius_sphere, internal_radius_disk, external_radius_disk,
			width_plane, height_plane, pixels_wide, pixels_high, number_of_sort_subdomains, unstructured_grid,
			density, absorp_coef, rad_en_loose_rate, out_file_grid_vtk, out_file_ray_vtk);
		clock_t end_clock = clock();
		std::cout << "\n Global step time: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";
	}

	return EXIT_SUCCESS;
}

int CallTrace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, const std::string out_file);

int main1() {

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\";

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkDataArray* foo1;
	vtkDataArray* foo2;
	vtkDataArray* foo3;

	clock_t start_clock = clock();
	if (ReadFileVtk(0, main_file_direction+"Sphere565.vtk", unstructured_grid, foo1, foo2, foo3, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	clock_t end_clock = clock();
	std::cout << "\n Reading time of the vtk file: " << ((Type)end_clock - start_clock) / CLOCKS_PER_SEC << "\n";

	CallTrace(unstructured_grid, main_file_direction + "trace.vtk");
	return 0;
}


size_t SetBasis(const Type* start_point, Eigen::Vector3d& normal, Eigen::Matrix3d& basis) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/
	Eigen::Vector3d vec_1;
	Eigen::Vector3d vec_2;

	if (fabs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Eigen::Matrix3d& local_basis, const Type* point, Eigen::Vector3d& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}
size_t NormalToFace(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, int number_face, Eigen::Vector3d& n) {

	vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(num_cell)->GetFace(number_face)->GetPointIds();

	Type P0[3], P1[3], P2[3];

	unstructured_grid->GetPoints()->GetPoint(idp->GetId(0), P0);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(1), P1);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];


	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructured_grid->GetCell(num_cell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++) {
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructured_grid->GetCell(num_cell)->GetPoints()->GetPoint(idp2->GetId(id), P3);


	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}
bool InTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type AA[3];
	Type BB[3];
	Type CC[3];
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(0, AA);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(1, BB);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(2, CC);


	 Eigen::Vector3d A, B, C, X;
	{
		 Type Xx[3] =  { XX[0], XX[1], XX[2]};
		Eigen::Vector3d n = { 1,0,0 };
		Eigen::Matrix3d basis;
		SetBasis(Xx, n, basis);
		//cout << basis << '\n';
		Make2dPoint(Xx, basis, AA, A);
		Make2dPoint(Xx, basis, BB, B);
		Make2dPoint(Xx, basis, CC, C);
		Make2dPoint(Xx, basis, Xx, X);
	}

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

size_t NewFindIdAndPoints(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const vtkSmartPointer<vtkUnstructuredGrid>& plane_grid, const vtkSmartPointer<vtkPoints>& points_of_plane, Type* normal_to_picture_plane, const Type* start_ray,
	std::vector< vtkIdType>& points_and_id) {

	const int n = plane_grid->GetNumberOfCells();

	// текущая грань треугольника (3 вершины, 3 координаты)
	std::vector<Type*> cur_face_of_triangle;
	cur_face_of_triangle.resize(3);
	for (size_t i = 0; i < 3; ++i)
		cur_face_of_triangle[i] = new Type[3];

	vtkSmartPointer<vtkIdList> idp =
		vtkSmartPointer<vtkIdList>::New();  // буферная переменная

	const Type ray_position_on_plane[3] = { 0,0,0 };

	for (vtkIdType cur_index = 0; cur_index < n; ++cur_index) {

		idp = plane_grid->GetCell(cur_index)->GetPointIds();

		// выбор текущего треугольника 2d
		for (size_t h = 0; h < 3; ++h)
			points_of_plane->GetPoint(idp->GetId(h), cur_face_of_triangle[h]);

		Type a[3] = { cur_face_of_triangle[0][0],cur_face_of_triangle[0][1] ,cur_face_of_triangle[0][2] };
		Type b[3] = { cur_face_of_triangle[1][0],cur_face_of_triangle[1][1] ,cur_face_of_triangle[1][2] };
		Type c[3] = { cur_face_of_triangle[2][0],cur_face_of_triangle[2][1] ,cur_face_of_triangle[2][2] };

		if (InTriangle(cur_face_of_triangle, ray_position_on_plane)/*локальные координаты*/) {

			points_and_id.push_back(cur_index / 4);  // точка пересечения и ячейка на ИСХОДНОЙ сетке
		}
	}

	for (size_t i = 0; i < 3; ++i)
		delete[] cur_face_of_triangle[i];

	return 0;
}

size_t WriteFileSolutionId(const std::string NameFileOut, const std::vector<int>& id, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
size_t FindTraceOfRay(const Type* start_ray, const Type* end_ray,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& id_cells) {

	// вектора локального базиса
	Type normal_to_picture_plane[3];
	Type vector_1[3];
	Type vector_2[3];

	Type** local_basis;
	SetMemory(3, local_basis);

	SetDirect(start_ray, end_ray, normal_to_picture_plane);
	SetBasis(start_ray, normal_to_picture_plane, vector_1, vector_2);

	for (size_t i = 0; i < 3; ++i) {
		local_basis[0][i] = vector_1[i];
		local_basis[1][i] = vector_2[i];
		local_basis[2][i] = normal_to_picture_plane[i];
	}

	// картинная плоскость
	vtkSmartPointer<vtkPoints> points_of_plane =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnstructuredGrid> plane_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// проецирование области на плоскость

	Build2dPlane(start_ray, local_basis, unstructured_grid, points_of_plane, plane_grid);

	ClearMemory(3, local_basis);

	// трассировка и построение картинной плоскости
	
	NewFindIdAndPoints(unstructured_grid, plane_grid, points_of_plane, normal_to_picture_plane, start_ray, id_cells);
	return 0;
}


size_t WriteFileSolutionId(const std::string NameFileOut, const std::vector<int>& id, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {

	int n = UGrid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	IllumArray->SetNumberOfTuples(n);
	for (size_t i = 0; i < n; i++)
		IllumArray->SetTuple1(i, 0);

	for (size_t i = 0; i < id.size(); i++)
		IllumArray->SetTuple1(id[i], 1);

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = UGrid;
	ungrid->GetCellData()->SetActiveScalars("trace");
	ungrid->GetCellData()->SetScalars(IllumArray);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(NameFileOut.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}

int CallTrace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, const std::string out_file) {
	Type start_ray[3] = { -1,0.1,0.1 };
	Type end_ray[3] = { 2,0.1,0.1 };
	std::vector<int> id_cells;

	//FindTraceOfRay(start_ray, end_ray, unstructured_grid, id_cells);
	for (size_t i = 0; i < unstructured_grid->GetNumberOfCells(); i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			if (InTriangle(i, unstructured_grid, unstructured_grid->GetCell(i), j, Eigen::Vector3d(-1, 0,0))){
				id_cells.push_back(i);
				break;
			}
		}
	}



	WriteFileSolutionId(out_file, id_cells, unstructured_grid);
	return 0;
}



