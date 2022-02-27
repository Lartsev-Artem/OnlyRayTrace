#include <fstream>
#include <iostream>

#include <set>
#include <string>
#include <vector>

#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>

#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>

#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>

#include <eigen3/Eigen/Dense>

typedef double Type;
typedef int IntId;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;

int main(int argc, char* argv[]) {

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\";

	std::string name_file_settings = main_file_direction + "FOOsettings_file.txt";
	std::string name_file_grid = main_file_direction + "Sphere2519.vtk";//"Sphere109.vtk";//"MyVTK1.vtk";//
	std::string name_file_sphere_direction = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph = main_file_direction + "TestClaster\\graph";//"GraphDirections\\graph";

	ReadStartSettings(name_file_settings, main_file_direction, name_file_sphere_direction, name_file_graph);

	std::vector<Type> directions;
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions);

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	ReadGridVtk(name_file_grid, unstructured_grid);
	std::cout << "Grid has cells: " << unstructured_grid->GetNumberOfCells() << '\n';


	std::vector<IntId> all_pairs_id;
	FindNeighborsPairFace(unstructured_grid, all_pairs_id);


	std::vector<IntId> faces_state;
	std::set<IntId> set_boundary_cells;
	std::set<IntId> set_graph;

	Vector3 direction;

	for (size_t cur_direction = 0; cur_direction < directions.size()/2; cur_direction++)
	{
		set_graph.clear();
		InitBoundarySet(unstructured_grid, set_boundary_cells);
		InitFacesState(all_pairs_id, faces_state);
		FromSphericalToDecart(cur_direction, directions, direction);

		//direction = Vector3(1, 0, 0);
	
		std::ofstream ofile;
		ofile.open(name_file_graph + to_string(cur_direction) + ".txt");
		if (!ofile.is_open()) {
			std::cout << "Error open file\n";
			std::cout << "file_graph is not opened for writing\n";
			return 1;
		}

		//std::vector<int> buf;
		while (set_boundary_cells.size()) {

			IntId id_cell = GetNextCell(direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id);

			set_graph.emplace(id_cell);
			//buf.push_back(id_cell);
			ofile << id_cell << '\n';

			ReBuildSetBondary(id_cell, direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id, set_graph);
		}
		ofile.close();

		std::cout << "Graph construction in the direction #" << cur_direction << " is completed\n";
		//WriteFileBoundary(main_file_direction + "Sphere565boundVTK.vtk", name_file_graph + "0.txt", unstructured_grid);
		//return -666;
	}
			
	return 0;
}

bool CheckCell(const IntId id_cell, const Vector3& direction, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id) {


	IntId face_state[4];
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {
			if (faces_state[id_cell * 4 + i] == 0){// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// грани не определены
				return false;
			}
		}
	return true;
}
IntId GetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id) {
	// возвращает (случайный)индекс ячейки со всеми определенными входящими гранями

	for (auto el : boundary_id) {
		if (CheckCell(el, direction, u_grid, faces_state, all_pairs_id)) {
			for (size_t i = 0; i < 4; i++) {
				faces_state[el * 4 + i] = 1;
				if (all_pairs_id[el * 4 + i] != -1)
					faces_state[all_pairs_id[el * 4 + i]] = 1;
			}
			return el;
		}
	}

	std::cout << "NextCell is -1\n";
	return -1;
}
int InitBoundarySet(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::set<IntId>& boundary_cells) {

	boundary_cells.clear();
	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {
				boundary_cells.emplace(i);
				break;
			}
		}
	}

	return 0;
}
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state) {
	
	faces_state.assign(all_pairs_id.size(), 0);
	for (size_t i = 0; i < all_pairs_id.size(); i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}
	return 0;
}
int ReadGridVtk(const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}
	return 0;
}
int ReBuildSetBondary(const IntId id_cell, const Vector3& direction, std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph) {

	IntId face_state[4];
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i] == 0) {
			int num_face = all_pairs_id[id_cell * 4 + i];

			int buf_s = set_graph.size();
			set_graph.emplace(num_face / 4);
			if (set_graph.size() != buf_s) {
				boundary_id.emplace(num_face / 4);
				set_graph.erase(num_face / 4);
			}
		}
	boundary_id.erase(id_cell);

	return 0;
}
size_t WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid) {

	int n = u_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);

	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	int i = 0, el;
	while (ifile >> el)
		bound_array->SetTuple1(el, i++);
	ifile.close();

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = u_grid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(bound_array);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}