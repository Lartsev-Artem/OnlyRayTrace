#pragma once

template<typename Type>
size_t SetMemory(size_t n, Type*& first)
{
	first = new Type[n];
	return 0;
}

template<typename Type, typename... Other>
size_t SetMemory(size_t n, Type*& first, Other*&... other)
{
	first = new Type[n];
	SetMemory(n, other...);
	return 0;
}

template<typename Type>
size_t SetMemory(size_t n, Type**& matrix)
{
	matrix = new Type * [n];   //Выделение памяти под массив указателей
	for (size_t i = 0; i < n; ++i) {
		matrix[i] = new Type[n]; // Выделение памяти под массив значений
	}
	return 0;
}

template<typename Type, typename... Other>
size_t SetMemory(size_t n, Type**& first, Other**&... other)
{
	first = new Type * [n];   //Выделение памяти под массив указателей
	for (size_t i = 0; i < n; ++i) {
		first[i] = new Type[n]; // Выделение памяти под массив значений
	}
	SetMemory(n, other...);
	return 0;
}

template<typename Type>
size_t SetMemory(size_t n, size_t m, Type**& matrix)
{
	matrix = new Type * [n];   //Выделение памяти под массив указателей
	for (size_t i = 0; i < n; ++i) {
		matrix[i] = new Type[m]; // Выделение памяти под массив значений
	}
	return 0;
}

template<typename Type, typename... Other>
size_t SetMemory(size_t n, size_t m, Type**& first, Other**&... other)
{
	first = new Type* [n];   //Выделение памяти под массив указателей
	for (size_t i = 0; i < n; ++i) {
		first[i] = new Type[m]; // Выделение памяти под массив значений
	}
	SetMemory(n, m, other...);
	return 0;
}

template<typename Type>
size_t ClearMemory(Type*& first)
{
	delete[]first;
	return 0;
}

template<typename Type, typename... Other>
size_t ClearMemory(Type*& first, Other*&... other)
{
	delete[]first;
	ClearMemory(other...);
	return 0;
}

template<typename Type>
size_t ClearMemory(size_t n, Type**& matrix)
{
	for (size_t i = 0; i < n; ++i)
		delete matrix[i];
	delete[] matrix;
	return 0;
}

template<typename Type, typename... Other>
size_t ClearMemory(size_t n, Type**& first, Other**&... other)
{
	for (size_t i = 0; i < n; ++i)
		delete first[i];
	delete[] first;
	ClearMemory(n, other...);
	return 0;
}


template<typename Type>
size_t SetIdentMatrix(size_t n,Type**& matrix)
{
	matrix = new Type * [n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new Type[n];
		for (int j = 0; j < n; ++j)
			i == j ? matrix[i][j] = 1 : matrix[i][j] = 0;
	}
	return 0;
}