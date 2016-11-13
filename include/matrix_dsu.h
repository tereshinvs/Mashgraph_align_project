#pragma once

#include "io.h"

class Matrix_DSU {
	private:
		struct DSU_Value {
			uint x, y, size;
			bool edge;
		};

		Matrix<DSU_Value> dsu;

	public:
		Matrix_DSU(uint rows, uint cols): dsu(rows, cols) {}

		void make_set(uint i, uint j, bool _edge);

		void union_set(uint i1, uint j1, uint i2, uint j2);

		std::tuple<uint, uint> find_set(uint i, uint j);

		bool is_edge(uint i, uint j);

		void set_edge(uint i, uint j);
};

#include "matrix_dsu.hpp"