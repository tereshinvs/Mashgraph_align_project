void Matrix_DSU::make_set(uint i, uint j, bool _edge) {
	dsu(i, j) = { i, j, 1, _edge };
}

void Matrix_DSU::union_set(uint i1, uint j1, uint i2, uint j2) {
	uint x1, y1, x2, y2;
	std::tie(x1, y1) = find_set(i1, j1);
	std::tie(x2, y2) = find_set(i2, j2);
	if (x1 != x2 || y1 != y2) {
		if (dsu(x1, y1).size < dsu(x2, y2).size) {
			std::swap(x1, x2);
			std::swap(y1, y2);
		}
		dsu(x2, y2).x = x1;
		dsu(x2, y2).y = y1;
		dsu(x1, y1).size += dsu(x2, y2).size;
		dsu(x1, y1).edge |= dsu(x2, y2).edge;
	}
}

std::tuple<uint, uint> Matrix_DSU::find_set(uint i, uint j) {
	uint &xp = dsu(i, j).x, &yp = dsu(i, j).y;
	if (xp == i && yp == j)
		return std::make_tuple(i, j);
	uint x_res, y_res;
	std::tie(x_res, y_res) = find_set(xp, yp);
	xp = x_res; yp = y_res;
	return std::make_tuple(x_res, y_res);
}

bool Matrix_DSU::is_edge(uint i, uint j) {
	uint x, y;
	std::tie(x, y) = find_set(i, j);
	return dsu(x, y).edge;
}

void Matrix_DSU::set_edge(uint i, uint j) {
	uint x, y;
	std::tie(x, y) = find_set(i, j);
	dsu(x, y).edge = true;
}