#include <tuple>

template<typename T>
double get_brightness(const std::tuple<T, T, T> &p) {
	double r, g, b;
	std::tie(r, g, b) = p;
	return 0.2125 * r + 0.7154 * g + 0.0721 * b;
}

bool is_correct_point(int i, int j, int rows, int cols) {
	return i >= 0 && i < rows && j >= 0 && j < cols;
}

double normalize_colour(double c, char norm_type) {
	return norm_type == 'S'
		? std::max(0.0, std::min(255.0, c))
		: std::min(255.0, std::abs(c));
}

const std::tuple<uint, uint, uint> &get_reflected(const Image &src, int i, int j) {
	if (i < 0) i = -i;
	if (j < 0) j = -j;
	if (i >= int(src.n_rows)) i = src.n_rows - (i - src.n_rows + 1);
	if (j >= int(src.n_cols)) j = src.n_cols - (j - src.n_cols + 1);
	return src(i, j);
}

uint get_median(const std::vector<int> &num, const int size) {
	int alr = 0;
	for (uint i = 0; i < num.size(); i++) {
		alr += num[i];
		if (alr > size / 2)
			return i;
	}
	return 0U;
}

double MSE::operator()(const Image &img1, const Image &img2, int shift_x, int shift_y) const {
	double res = 0;
	uint cor_p = 0;
	int rows = img1.n_rows, cols = img1.n_cols;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			if (is_correct_point(i + shift_x, j + shift_y, img2.n_rows, img2.n_cols)) {
				res += std::pow(get_brightness<uint>(img1(i, j)) -
					get_brightness<uint>(img2(i + shift_x, j + shift_y)), 2);
				cor_p++;
			}
	return res / cor_p;
}

double Cross_Correlation::operator()(const Image &img1, const Image &img2, int shift_x, int shift_y) const {
	double res = 0;
	int rows = img1.n_rows, cols = img1.n_cols;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			if (is_correct_point(i + shift_x, j + shift_y, img2.n_rows, img2.n_cols))
				res += get_brightness<uint>(img1(i, j)) *
					get_brightness<uint>(img2(i + shift_x, j + shift_y));
	return res;
}

double Gradient_Module_OP::operator()(
	const Matrix<std::tuple<double, double, double>> &m1,
	const Matrix<std::tuple<double, double, double>> &m2) const {
	double br_x = get_brightness<double>(m1(0, 0)), br_y = get_brightness<double>(m2(0, 0));
	return std::sqrt(br_x * br_x + br_y * br_y);
}

double Gradient_Direction_OP::operator()(
	const Matrix<std::tuple<double, double, double>> &m1,
	const Matrix<std::tuple<double, double, double>> &m2) const {
	double br_x = get_brightness<double>(m1(0, 0)), br_y = get_brightness<double>(m2(0, 0));
	return std::atan2(br_y, br_x);
}

double Non_Maximum_Supression_OP::operator()(const Matrix<double> &G, const Matrix<double> &O) const {
	uint x1 = 0, y1 = 0, x2 = 0, y2 = 0;
	std::tie(x1, y1) = get_neighbour(O(0, 0));
//	std::tie(x2, y2) = get_neighbour(O(0, 0) > 0 ? O(0, 0) - PI : (O(0, 0) + PI));
	x2 = 2 - x1; y2 = 2 - y1;
	return (G(1, 1) <= G(x1, y1) || G(1, 1) <= G(x2, y2)) ? -1 : G(1, 1);
}

std::tuple<uint, uint> Non_Maximum_Supression_OP::get_neighbour(const double o) const {
	if (-7 * PI / 8.0 <= o && o <= -5 * PI / 8.0)
		return std::make_tuple(2, 0);
	if (-5 * PI / 8.0 <= o && o <= -3 * PI / 8.0)
		return std::make_tuple(2, 1);
	if (-3 * PI / 8.0 <= o && o <= -PI / 8.0)
		return std::make_tuple(2, 2);
	if (-PI / 8.0 <= o && o <= PI / 8.0)
		return std::make_tuple(1, 2);
	if (PI / 8.0 <= o && o <= 3 * PI / 8.0)
		return std::make_tuple(0, 2);
	if (3 * PI / 8.0 <= o && o <= 5 * PI / 8.0)
		return std::make_tuple(0, 1);
	if (5 * PI / 8.0 <= o && o <= 7 * PI / 8.0)
		return std::make_tuple(0, 0);
	return std::make_tuple(1, 0);
}

uint Tresholding_OP::operator()(const Matrix<double> &m) const {
	if (m(0, 0) < treshold1)
		return 0U;
	else if (m(0, 0) <= treshold2)
		return 1U;
	else
		return 2U;
}

std::tuple<uint, uint, uint> R_plus_G::operator()(const Image &m1, const Image &m2) const {
	return std::make_tuple(get_brightness<uint>(m1(0, 0)), get_brightness<uint>(m2(0, 0)), 0);
}

std::tuple<uint, uint, uint> RG_plus_B::operator()(const Image &m1, const Image &m2) const {
	double r, g, b;
	std::tie(r, g, b) = m1(0, 0);
	return std::make_tuple(r, g, get_brightness<uint>(m2(0, 0)));
}

std::tuple<uint, uint, uint> Make_RGB_OP::operator()(const Matrix<std::tuple<double, double, double>> &m) const {
	double r, g, b;
	std::tie(r, g, b) = m(0, 0);
	return std::make_tuple(normalize_colour(r, norm_type),
		normalize_colour(g, norm_type),
		normalize_colour(b, norm_type));
}