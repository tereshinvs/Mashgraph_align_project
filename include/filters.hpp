#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <functional>
#include <climits>

#include "matrix_dsu.h"
#include "image_tools.h"

class Custom_Filter_OP {
	public:
		const int radius;
		const bool need_div;

		Custom_Filter_OP(const Matrix<double> &_kernel, bool _need_div = false):
			radius(std::max(_kernel.n_rows, _kernel.n_cols) / 2), need_div(_need_div),
			kernel(_kernel.deep_copy()), sum_k(0.0),
			c_i_r(radius), c_i_c(radius), c_k_r(kernel.n_rows / 2), c_k_c(kernel.n_cols / 2) {
			for (uint i = 0; i < kernel.n_rows; i++)
				for (uint j = 0; j < kernel.n_cols; j++)
					sum_k += kernel(i, j);
		}

		std::tuple<double, double, double> operator()(const Image &m) const {
			double sum_r = 0, sum_g = 0, sum_b = 0, r, g, b;
			for (int i = 0; i < int(kernel.n_rows); i++)
				for (int j = 0; j < int(kernel.n_cols); j++) {
					std::tie(r, g, b) = m(c_i_r + i - c_k_r, c_i_c + j - c_k_c);
					sum_r += r * kernel(i, j);
					sum_g += g * kernel(i, j);
					sum_b += b * kernel(i, j);
				}
			if (need_div) {
				sum_r /= sum_k; sum_g /= sum_k; sum_b /= sum_k;
			}
			return std::make_tuple(sum_r, sum_g, sum_b);
		}

	private:
		Matrix<double> kernel;
		double sum_k;
		int c_i_r, c_i_c, c_k_r, c_k_c; // Centres of input image-part and kernel
};

Matrix<std::tuple<double, double, double>>
custom_double(const Image &src, const Matrix<double> &kernel, const bool need_div = false) {
	int radius = std::max(kernel.n_rows, kernel.n_cols) / 2;
	Image big = src.edge_reflection(radius);
	return big.unary_map(Custom_Filter_OP(kernel, need_div)).submatrix(radius, radius, src.n_rows, src.n_cols);	
}

Image custom(const Image &src, const Matrix<double> &kernel, const char norm_type, const bool need_div) {
	return custom_double(src, kernel, need_div).unary_map(Make_RGB_OP(norm_type));
}

Matrix<std::tuple<double, double, double>> sobel_x_double(const Image &src) {
	return custom_double(src, { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} });
}

Matrix<std::tuple<double, double, double>> sobel_y_double(const Image &src) {
	return custom_double(src, { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} });
}

Image sobel_x(const Image &src, const char norm_type) {
	return custom(src, { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} }, norm_type);
}

Image sobel_y(const Image &src, const char norm_type) {
	return custom(src, { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} }, norm_type);
}

Image unsharp(const Image &src) {
	return custom(src, { {-1/6.0, -2/3.0, -1/6.0}, {-2/3.0, 13/3.0, -2/3.0}, {-1/6.0, -2/3.0, -1/6.0} }, 'S', true);
}

Image gray_world(const Image &src) {
	unsigned long long sum_r = 0, sum_g = 0, sum_b = 0;
	for (uint i = 0; i < src.n_rows; i++)
		for (uint j = 0; j < src.n_cols; j++) {
			int r, g, b;
			std::tie(r, g, b) = src(i, j);
			sum_r += r; sum_g += g; sum_b += b;
		}

	double N = src.n_rows * src.n_cols;
	double avg = (sum_r + sum_g + sum_b) / (3.0 * N);
	double r_m = avg * N / sum_r, g_m = avg * N / sum_g, b_m = avg * N / sum_b;

	Image res(src.n_rows, src.n_cols);
	for (uint i = 0; i < src.n_rows; i++)
		for (uint j = 0; j < src.n_cols; j++) {
			int r, g, b;
			std::tie(r, g, b) = src(i, j);
			res(i, j) = std::make_tuple(
				normalize_colour(std::isfinite(r_m) ? r * r_m : avg),
					normalize_colour(std::isfinite(g_m) ? g * g_m : avg),
					normalize_colour(std::isfinite(b_m) ? b * b_m : avg));
		}
	return res;
}

Image resize(const Image &src, const double scale) {
	auto get_colour = [](const double C[4], // 1st point is left-up, 2nd is right-up, 3rd is left-down, 4th is right-down
			double x, double y) {
			return C[0] * (1 - x) * (1 - y) + C[1] * x * (1 - y) + C[2] * (1 - x) * y + C[3] * x * y;
		};

	Image big = src.edge_reflection(1);
	Image res(src.n_rows * scale, src.n_cols * scale);
	for (uint i = 0; i < res.n_rows; i++)
		for (uint j = 0; j < res.n_cols; j++) {
			double old_y = i / scale;
			double old_x = j / scale;

			uint old_col = std::floor(old_x);
			uint old_row = std::floor(old_y);

			old_x -= old_col; old_y -= old_row;
			old_col += 1; old_row += 1;

			// 1st point is left-up, 2nd is right-up,
			// 3rd is left-down, 4th is right-down
			double R[4], G[4], B[4];
			std::tie(R[0], G[0], B[0]) = big(old_row, old_col);
			std::tie(R[1], G[1], B[1]) = big(old_row, old_col + 1);
			std::tie(R[2], G[2], B[2]) = big(old_row + 1, old_col);
			std::tie(R[3], G[3], B[3]) = big(old_row + 1, old_col + 1);

			double r = get_colour(R, old_x, old_y);
			double g = get_colour(G, old_x, old_y);
			double b = get_colour(B, old_x, old_y);
			res(i, j) = std::make_tuple(r, g, b);
		}
	return res;
}

Image autocontrast(const Image &src, const double fraction) {
	auto get_colour = [](double clr, double Y_max, double Y_min) -> double {
		return (clr - Y_min) * 255 / (Y_max - Y_min);
	};
	std::vector<double> Y; // Vector for counting brightness
	for (uint i = 0; i < src.n_rows; i++)
		for (uint j = 0; j < src.n_cols; j++)
			Y.push_back(get_brightness<uint>(src(i, j)));
	std::sort(Y.begin(), Y.end());
	double Y_min = Y[Y.size() * fraction];
	double Y_max = Y[Y.size() - 1 - Y.size() * fraction];
	bool good_hist = std::abs(Y_max - Y_min) > 10e-6;
	Image res(src.n_rows, src.n_cols);
	for (uint i = 0; i < src.n_rows; i++)
		for (uint j = 0; j < src.n_cols; j++) {
			double r, g, b;
			std::tie(r, g, b) = src(i, j);
			double y = get_brightness<uint>(src(i, j));
			if (y < Y_min)
				res(i, j) = std::make_tuple(0, 0, 0);
			else if (y > Y_max)
				res(i, j) = std::make_tuple(255, 255, 255);
			else if (good_hist)
				res(i, j) = std::make_tuple(normalize_colour(get_colour(r, Y_max, Y_min)),
					normalize_colour(get_colour(g, Y_max, Y_min)), 
					normalize_colour(get_colour(b, Y_max, Y_min)));
			else
				res(i, j) = std::make_tuple(r, g, b);
		}
	return res;
}

Image gaussian(const Image &src, const double sigma, const int radius) {
	Matrix<double> G(2 * radius + 1, 2 * radius + 1);
	for (int i = -radius; i <= radius; i++)
		for (int j = -radius; j <= radius; j++)
			G(i + radius, j + radius) = std::exp(-(std::pow(i, 2) + std::pow(j, 2)) /
				(2 * std::pow(sigma, 2))) / (2 * PI * sigma * sigma);
	return custom(src, G, 'S', true);
}

Image gaussian_separable(const Image &src, const double sigma, const int radius) {
	Matrix<double> GH(1, 2 * radius + 1), GV(2 * radius + 1, 1);
	for (int i = -radius; i <= radius; i++)
		GH(0, i + radius) = GV(i + radius, 0) = std::exp(-std::pow(i, 2) /
			(2 * std::pow(sigma, 2))) /	(2 * PI * sigma * sigma);
	return custom(custom(src, GH, 'S', true), GV, 'S', true);
}

Image median(const Image &src, const int radius) {
	Image big = src.edge_reflection(radius);
	Image res(src.n_rows, src.n_cols);
	uint size = (2 * radius + 1) * (2 * radius + 1);
	std::vector<int> r_num(256), g_num(256), b_num(256);
	for (uint i = radius; i < big.n_rows - radius; i++)
		for (uint j = radius; j < big.n_cols - radius; j++) {
			r_num.assign(256, 0); g_num.assign(256, 0); b_num.assign(256, 0);
			for (int i1 = -radius; i1 <= radius; i1++)
				for (int j1 = -radius; j1 <= radius; j1++) {
					uint r, g, b;
					std::tie(r, g, b) = big(i + i1, j + j1);
					r_num[r]++; g_num[g]++; b_num[b]++;
				}
			res(i - radius, j - radius) = std::make_tuple(
				get_median(r_num, size), get_median(g_num, size), get_median(b_num, size));
		}
	return res;
}

Image median_linear(const Image &src, const int radius) {
	Image big = src.edge_reflection(radius);
	Image res(src.n_rows, src.n_cols);

	std::vector<int> l_r_hist(256), l_g_hist(256), l_b_hist(256);
	int size = 2 * radius + 1;
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			uint r, g, b;
			std::tie(r, g, b) = big(i, j);
			l_r_hist[r]++; l_g_hist[g]++; l_b_hist[b]++;
		}

	for (int i = radius; i < int(big.n_rows) - radius; i++) {
		std::vector<int> r_hist = l_r_hist, g_hist = l_g_hist, b_hist = l_b_hist;

		for (int j = radius; j < int(big.n_cols) - radius; j++) {
			res(i - radius, j - radius) = std::make_tuple(
				get_median(r_hist, size * size), get_median(g_hist, size * size), get_median(b_hist, size * size));

			if (j < int(big.n_cols) - radius - 1)
				for (int k = -radius; k <= radius; k++) {
					uint r, g, b;
					std::tie(r, g, b) = big(i + k, j - radius);
					r_hist[r]--; g_hist[g]--; b_hist[b]--;
					std::tie(r, g, b) = big(i + k, j + radius + 1);
					r_hist[r]++; g_hist[g]++; b_hist[b]++;
				}
		}

		if (i < int(big.n_rows) - radius - 1)
			for (int k = 0; k < size; k++) {
				uint r, g, b;
				std::tie(r, g, b) = big(i - radius, k);
				l_r_hist[r]--; l_g_hist[g]--; l_b_hist[b]--;
				std::tie(r, g, b) = big(i + radius + 1, k);
				l_r_hist[r]++; l_g_hist[g]++; l_b_hist[b]++;
			}
	}
	return res;
}

Image canny(const Image &src, const double treshold1, const double treshold2) {
	// 1st stage
	Image img = gaussian_separable(src, 1.4, 2);

	// 2nd stage - counting brightness derivative params
	Matrix<std::tuple<double, double, double>> Ix = sobel_x_double(img), Iy = sobel_y_double(img);
	//save_image(Ix, "Ix.bmp"); save_image(Iy, "Iy.bmp");
	auto G = binary_map<std::tuple<double, double, double>, Gradient_Module_OP>(Gradient_Module_OP(), Ix, Iy);
	auto O = binary_map<std::tuple<double, double, double>, Gradient_Direction_OP>(Gradient_Direction_OP(), Ix, Iy);

	// 3rd stage - non-maximum supression
	G = binary_map<double, Non_Maximum_Supression_OP>(Non_Maximum_Supression_OP(), G, O);

	// 4th stage - supression with tresholds
	auto edge = G.unary_map(Tresholding_OP(treshold1, treshold2));

	// 5th stage - hysteresis tresholding
	Matrix_DSU dsu(img.n_rows, img.n_cols);
	int n_rows = int(img.n_rows), n_cols = int(img.n_cols);
	for (int i = 0; i < n_rows; i++)
		for (int j = 0; j < n_cols; j++) {
			if (edge(i, j) == 0)
				continue;
			dsu.make_set(i, j, edge(i, j) == 2);
			if (is_correct_point(i, j - 1, n_rows, n_cols) && edge(i, j - 1) > 0)
				dsu.union_set(i, j, i, j - 1);
			if (is_correct_point(i - 1, j - 1, n_rows, n_cols) && edge(i - 1, j - 1) > 0)
				dsu.union_set(i, j, i - 1, j - 1);
			if (is_correct_point(i - 1, j, n_rows, n_cols) && edge(i - 1, j) > 0)
				dsu.union_set(i, j, i - 1, j);
			if (is_correct_point(i - 1, j + 1, n_rows, n_cols) && edge(i - 1, j + 1) > 0)
				dsu.union_set(i, j, i - 1, j + 1);
		}

	Image res(img.n_rows, img.n_cols), res1(img.n_rows, img.n_cols);
	for (uint i = 0; i < res.n_rows; i++)
		for (uint j = 0; j < res.n_cols; j++) {
			res1(i, j) = std::make_tuple(125 * edge(i, j), 125 * edge(i, j), 125 * edge(i, j));
			if (edge(i, j) > 0 && dsu.is_edge(i, j))
				res(i, j) = std::make_tuple(255, 255, 255);
			else
				res(i, j) = std::make_tuple(0, 0, 0);
		}
	save_image(res1, "edge.bmp");
	return res;
}

std::tuple<uint, uint, uint, uint> get_frames_value(Image &src, const double near) { // up, right, down, left
	uint res1, res2, res3, res4;

	uint border_rows = src.n_rows * near;
	uint border_cols = src.n_cols * near;

	auto find_maximum_row = [&](const Image &img, uint start, uint end) {
		uint res_num = start == 0 ? start : end, res = 0;
		for (uint i = start; i < end; i++) {
			uint tmp = 0;
			for (uint j = 0; j < img.n_cols; j++)
				if (get_brightness<uint>(img(i, j)) >= 250)
					tmp++;
			if (tmp > res) {
				res = tmp;
				res_num = i;
			}
		}
		return res_num;
	};

	auto clear_near_rows = [](Image &img, uint row_num) {
		for (uint j = 0; j < img.n_cols; j++)
			img(row_num, j) = std::make_tuple(0, 0, 0);
		if (row_num > 0)
			for (uint j = 0; j < img.n_cols; j++)
				img(row_num - 1, j) = std::make_tuple(0, 0, 0);
		if (row_num < img.n_rows - 1)
			for (uint j = 0; j < img.n_cols; j++)
				img(row_num + 1, j) = std::make_tuple(0, 0, 0);
	};

	auto find_maximum_col = [&](const Image &img, uint start, uint end) {
		uint res_num = start == 0 ? end : start, res = 0;
		for (uint i = start; i < end; i++) {
			uint tmp = 0;
			for (uint j = 0; j < img.n_rows; j++)
				if (get_brightness<uint>(img(j, i)) >= 250)
					tmp++;
			if (tmp > res) {
				res = tmp;
				res_num = i;
			}
		}
		return res_num;
	};

	auto clear_near_cols = [](Image &img, uint col_num) {
		for (uint j = 0; j < img.n_rows; j++)
			img(j, col_num) = std::make_tuple(0, 0, 0);
		if (col_num > 0)
			for (uint j = 0; j < img.n_rows; j++)
				img(j, col_num - 1) = std::make_tuple(0, 0, 0);
		if (col_num < img.n_cols - 1)
			for (uint j = 0; j < img.n_rows; j++)
				img(j, col_num + 1) = std::make_tuple(0, 0, 0);
	};

	// After finding the first maximum we must clear rows/cols near it to avoid getting
	// maximum near themselves. After that, find the second maximum and return that
	// is more near to centre of image.

	uint max1_num = find_maximum_row(src, 0, border_rows);
	clear_near_rows(src, max1_num);
	uint max2_num = find_maximum_row(src, 0, border_rows);
	res1 = std::max(max1_num, max2_num);

	max1_num = find_maximum_col(src, src.n_cols - border_cols, src.n_cols);
	clear_near_cols(src, max1_num);
	max2_num = find_maximum_col(src, src.n_cols - border_cols, src.n_cols);
	res2 = std::min(max1_num, max2_num);

	max1_num = find_maximum_row(src, src.n_rows - border_rows, src.n_rows);
	clear_near_rows(src, max1_num);
	max2_num = find_maximum_row(src, src.n_rows - border_rows, src.n_rows);
	res3 = std::min(max1_num, max2_num);

	max1_num = find_maximum_col(src, 0, border_cols);
	clear_near_cols(src, max1_num);
	max2_num = find_maximum_col(src, 0, border_cols);
	res4 = std::max(max1_num, max2_num);

	return std::make_tuple(res1, res2, res3, res4);
}

template<typename Metric>
std::pair<int, int> find_best_alignment(const Image &img1, const Image &img2, int shift_size, const Metric &op) {
	std::pair<int, int> res = std::make_pair(0, 0);
	if (shift_size > 1 &&
		(std::max(img1.n_rows, img1.n_cols) > 200 ||
		std::max(img2.n_rows, img2.n_cols) > 200)) {
		res = find_best_alignment(resize(img1, 0.5), resize(img2, 0.5), shift_size / 2, op);
		res = std::make_pair(res.first * 2, res.second * 2);
		shift_size = 1;
	}

	std::pair<int, int> best_pos = res, worst_pos = res;
	double best = LLONG_MIN, worst = LLONG_MAX;
	for (int i = -shift_size; i <= shift_size; i++)
		for (int j = -shift_size; j <= shift_size; j++) {
			double tmp = op(img1, img2, res.first + i, res.second + j);
			//std::cout << res.first + i << " " << res.second + j << "   " << tmp << std::endl;
			if (tmp > best) {
				best = tmp;
				best_pos = std::make_pair(res.first + i, res.second + j);
			}
			if (tmp < worst) {
				worst = tmp;
				worst_pos = std::make_pair(res.first + i, res.second + j);
			}
		}
	return op.param == 'B' ? best_pos : worst_pos;
}

std::pair<Image, Image> align_two_images(const Image &img1, const Image &img2, std::pair<int, int> shift) {
	auto align_direction = [](uint &i1_s, uint &i1_l, uint &i2_s, uint &i2_l, uint rows1, uint rows2, int sh) {
		if (sh < 0) {
			i1_s = -sh; i1_l = std::min(rows1 + sh, rows2);
			i2_s = 0; i2_l = i1_l;
		} else {
			i1_s = 0; i1_l = std::min(rows1, rows2 - sh);
			i2_s = sh; i2_l = i1_l;
		}
	};

	uint i1_x_s = 0, i1_y_s = 0, i1_x_l = 0, i1_y_l = 0;
	uint i2_x_s = 0, i2_y_s = 0, i2_x_l = 0, i2_y_l = 0;
	align_direction(i1_x_s, i1_x_l, i2_x_s, i2_x_l, img1.n_rows, img2.n_rows, shift.first);
	align_direction(i1_y_s, i1_y_l, i2_y_s, i2_y_l, img1.n_cols, img2.n_cols, shift.second);

	return std::make_pair(
		img1.submatrix(i1_x_s, i1_y_s, i1_x_l, i1_y_l), img2.submatrix(i2_x_s, i2_y_s, i2_x_l, i2_y_l));
}

Image align(const Image &src, const std::string &postprocessing, const double fraction) {
	uint height = src.n_rows / 3;

	Image B = src.submatrix(0, 0, height, src.n_cols);
	Image G = src.submatrix(height, 0, height, src.n_cols);
	Image R = src.submatrix(2 * height, 0, src.n_rows - 2 * height, src.n_cols);

	const double treshold1 = 150, treshold2 = 260;
	Image R1 = canny(R, treshold1, treshold2);
	Image G1 = canny(G, treshold1, treshold2);
	Image B1 = canny(B, treshold1, treshold2);

	const double near = 0.05;
	uint r_fr[4], g_fr[4], b_fr[4]; // up, right, down, left
	std::tie(r_fr[0], r_fr[1], r_fr[2], r_fr[3]) = get_frames_value(R1, near);
	std::tie(g_fr[0], g_fr[1], g_fr[2], g_fr[3]) = get_frames_value(G1, near);
	std::tie(b_fr[0], b_fr[1], b_fr[2], b_fr[3]) = get_frames_value(B1, near);

	Image Rsub = R.submatrix(r_fr[0], r_fr[3], r_fr[2] - r_fr[0], r_fr[1] - r_fr[3]);
	Image Gsub = G.submatrix(g_fr[0], g_fr[3], g_fr[2] - g_fr[0], g_fr[1] - g_fr[3]);
	Image Bsub = B.submatrix(b_fr[0], b_fr[3], b_fr[2] - b_fr[0], b_fr[1] - b_fr[3]);

	Image R2, G2, B2;
	if (postprocessing == "--subpixel") {
		R2 = resize(Rsub, fraction);
		G2 = resize(Gsub, fraction);
		B2 = resize(Bsub, fraction);
	} else {
		R2 = Rsub;
		G2 = Gsub;
		B2 = Bsub;
	}

	int shift_size1 = std::min(100, int(std::min(R2.n_rows, R2.n_cols)) / 3);
	std::pair<int, int> shift1 = find_best_alignment(R2, G2, shift_size1, MSE());
	//std::pair<int, int> shift1 = find_best_alignment(R2, G2, shift_size, Cross_Correlation());
	Image r_tmp, g_tmp;
	std::tie(r_tmp, g_tmp) = align_two_images(R2, G2, shift1);
	Image res1 = binary_map<std::tuple<uint, uint, uint>, R_plus_G>(R_plus_G(), r_tmp, g_tmp);

	int shift_size2 = std::min(100, int(std::min(res1.n_rows, res1.n_cols)) / 3);
	std::pair<int, int> shift2 = find_best_alignment(g_tmp, B2, shift_size2, MSE()); // Hack for very-very blue images
	//std::pair<int, int> shift2 = find_best_alignment(res1, B2, shift_size, Cross_Correlation());
	Image rg_tmp, b_tmp;
	std::tie(rg_tmp, b_tmp) = align_two_images(res1, B2, shift2);
	Image res = binary_map<std::tuple<uint, uint, uint>, RG_plus_B>(RG_plus_B(), rg_tmp, b_tmp);

	if (postprocessing == "--subpixel")
		return resize(res, 1.0 / fraction);
	if (postprocessing == "--gray_world")
		return gray_world(res);
	if (postprocessing == "--unsharp")
		return unsharp(res);
	if (postprocessing == "--autocontrast")
		return autocontrast(res, fraction);
	return res; // postprocessing == "--no" - no postprocessing
}
