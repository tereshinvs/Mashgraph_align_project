#pragma once

#include "io.h"

const double PI = std::acos(-1.0);

double get_brightness(const std::tuple<uint, uint, uint> &);

bool is_correct_point(int, int, int, int);

double normalize_colour(double, char = 'S');

const std::tuple<uint, uint, uint> &get_reflected(const Image &, int, int);

uint get_median(const std::vector<int> &, const int);

class MSE {
	// MSE = 1 / (weight * height) * SUM (I1 - I2) ^ 2

	public:
		const char param = 'W';
		double operator()(const Image &img1, const Image &img2, int shift_x, int shift_y) const;
};

class Cross_Correlation {
	// CC = SUM I1 * I2
	public:
		const char param = 'B';
		double operator()(const Image &img1, const Image &img2, int shift_x, int shift_y) const;
};

class Gradient_Module_OP {
	public:
		double operator()(
			const Matrix<std::tuple<double, double, double>> &m1,
			const Matrix<std::tuple<double, double, double>> &m2) const;
		static const int radius1 = 0, radius2 = 0;
};

class Gradient_Direction_OP {
	public:
		double operator()(
			const Matrix<std::tuple<double, double, double>> &m1,
			const Matrix<std::tuple<double, double, double>> &m2) const;
		static const int radius1 = 0, radius2 = 0;
};

class Non_Maximum_Supression_OP {
	public:
		double operator()(const Matrix<double> &G, const Matrix<double> &O) const;
		static const int radius1 = 1, radius2 = 0;
	private:
		std::tuple<uint, uint> get_neighbour(const double o) const;
};

class Tresholding_OP {
	public:
		Tresholding_OP(const double _treshold1, const double _treshold2):
			treshold1(_treshold1), treshold2(_treshold2) {}
		uint operator()(const Matrix<double> &m) const;
		static const int radius = 0;
		const double treshold1, treshold2;
};

class R_plus_G {
	public:
		std::tuple<uint, uint, uint> operator()(const Image &m1, const Image &m2) const;
		static const int radius1 = 0, radius2 = 0;
};

class RG_plus_B {
	public:
		std::tuple<uint, uint, uint> operator()(const Image &m1, const Image &m2) const;
		static const int radius1 = 0, radius2 = 0;
};

class Make_RGB_OP {
	public:
		Make_RGB_OP(char _norm_type): norm_type(_norm_type) {}
		std::tuple<uint, uint, uint> operator()(const Matrix<std::tuple<double, double, double>> &m) const;
		static const int radius = 0;
	private:
		const char norm_type;
};

#include "image_tools.hpp"