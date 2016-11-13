#pragma once

#include "io.h"

Image custom(const Image &, const Matrix<double> &, const char = 'S', const bool = false);

Image sobel_x(const Image &, const char = 'S');
Image sobel_y(const Image &, const char = 'S');

Image unsharp(const Image &);

Image gray_world(const Image &);

Image resize(const Image &, const double);

Image autocontrast(const Image &, const double);

Image gaussian(const Image &, const double, const int);
Image gaussian_separable(const Image &, const double, const int);

Image median(const Image &, const int);
Image median_linear(const Image &, const int);

Image canny(const Image &, const double, const double);

Image align(const Image &, const std::string & = "--no", const double = 0.0);

#include "filters.hpp"