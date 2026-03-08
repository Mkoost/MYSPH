#pragma once

#include<cmath>
#include<algorithm>

#ifndef MYSPH_MATHEMATICS
#define MYSPH_MATHEMATICS


namespace MYSPH{
	constexpr double pi =  3.14159265358979323;
	constexpr double e  =  2.71828182845904523;


	template<class T>
	static double soundspeed(const T& pi) {
		return std::sqrt((5. / 3.) * pi.pressure / pi.density);
	};

	double sphere_metric(const dvec3& r1, const dvec3& r2) {
		double res = 0;
		for (size_t i = 0; i < 3; ++i)
			res += std::pow(r1[i] - r2[i], 2);
		return std::sqrt(res);
	}

	double sphere_norm(const dvec3& r1) {
		double res = 0;
		for (size_t i = 0; i < 3; ++i)
			res += r1[i] * r1[i];
		return std::sqrt(res);
	}

	double cube_norm(const dvec3& r1) {
		double res = 0;
		for (size_t i = 0; i < 3; ++i)
			res = std::max(std::abs(r1[i]), res);
		return (res);
	}

	double sphere_metric(const double& r1, const double& r2) {
		return std::fabs(r1 - r2);
	}

	double sphere_norm(const double& r1) {
		return std::fabs(r1);
	}

	double cube_norm(const double& r1) {
		return std::fabs(r1);
	}

	double binary_gravitational_potential(double r1, double r2, double proj1, double proj2, double mu = 1., double eps = 0) {
		r1 = std::max(r1 * r1 * r1, eps);
		r2 = std::max(r2 * r2 * r2, eps);
		return -(mu * proj1 / (r1) + (1. - mu) * proj2 / (r2));
	}

	// TODO: optimize f(x) calls number
	template<class Func>
	double secant_method_1d(Func f, double l, double r,  double step = 1e-7, double eps = 1e-10) {
		l += step;
		r -= step;

		double err = f(r);

		while (err > eps) {
			double xi = l - f(l) * (r - l) / (f(r) - f(l));
			err = std::fabs(f(xi));
			l = r;
			r = xi;
		}

		return r;

	}

	double find_lagrange_point(double mu, double padding = 1e-7, double eps = 1e-10) {
	    double x1 = -mu;
        double x2 = 1. - mu;
		return secant_method_1d([x1, x2, mu](double x) -> double {return MYSPH::binary_gravitational_potential(std::fabs(x1-x), std::fabs(x2-x), x - x1, x - x2, mu, 1e-10); },
								x2, x1, padding, eps);
	};


}
#endif
