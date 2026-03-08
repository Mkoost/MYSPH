#include <cmath>

#ifndef MYSPH_KERNALS
#define MYSPH_KERNALS
namespace MYSPH {

	struct SplineKernal3D {

		static double W(double r, double h) {
			double R = r / h;
			double alpha = 3. / 2. * pi * h * h * h;
			if (R >= 2.0) return 0;
			if (R <= 1.0) return alpha * (1 - 3. * R * R / 2. + 3. * R * R * R / 4.);
			return alpha * (2. - R) * (2. - R) * (2. - R) / 4.;
		};

		static double dW(double r, double h) {
			double R = r / h;
			double alpha = 3. / 2. * pi * h * h * h;
			if (R >= 2.0) return 0;
			if (R <= 1.0) return alpha* (-3. * R + 9. * R * R / 4.) / h;
			return -3. * alpha * (2. - R) * (2. - R) / 4. / h;
		};

	};

	struct ForthOrderPolynomialKernal3D {
		static double W(double r, double h) {
			double R = r / h;
			double alpha = 315. / 208. * pi * h * h * h;

			if (R >= 2.0) return 0;
			return alpha * (-9. / 4. * R + 19. / 8 * R * R - 5. / 8. * R * R * R) / h;
		}

		static double dW(double r, double h) {
			double R = r / h;
			double alpha = 315. / 208. * pi * h * h * h;
			if (R >= 2.0) return 0;
			return alpha * (-9. / 4. * R + 19. / 8 * R * R - 5. / 8. * R * R * R) / h;
		}
	};

	struct ForthOrderPolynomialKernal1D {
		static double W(double r, double h) {
			double R = r / h;
			double alpha = 1. / h;
			if (R >= 2.0) return 0;
			return alpha * (2. / 3. - 9. / 8. * R * R + 19. / 24. * R * R * R - 5. / 32. * R * R * R * R);
		}

		static double dW(double r, double h) {
			double R = r / h;
			double alpha = 1./ h;
			if (R >= 2.0) return 0;
			return alpha * (-9. / 4. * R + 19. / 8 * R * R - 5. / 8. * R * R * R) / h;
		}
	};

	struct GaussKernal3D {
		static double W(double r, double h) {
			double R = r / h;
			double alpha = 1. / (std::pow(pi, 1.5) * h * h * h);
			return alpha * std::exp(- R * R);
		}

	    static double dW(double r, double h) {
			double R = r / h;
			double alpha = 1. / (std::pow(pi, 1.5) * h * h * h);
			return -2. * alpha * R * std::exp(- R * R);
		}
	};

}

#endif
