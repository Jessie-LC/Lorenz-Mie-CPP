#include "Utility.h"

using namespace std;

double Pn(double x, int n) {
	if (n == 0) {
		return 1.0;
	}
	else {
		double currentR2 = 1.0;
		double currentR1 = x;
		for (int currentN = 2; currentN <= n; currentN++) {
			double currentR = ((2 * currentN - 1) * x * currentR1 - (currentN - 1) * currentR2) / currentN;
			currentR2 = currentR1;
			currentR1 = currentR;
		}
		return currentR1;
	}
}

double derivativeLegendre(int n, double x) {
	if (n == 0) {
		return 0.0;
	}
	else {
		double currentR2 = 1.0;
		double dcurrentR2dx = 0.0;
		double currentR1 = x;
		double dcurrentR1dx = 1.0;
		for (int currentN = 2; currentN <= n; currentN++) {
			double currentR = ((2 * currentN - 1) * x * currentR1 - (currentN - 1) * currentR2) / currentN;
			double dcurrentRdx = ((2 * currentN - 1) * (currentR1 + x * dcurrentR1dx) - (currentN - 1) * dcurrentR2dx) / currentN;
			currentR2 = currentR1;
			dcurrentR2dx = dcurrentR1dx;
			currentR1 = currentR;
			dcurrentR1dx = dcurrentRdx;
		}
		return dcurrentR1dx;
	}
}

double computePi(double mu, int n) {
	return derivativeLegendre(n, mu);
}

double derivativePi(double x, int n) {
	if (n == 0) {
		return 0.0;
	}
	else {
		double currentR2 = 1.0;
		double dcurrentR2dx = 0.0;
		double d2currentR2dx = 0.0;
		double currentR1 = x;
		double dcurrentR1dx = 1.0;
		double d2currentR1dx = 0.0;
		for (int currentN = 2; currentN <= n; currentN++) {
			double currentR = ((2 * currentN - 1) * x * currentR1 - (currentN - 1) * currentR2) / currentN;
			double dcurrentRdx = ((2 * currentN - 1) * (currentR1 + x * dcurrentR1dx) - (currentN - 1) * dcurrentR2dx) / currentN;
			double d2currentRdx = ((2 * currentN - 1) * (dcurrentR1dx + (dcurrentR1dx + x * d2currentR1dx)) - (currentN - 1) * d2currentR2dx) / currentN;
			currentR2 = currentR1;
			dcurrentR2dx = dcurrentR1dx;
			d2currentR2dx = d2currentR1dx;
			currentR1 = currentR;
			dcurrentR1dx = dcurrentRdx;
			d2currentR1dx = d2currentRdx;
		}
		return d2currentR1dx;
	}
}

double computeTau(double mu, int n) {
	double mu_sin = sin(acos(mu));
	return mu * computePi(mu, n) - pow(mu_sin, 2.0) * derivativePi(mu, n);
}

//Everything past here can be optimized most likely.
double Factorial(int n) {
	double x = 1.0;
	for (int i = 1; i <= n; ++i) {
		x *= i;
	}
	return x;
}

double T_Gamma(double z) {
	//This is the Gamma function.
	double product = 1.0 / z;
	int n = 1;
	while (true) {
		double value = pow(1.0 + (1.0 / n), z) / (1.0 + (z / (double)n));
		if (max(abs(value), 1.0 / abs(value)) < 1.000000001) {
			break;
		}
		product *= value;
		++n;
	}
	return product;
}

complex<double> Jn(double n, complex<double> z) {
	//This is the Bessel function of the first kind.
	complex<double> sum = 0.0;
	for (int m = 0; m < 20; ++m) {
		//Apparently the gamma function is equal to the factorial if you do "n + 1" as the input.
		sum += (((m & 1) ? -1.0 : 1.0) * pow(z / 2.0, 2 * m + n)) / (tgamma(m + 1) * tgamma(m + n + 1));
	}
	return sum;
}

complex<double> SphJn(int n, complex<double> z) {
	//This is the Spherical Bessel function of the first kind.
	if (isnan(z.real()) || isnan(z.imag())) {
		return z;
	}
	if (isinf(z.real()) || -isinf(z.imag())) {
		return complex<double>(0.0, 0.0);
	}
	if (z.real() == 0.0 && z.imag() == 0.0) {
		if (n == 0) {
			return complex<double>(1.0, 0.0);
		}
		else {
			return complex<double>(0.0, 0.0);
		}
	}
	complex<double> j = sqrt(pi / (2.0 * z)) * Jn((double)n + 0.5, z);
	if (isnan(j.real()) || isnan(j.imag())) {
		//If the bessel function produces a NaN, return 1.0+0.0i.
		return complex<double>(1.0, 0.0);
	}
	return j;
}

complex<double> Yn(double n, complex<double> z) {
	//This is the Bessel function of the second kind.
	return (Jn(n, z) * cos(n * pi) - Jn(-n, z)) / sin(n * pi);
}

complex<double> SphYn(int n, complex<double> z) {
	//This is the Spherical Bessel function of the second kind.
	if (isnan(z.real()) || isnan(z.imag())) {
		return z;
	}
	if (isinf(z.real()) || -isinf(z.imag())) {
		return complex<double>(0.0, 0.0);
	}
	if (z.real() == 0.0 && z.imag() == 0.0) {
		if (n == 0) {
			return complex<double>(1.0, 0.0);
		}
		else {
			return complex<double>(0.0, 0.0);
		}
	}
	complex<double> y = sqrt(pi / (2.0 * z)) * Yn((double)n + 0.5, z);
	if (isnan(y.real()) || isnan(y.real())) {
		return complex<double>(1.0, 0.0);
	}
	return y;
}

double sqr(double x) {
	return x * x;
}

glm::dvec3 xyzToRGB(glm::dvec3 xyz) {
	return xyz * xyzToRGBMatrix;
}

glm::dvec3 SpectrumToXYZ(glm::dvec3 spectrum, float w) {
	float n = (w - 390.0);
	int i = int(n);
	if (i < 0 || i >= (830 - 390)) {
		return glm::dvec3(0.0);
	}
	else {
		int n0 = min(i - 1, int(n));
		int n1 = min(i - 1, n0 + 1);
		glm::dvec3 c0 = cie[n0];
		glm::dvec3 c1 = cie[n1];

		glm::dvec3 xyz = mix(c0, c1, n);

		xyz = xyz * spectrum;

		//Normalize the conversion
		return xyz * 441.0 / glm::dvec3(113.042);
	}
}