#include "Utility.h"

using namespace std;

void LegendreBulk(double x, int n, double out0[], double out1[], double out2[]) {
	out0[0] = 1.0;
	out1[0] = 0.0;
	out2[0] = 0.0;
	if (n == 0) return;
	out0[1] = x;
	out1[1] = 1.0;
	out2[1] = 0.0;
	for (int currentN = 2; currentN <= n; currentN++) {
		double rcpN = 1.0 / currentN;
		double n2minus1 = (currentN << 1) - 1;
		int prevNInt = currentN - 1;
		int prev2NInt = prevNInt - 1;
		double prevNDouble = prevNInt;
		out0[currentN] = (n2minus1 * x * out0[prevNInt] - prevNDouble * out0[prev2NInt]) * rcpN;
		out1[currentN] = (n2minus1 * (out0[prevNInt] + x * out1[prevNInt]) - prevNDouble * out1[prev2NInt]) * rcpN;
		out2[currentN] = (n2minus1 * (out1[prevNInt] + (out1[prevNInt] + x * out2[prevNInt])) - prevNDouble * out2[prev2NInt]) * rcpN;
	}
}

LegendreDerivatives LegendreAll(double x, int n) {
	if (n == 0) {
		return LegendreDerivatives{ 1.0, 0.0, 0.0 };
	}
	else {
		double currentR2 = 1.0;
		double dcurrentR2dx = 0.0;
		double d2currentR2dx = 0.0;
		double currentR1 = x;
		double dcurrentR1dx = 1.0;
		double d2currentR1dx = 0.0;
		for (int currentN = 2; currentN <= n; currentN++) {
			double rcpN = 1.0 / currentN;
			double n2minus1 = 2 * currentN - 1;
			double prevN = currentN - 1;
			double currentR = (n2minus1 * x * currentR1 - prevN * currentR2) * rcpN;
			double dcurrentRdx = (n2minus1 * (currentR1 + x * dcurrentR1dx) - prevN * dcurrentR2dx) * rcpN;
			double d2currentRdx = (n2minus1 * (dcurrentR1dx + (dcurrentR1dx + x * d2currentR1dx)) - prevN * d2currentR2dx) * rcpN;
			currentR2 = currentR1;
			dcurrentR2dx = dcurrentR1dx;
			d2currentR2dx = d2currentR1dx;
			currentR1 = currentR;
			dcurrentR1dx = dcurrentRdx;
			d2currentR1dx = d2currentRdx;
		}
		return LegendreDerivatives{ currentR1, dcurrentR1dx, d2currentR1dx };
	}
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

glm::vec3 SpectrumToXYZ(float spectrum, float w) {
	float n = (w - 390.0);
	int i = int(n);
	if (i < 0 || i >= (830 - 390)) {
		return glm::vec3(0.0);
	}
	else {
		int n0 = glm::min(i - 1, int(n));
		int n1 = glm::min(i - 1, n0 + 1);
		glm::vec3 c0 = glm::vec3(cie[n0]);
		glm::vec3 c1 = glm::vec3(cie[n1]);

		glm::vec3 xyz = mix(c0, c1, n);

		xyz = xyz * spectrum;

		glm::vec3 cieIntegral = glm::vec3(0.0);
		for (int i = 0; i < 441; ++i) {
			cieIntegral += cie[i];
		}

		//Normalize the conversion
		return xyz * 441.0f / cieIntegral;
	}
}