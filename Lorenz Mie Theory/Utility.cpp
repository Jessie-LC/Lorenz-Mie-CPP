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

double sqr(double x) {
	return x * x;
}

glm::dvec3 xyzToRGB(glm::dvec3 xyz) {
	return xyz * xyzToRGBMatrix;
}

glm::dvec3 spectrumToXYZ(glm::dvec3 spectrum, float w) {
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