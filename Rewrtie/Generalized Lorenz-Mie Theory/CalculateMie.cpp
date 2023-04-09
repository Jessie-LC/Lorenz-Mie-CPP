#include "CalculateMie.h"

struct LegendreDerivatives {
	double value;
	double derivative1;
	double derivative2;
};

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

void ComputeParticleProperties(std::complex<double> hostIOR, std::complex<double> particleIOR, double theta, double radius, double lambda, ParticlePhase& phase, ParticleProperties& particle) {
	double x = 2.0 * 3.14159 * radius / lambda;
	std::complex<double> k = 2.0 * 3.14159 * hostIOR / lambda;

	std::complex<double> hostZ = x * hostIOR;
	std::complex<double> particleZ = x * particleIOR;
	double size = abs(hostZ);
	unsigned int MAX_N = static_cast<unsigned int>(ceil(size + 4.3 * cbrt(size) + 1.0));
		
	std::valarray<std::complex<double>> hostA;
	std::valarray<std::complex<double>> particleA;
	hostA.resize(MAX_N + 2u);
	particleA.resize(MAX_N + 2u);
	hostA[MAX_N + 1] = std::complex<double>(1.0, 0.0);
	particleA[MAX_N + 1] = std::complex<double>(1.0, 0.0);
	for (unsigned int n = MAX_N; n <= MAX_N; --n) {
		std::complex<double> tmpHost = std::complex<double>(n + 1.0, 0.0) / hostZ;
		std::complex<double> tmpParticle = std::complex<double>(n + 1.0, 0.0) / particleZ;

		hostA[n] = tmpHost - (std::complex<double>(1.0, 0.0) / (tmpHost + hostA[n + 1]));
		particleA[n] = tmpParticle - (std::complex<double>(1.0, 0.0) / (tmpParticle + particleA[n + 1]));
	}

	std::complex<double> B = std::complex<double>(0.0, 1.0);
	std::complex<double> psiZeta = std::complex<double>(0.5, 0.0) * (1.0 - exp(std::complex<double>(-2.0 * hostZ.imag(), 2.0 * hostZ.real())));
	std::complex<double> ratio   = std::complex<double>(0.5, 0.0) * (1.0 - exp(std::complex<double>(2.0 * hostZ.imag(), -2.0 * hostZ.real())));
	std::valarray<std::complex<double>> a;
	std::valarray<std::complex<double>> b;
	a.resize(MAX_N + 1u);
	b.resize(MAX_N + 1u);
	for (unsigned int n = 1u; n < MAX_N; ++n) {
		std::complex<double> n_z = std::complex<double>(n, 0.0) / hostZ;
		psiZeta *= (n_z - hostA[n - 1]) * (n_z - B);
		B = hostA[n] + (std::complex<double>(0.0, 1.0) / psiZeta);
		ratio *= (B + n_z) / (hostA[n] + n_z);

		a[n] = ratio * ((hostIOR * particleA[n] - particleIOR * hostA[n]) / (hostIOR * particleA[n] - particleIOR * B));
		b[n] = ratio * ((particleIOR * particleA[n] - hostIOR * hostA[n]) / (particleIOR * particleA[n] - hostIOR * B));
	}

	double cosTheta = cos(theta);
	double sinTheta = sin(theta);

	particle.S1 = std::complex<double>(0.0);
	particle.S2 = std::complex<double>(0.0);
	particle.scattering = 0.0;
	std::complex<double> ext = std::complex<double>(0.0);
	double sum = 0.0;
	double g_sum = 0.0;
	for (unsigned int n = 1u; n < MAX_N; ++n) {
		LegendreDerivatives legendrePolynomial = LegendreAll(cosTheta, n);
		double Pi = legendrePolynomial.derivative1;
		double Tau = cosTheta * legendrePolynomial.derivative1 - pow(sinTheta, 2.0) * legendrePolynomial.derivative2;

		double tmp = (2.0 * n + 1.0) / (n * (n + 1.0));
		particle.S1 += tmp * (a[n] * Pi + b[n] * Tau);
		particle.S2 += tmp * (a[n] * Tau + b[n] * Pi);

		particle.scattering += (2.0 * n + 1.0) * (pow(abs(a[n]), 2.0) + pow(abs(b[n]), 2.0));
		ext += (2.0 * n + 1.0) * (a[n] + b[n]) / (hostIOR * hostIOR);

		g_sum += ((n * (n + 2)) / (n + 1)) * real(a[n] * conj(a[n + 1]) + b[n] * conj(b[n + 1])) + tmp * real(a[n] * conj(b[n]));
	}
	double alpha = 4.0 * 3.14159 * radius * imag(hostIOR) / lambda;
	double gamma = alpha < 1e-6 ? 1.0 : (2.0 * (1.0 + (alpha - 1.0) * exp(alpha))) / pow(alpha, 2.0);
	double term1 = pow(lambda, 2.0) * exp(-alpha);
	double term2 = 2.0 * 3.14159 * gamma * pow(abs(hostIOR), 2.0);

	particle.scattering *= pow(lambda, 2.0) / (2.0 * M_PI) * exp(-alpha) * (alpha > 1e-6 ? pow(alpha, 2.0) / (2.0 * (1.0 + (alpha - 1.0) * exp(alpha))) : 1.0);
	particle.extinction = (pow(lambda, 2.0) / (2.0 * M_PI)) * ext.real();
	particle.absorption = particle.extinction - particle.scattering;

	phase.unpolarized = (pow(abs(particle.S1), 2.0) + pow(abs(particle.S2), 2.0)) / (2.0 * pow(abs(k), 2.0) * particle.scattering);
	phase.s_polarized = pow(abs(particle.S1), 2.0) / (pow(abs(k), 2.0) * particle.scattering);
	phase.p_polarized = pow(abs(particle.S2), 2.0) / (pow(abs(k), 2.0) * particle.scattering);

	if (std::isnan(phase.unpolarized)) {
		phase.unpolarized = 0.0;
		phase.s_polarized = 0.0;
		phase.p_polarized = 0.0;
	}
	if (std::isnan(particle.scattering)) {
		particle.scattering = 0.0;
	}
	if (std::isnan(particle.extinction)) {
		particle.extinction = 0.0;
	}

	phase.g = g_sum / (0.5 * sum);
}

void ComputeBulkMediumProperties(std::complex<double> hostIOR, double theta, double lambda, MediumParticles mParticle, MediumProperties& medium) {
	/*
		iParticle = individual particle
		mParticle = medium/multiple particle(s)
	*/
	medium.absorption = 0.0;
	medium.scattering = 0.0;
	medium.extinction = 0.0;

	medium.phase.unpolarized = 0.0;
	medium.phase.s_polarized = 0.0;
	medium.phase.p_polarized = 0.0;
	medium.phase.g = 0.0;

	int counter = 0;
	for (double r = mParticle.rMin + mParticle.rStep * 0.5; r < mParticle.rMax; r += mParticle.rStep) {
		ParticlePhase phase;
		ParticleProperties iParticle;
		ComputeParticleProperties(hostIOR, mParticle.ior, theta, r, lambda, phase, iParticle);

		double sigmaS = iParticle.scattering * mParticle.N[counter] * mParticle.rStep;

		medium.scattering += sigmaS;
		medium.extinction += iParticle.extinction * mParticle.N[counter] * mParticle.rStep;
		medium.absorption += iParticle.absorption * mParticle.N[counter] * mParticle.rStep;

		medium.phase.unpolarized += sigmaS * phase.unpolarized;
		medium.phase.s_polarized += sigmaS * phase.s_polarized;
		medium.phase.p_polarized += sigmaS * phase.p_polarized;
		medium.phase.g += sigmaS * phase.g;

		++counter;
	}
	medium.phase.unpolarized = (1.0 / medium.scattering) * medium.phase.unpolarized;
	medium.phase.s_polarized = (1.0 / medium.scattering) * medium.phase.s_polarized;
	medium.phase.p_polarized = (1.0 / medium.scattering) * medium.phase.p_polarized;
	medium.phase.g = (1.0 / medium.scattering) * medium.phase.g;
}