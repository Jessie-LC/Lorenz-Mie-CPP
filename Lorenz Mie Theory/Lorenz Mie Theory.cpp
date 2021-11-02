#include "Utility.h"
#include "MieTheory.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 1800;

double MiePhase(double wl, complex<double> iorHost, complex<double> iorParticle, double radius) {
	//Based directly on https://github.com/nmoteki/gmie_cpp/blob/master/gmie.cpp
	float k0 = 2.0 * pi / wl;
	complex<double> k = iorHost * k;
	complex<double> x = k * radius;
	complex<double> iorRelative = iorParticle / iorHost;
	complex<double> y = x * iorRelative;

	complex<double> e_host = complex<double>(pow(iorHost.real(), 2.0) - pow(iorHost.imag(), 2.0), 2.0 * iorHost.real() * iorHost.imag());
	complex<double> e_particle = complex<double>(pow(iorParticle.real(), 2.0) - pow(iorParticle.imag(), 2.0), 2.0 * iorParticle.real() * iorParticle.imag());

	complex<double> u_host = complex<double>(1.0, 0.0);
	complex<double> u_particle = complex<double>(1.0, 0.0);

	int M = static_cast<unsigned int>(abs(x) + 4.3 * cbrt(abs(x)) + 2);

	vector<complex<double>> DD(M, { 0.0, 0.0 });
	for (unsigned int n = M - 1; n <= M; --n) {
		DD[n - 2] = (1.0 * n) / y - 1.0 / (DD[n - 1] + (1.0 * n) / y);
	}

	vector<complex<double>> RR(M + 1, { 0.0, 0.0 });
	RR[M] = x / (2.0 * M + 1.0);
	for (unsigned int n = M - 1; n > -1; --n) {
		RR[n] = 1.0 / ((2.0 * n + 1.0) / x - RR[n + 1]);
	}

	vector<complex<double>> psi(M + 1, { 0.0, 0.0 });
	psi[0] = RR[0] * cos(x);
	for (unsigned int n = 1; n <= M; ++n) {
		psi[n] = RR[n] * psi[n - 1];
	}

	vector<complex<double>> chi(M + 1, { 0.0, 0.0 });
	for (unsigned int n = 2; n <= M; ++n) {
		//This code probably does not work.
		if (n == 0) {
			chi[n] = -cos(x);
		} else if(n == 1) {
			chi[n] = (1.0 / x) * chi[n - 1] - sin(x);
		}
		chi[n] = ((2.0 * n - 1.0) / x) * chi[n - 1] - chi[n - 2];
	}

	vector<complex<double>> qsi(M + 1, { 0.0, 0.0 });
	for (unsigned int n = 0; n <= M; ++n) {
		qsi[n] = psi[n] + chi[n] * complex<double>(0.0, 1.0);
	}

	fill(RR.begin(), RR.end(), 0.0);
	RR[M] = y / (2.0 * M + 1.0);
	for (unsigned int n = M - 1; n > -1; --n) {
		RR[n] = 1.0 / ((2.0 * n + 1.0) / y - RR[n + 1]);
	}

	vector<complex<double>> psim(M + 1, { 0.0, 0.0 });
	psim[0] = RR[0] * cos(y);
	for (unsigned int n = 1; n <= M; ++n) {
		psim[n] = RR[n] * psim[n - 1];
	}

	vector<complex<double>> a(M, { 0.0, 0.0 }), b(M, { 0.0, 0.0 });
	for (unsigned int n = 0; n < M; ++n) {
		a[n] = ((e_particle / (e_host * iorRelative) * DD[n] + (n + 1.0) / x) * psi[n + 1] - psi[n]) / ((e_particle / (e_host * iorRelative) * DD[n] + (n + 1.0) / x) * qsi[n + 1] - qsi[n]); // BH83, Eq.4.88
		b[n] = ((e_host * iorRelative / e_particle * DD[n] + (n + 1.0) / x) * psi[n + 1] - psi[n]) / ((e_host * iorRelative / e_particle * DD[n] + (n + 1.0) / x) * qsi[n + 1] - qsi[n]); // BH83, Eq.4.88
	}

	vector<double> tmp(M, 0.0);
	for (unsigned int n = 1; n < M; ++n) {
		tmp[n] = (2.0 * n + 1.0) / (n * (n + 1.0));
	}

	vector<complex<double>> S1(angles, {0.0, 0.0}), S2(angles, {0.0, 0.0});
	vector<double> phase(angles, 0.0);
	for (int an = 0; an < angles; ++an) {
		double dtheta{ pi / angles };
		double theta = an * dtheta;
		double cosTheta = cos(theta);
		double sum = 0.0;
		for (unsigned int n = 1; n < M; ++n) {
			double Pi = computePi(cosTheta, n);
			double Tau = computeTau(cosTheta, n);
			S1[an] += tmp[n] * (a[n] * Pi + b[n] * Tau);
			S2[an] += tmp[n] * (a[n] * Tau + b[n] * Pi);
			sum += (2.0 * n + 1.0) * (pow(abs(a[n]), 2.0) + pow(abs(b[n]), 2.0));
		}
		phase[an] = (pow(abs(S1[an]), 2.0) + pow(abs(S2[an]), 2.0)) / ((4.0 * pi) * sum);
	}

	double Qsca{ 0.0 }, Qext{ 0.0 }, Qabs{ 0.0 };
	for (unsigned int n = 1; n < M; ++n) {
		Qsca += (2.0 * n + 1.0) * (pow(abs(a[n]), 2.0) + pow(abs(b[n]), 2.0));
		Qext += (2.0 * n + 1.0) * real((a[n] + b[n]) / (iorHost * iorHost));
	}
	double alpha = 4.0 * pi * radius * imag(iorHost) / wl;
	double _y = (2.0 * (1.0 + (alpha - 1.0) * exp(alpha))) / pow(alpha, 2.0);
	double term1 = pow(wl, 2.0) * exp(-(4.0 * pi * radius * imag(iorHost) / pow(wl, 2.0)));
	double term2 = 2.0 * pi * _y * pow(abs(iorHost), 2.0);

	Qsca = (term1 / term2) * Qsca;
	Qext = (pow(wl, 2.0) / tau) * Qext;

	Qabs = Qext - Qsca;

	return Qext;
}

/*
double ComputeMediumPhase(complex<double> iorHost, double theta, double lambda, ParticleDistribution& particle) {
	float phase = 0.0;
	int counter = 0;
	for (double r = particle.rMin + particle.stepSize * 0.5; r < particle.rMax; r += particle.stepSize) {
		complex<double> S1;
		complex<double> S2;
		double Qsca;
		double Qabs;
		double Qext;
		double particlePhase = ComputeParticlePhase(iorHost, particle.ior, theta, r, lambda, S1, S2, Qabs, Qsca, Qext);

		double sigmaS = Qsca * particle.N[counter] * particle.stepSize;

		double phaseI = Qsca * particlePhase * particle.stepSize;
		phaseI = (1.0 / sigmaS) * phaseI;

		particle.scattering += sigmaS;
		particle.absorption += Qabs * particle.N[counter] * particle.stepSize;
		particle.extinction += Qext * particle.N[counter] * particle.stepSize;
		phase += sigmaS * phaseI;

		++counter;
	}
	particle.scattering /= counter;
	phase /= particle.scattering;

	return phase;
}
*/

int main() {
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const vec3 mieData[] = vec3[](" << endl;

	std::cout << "Start generating Mie Phase:" << endl;

	double radius = 10.0e-6;

	complex<double> iorHost = { 1.34, 1.4e-9 };

	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.01e-6;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;

	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / 1e-6) + 1u);
	int counter = 0;
	for (double r = rMin_mineral + 1e-6 * 0.5; r < rMax_mineral; r += 1e-6) {
		N_mineral[counter] = 5.20019157e-8 * pow(2.0 * r, -3.4);
		++counter;
	}

	valarray<double> N_algae;
	counter = 0;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / 1e-6) + 1u);
	for (double r = rMin_algae + 1e-6 * 0.5; r < rMax_algae; r += 1e-6) {
		N_algae[counter] = 3.82828994e-7 * pow(2.0 * r, -3.6);
		++counter;
	}

	ParticleDistribution mineral{
		0.0,
		0.0,
		0.0,
		complex<double>{ 1.58, 3.61e-4 },
		rMin_mineral,
		rMax_mineral,
		1e-6,
		N_mineral,

		LogNormalParticleDistribution{
			5.429e-7,
			0.25
		}
	};

	ParticleDistribution algae{
		0.0,
		0.0,
		0.0,
		complex<double>{ 1.41, 8.19e-5 },
		rMin_algae,
		rMax_algae,
		1e-6,
		N_algae,

		LogNormalParticleDistribution{
			1.904e-6,
			0.15
		}
	};

	vector<complex<double>> S1, S2;
	vector<double> phase;
	double scatteringCoefficient, extinctionCoefficient, absorptionCoefficient;
	//tie(extinctionCoefficient, scatteringCoefficient, absorptionCoefficient, S1, S2) = MiePhase(550.0, iorHost, algae.ior, radius);

	for (int n = 0; n < angles; ++n) {
		double dtheta{ pi / angles };
		double theta = n * dtheta;
		BulkMedium bulk;
		ComputeMediumPhase(iorHost, theta, 550e-9, algae, bulk);
		std::cout << algae.scattering << endl;
		//std::cout << S1[n] << endl;
	}

	std::cout << "Finished generating Mie Phase" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}