#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include "glm.hpp"

using namespace std;

struct LogNormalParticleDistribution {
	double v;
	double c;
};

struct ParticleDistribution {
	complex<double> ior;
	double rMin, rMax;
	double stepSize;
	valarray<double> N;

	LogNormalParticleDistribution logNormalDist;
};

struct BulkMedium {
	double extinction, scattering, absorption;
	double phase;
};

void ComputeParticleProperties(complex<double> iorHost, complex<double> iorParticle, double theta, double r, double lambda, complex<double>& S1, complex<double>& S2, double& Qabs, double& Qsca, double& Qext, double& phase);

void ComputeBulkOpticalProperties(complex<double> iorHost, double theta, double lambda, ParticleDistribution& particle, BulkMedium& bulk);