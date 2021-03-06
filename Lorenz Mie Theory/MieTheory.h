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
	double mean;
	double standardDeviation;
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
	double phase, phaseAsymmetry;
	double phaseCS, phaseHG;
};

unsigned int TermsToSum(const complex<double> z);

void ComputeParticleProperties(complex<double> iorHost, complex<double> iorParticle, double theta, double r, double lambda, complex<double>& S1, complex<double>& S2, double& Qabs, double& Qsca, double& Qext, double& phase, double& asymmetry);

void ComputeBulkOpticalProperties(complex<double> iorHost, double theta, double lambda, ParticleDistribution& particle, BulkMedium& bulk);