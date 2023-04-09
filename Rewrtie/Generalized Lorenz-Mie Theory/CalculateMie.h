#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include <valarray>
#include <chrono>

#define _USE_MATH_DEFINES
#include <math.h>

struct ParticlePhase {
	double s_polarized;
	double p_polarized;
	double unpolarized;
	double cornetteShanks;
	double henyeyGreenstein;
	double g;
};

struct ParticleProperties {
	double extinction, scattering, absorption;
	std::complex<double> S1, S2;
};

struct MediumProperties {
	double extinction, scattering, absorption;

	ParticlePhase phase;
};

struct MediumParticles {
	std::complex<double> ior;
	double rMin, rMax;
	double rStep;
	std::valarray<double> N;
};

void ComputeParticleProperties(std::complex<double> hostIOR, std::complex<double> particleIOR, double theta, double r, double lambda, ParticlePhase& phase, ParticleProperties& particle);

void ComputeBulkMediumProperties(std::complex<double> hostIOR, double theta, double lambda, MediumParticles particles, MediumProperties& medium);