// Generalized Lorenz-Mie Theory.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "CalculateMie.h"
#include "RefractiveIndex.h"
#include <glm.hpp>
#include "Constants.h"

const int angles = 180;

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

		glm::vec3 xyz = glm::mix(c0, c1, n);

		xyz = xyz * spectrum;

		glm::vec3 cieIntegral = glm::vec3(0.0);
		for (int i = 0; i < 441; ++i) {
			cieIntegral += cie[i];
		}

		//Normalize the conversion
		return xyz * 441.0f / cieIntegral;
	}
}

int main()
{
	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.01e-6;
	double mineral_stepSize = 0.05e-6;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;
	double algae_stepSize = 0.1e-6;

	double vAlgae = 2.171e-6;
	double vMineral = 2.077e-6;

	int counter = 0;
	std::valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);
	for (double r = rMin_mineral + mineral_stepSize * 0.5; r < rMax_mineral; r += mineral_stepSize) {
		N_mineral[counter] = (vMineral / 10.44) * pow(2.0 * r, -3.4);
		++counter;
	}

	counter = 0;
	std::valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);
	for (double r = rMin_algae + algae_stepSize * 0.5; r < rMax_algae; r += algae_stepSize) {
		N_algae[counter] = (vAlgae / 4.97) * pow(2.0 * r, -3.6);
		++counter;
	}

	MediumParticles mineral{
		std::complex<double>{ 1.58, 3.61e-4 },
		rMin_mineral,
		rMax_mineral,
		mineral_stepSize,
		N_mineral
	};

	MediumParticles algae{
		std::complex<double>{ 1.41, 8.19e-5 },
		rMin_algae,
		rMax_algae,
		algae_stepSize,
		N_algae
	};

	MediumParticles particles[2];
	particles[0] = algae;
	particles[1] = mineral;

	const int wavelengths = 88;
	glm::vec3 rgbScattering = glm::vec3(0.0f);
	glm::vec3 rgbExtinction = glm::vec3(0.0f);
	for (int i = 0; i < wavelengths; ++i) {
		double r = (double)i / wavelengths;
		double lambda = (390.0 + (441.0 * r)) * 1e-9;

		std::cout << lambda << ", " << i << std::endl;

		particles[0].ior = std::complex<double>(particles[0].ior.real(), InterpAlgaeK(lambda * 1e9));
		particles[1].ior = std::complex<double>(particles[1].ior.real(), InterpMineralK(lambda * 1e9));

		double absorptionMedium = 4.0 * M_PI * WaterIOR(lambda * 1e9).imag() / lambda;

		double scattering = 0.0;
		double extinction = 0.0;
		double absorption = 0.0;
		for (int n = 0; n < 2; ++n) {
			MediumProperties medium;
			ComputeBulkMediumProperties(WaterIOR(lambda * 1e9), 0.0, lambda, particles[n], medium);

			scattering += medium.scattering;
			extinction += medium.extinction;
		}
		extinction = absorptionMedium + extinction;
		absorption = extinction - scattering;

		rgbScattering += SpectrumToXYZ(scattering, lambda * 1e9) / (float)wavelengths;
		rgbExtinction += SpectrumToXYZ(extinction, lambda * 1e9) / (float)wavelengths;
	}

	std::cout << "vec3(" << rgbScattering.r << ", " << rgbScattering.g << ", " << rgbScattering.b << ")," << std::endl;
	std::cout << "vec3(" << rgbExtinction.r << ", " << rgbExtinction.g << ", " << rgbExtinction.b << ")," << std::endl;

	/*
	for (int i = 0; i < angles; ++i) {
		double dtheta{ 3.14159 / (angles - 1) };
		double theta = i * dtheta;

		double phase = 0.0;
		double scattering = 0.0;
		double extinction = 0.0;
		double absorption = 0.0;
		for (int n = 0; n < 2; ++n) {
			MediumProperties medium;
			ComputeBulkMediumProperties(std::complex<double>(1.333, 0.0), theta, 550e-9, particles[n], medium);

			scattering += medium.scattering;
			extinction += medium.extinction;
			absorption += medium.absorption;

			phase += medium.phase.unpolarized * medium.scattering;
		}
		phase = (1.0 / scattering) * phase;

		extinction = scattering + absorption;

		std::cout << phase << ", " << scattering << ", " << extinction << std::endl;
	}
	*/
}