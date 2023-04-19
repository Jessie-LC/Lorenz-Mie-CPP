// Generalized Lorenz-Mie Theory.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "CalculateMie.h"
#include "RefractiveIndex.h"
#include <glm.hpp>
#include "Constants.h"

//#define OCEAN
#define CALCULATE_PHASE_FUNCTION

const int angles = 3600;
#if defined OCEAN && defined CALCULATE_PHASE_FUNCTION
	const int wavelengths = 1;
#else
	const int wavelengths = 100;
#endif
const int threadCount = 10;

glm::vec3 SpectrumToXYZ(float spectrum, float w) {
	float n = w - 390.0f;
	int n0 = (int)n;
	if (n0 < 0 || n >= 830.0f - 390.0f) {
		return glm::vec3(0.0);
	}

	int n1 = glm::min(n0 + 1, 440);

	glm::vec3 xyz = glm::mix(cie[n0], cie[n1], glm::mod(n, 1.0f));

	xyz = xyz * spectrum;

	return xyz * 441.0f / 113.042f;
}

float Plancks(float t, float lambda) {
	const float h = 6.62607015e-16;
	const float c = 2.9e17;
	const float k = 1.38e-5;
	const float c_1L = 2.0 * h * c * c;
	const float c_2 = h * c / k;

	float p1 = 2.0 * h * pow(c, 2.0) * pow(lambda, -5.0);
	float p2 = exp((h * c) / (lambda * k * t)) - 1.0;

	return p1 / p2;
}

float maxof(glm::vec3 p) {
	return glm::max(p.r, glm::max(p.g, p.b));
}

void CalculateMediumCoefficients(int s, glm::vec3 *xyzCoefficients, MediumParticles particle0, MediumParticles particle1, int particleTypes) {
	glm::vec3 D65 = glm::vec3(0.0f);
	for (int i = s; i < wavelengths; i += threadCount) {
		double r = (double)i / wavelengths;
		double lambda = (390.0 + (441.0 * r)) * 1e-9;

#ifdef OCEAN
		MediumParticles particles[2];
		particles[0] = particle0;
		particles[1] = particle1;
		std::complex<double> hostIOR = BrineIOR(lambda * 1e9);
		particles[0].ior = std::complex<double>(particles[0].ior.real(), InterpAlgaeK(lambda * 1e9));
		particles[1].ior = std::complex<double>(particles[1].ior.real(), InterpMineralK(lambda * 1e9));
#else
		MediumParticles particles = particle0;
		std::complex<double> hostIOR = AirIOR(lambda * 1e9);
		particles.ior = WaterIOR(lambda * 1e9);
#endif

		double absorptionMedium = 4.0 * M_PI * imag(hostIOR) / lambda;

		double scattering = 0.0;
		double extinction = 0.0;
		double absorption = 0.0;
		double g = 0.0;
		for (int n = 0; n < particleTypes; ++n) {
#ifdef OCEAN
			MediumParticles particle = particles[n];
#else
			MediumParticles particle = particles;
#endif
			MediumProperties medium;
			ComputeBulkMediumProperties(hostIOR, 0.0, lambda, particle, medium);

			scattering += medium.scattering;
			extinction += medium.extinction;

			g += medium.scattering * medium.phase.g;
		}
		double pExtinction = extinction;
		extinction = absorptionMedium + extinction;
		absorption = extinction - scattering;
		g = (1.0 / scattering) * g;

		if (std::isnan(g)) {
			g = 0.0;
		}
		if (std::isinf(g)) {
			g = 1.0;
		}

		float D65_spectrum = Plancks(6504.0, lambda * 1e9);

		xyzCoefficients[0] += SpectrumToXYZ(scattering * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[1] += SpectrumToXYZ(extinction * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[2] += SpectrumToXYZ(absorption * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[3] += SpectrumToXYZ(absorptionMedium * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[4] += SpectrumToXYZ(pExtinction * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[5] += SpectrumToXYZ((pExtinction - scattering) * D65_spectrum, lambda * 1e9) / (float)wavelengths;
		xyzCoefficients[6] += SpectrumToXYZ(D65_spectrum, lambda * 1e9) / (float)wavelengths;
	}
}

void CalculatePhaseFunction(int s, float* phaseFunction, MediumParticles particle0, MediumParticles particle1, int particleTypes) {
	for (int i = 0; i < wavelengths; ++i) {
		double r = (double)i / wavelengths;
		double lambda = (390.0 + (441.0 * r)) * 1e-9;
		for (int n = s; n < angles; n += threadCount) {
			double dtheta{ M_PI / (angles - 1) };
			double theta = n * dtheta;

	#ifdef OCEAN
			MediumParticles particles[2];
			particles[0] = particle0;
			particles[1] = particle1;
			std::complex<double> hostIOR = BrineIOR(lambda * 1e9);
			particles[0].ior = std::complex<double>(particles[0].ior.real(), InterpAlgaeK(lambda * 1e9));
			particles[1].ior = std::complex<double>(particles[1].ior.real(), InterpMineralK(lambda * 1e9));
	#else
			MediumParticles particles = particle0;
			std::complex<double> hostIOR = AirIOR(lambda * 1e9);
			particles.ior = WaterIOR(lambda * 1e9);
	#endif

			double absorptionMedium = 4.0 * M_PI * imag(hostIOR) / lambda;

			double scattering = 0.0;
			double extinction = 0.0;
			double absorption = 0.0;
			double phase = 0.0;
			double g = 0.0;
			for (int m = 0; m < particleTypes; ++m) {
	#ifdef OCEAN
				MediumParticles particle = particles[m];
	#else
				MediumParticles particle = particles;
	#endif
				MediumProperties medium;
				ComputeBulkMediumProperties(hostIOR, theta, lambda, particle, medium);

				scattering += medium.scattering;
				extinction += medium.extinction;

				phase += medium.scattering * medium.phase.unpolarized;
				g += medium.scattering * medium.phase.g;
			}
			extinction = absorptionMedium + extinction;
			absorption = extinction - scattering;
			phase = (1.0 / scattering) * phase;
			g = (1.0 / scattering) * g;
			if (std::isnan(g)) {
				g = 0.0;
			}
			if (std::isinf(g)) {
				g = 1.0;
			}

			phaseFunction[n + i * angles] = phase;
		}
	}
}

int main()
{
#ifdef OCEAN
	const int particleTypes = 2;
	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.01e-6;
	double mineral_stepSize = rMin_mineral * 1.0;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.2e-6;
	double algae_stepSize = rMin_algae;

	const int waterBody = 4;
	double vAlgae = 2.171e-6;
	double vMineral = 2.077e-6;
	switch (waterBody) {
	case 0: {
		break;
	}
	case 1: {
		vAlgae = 3.878e-7;
		vMineral = 3.075e-7;
		break;
	}
	case 2: {
		vAlgae = 4.999e-7;
		vMineral = 2.300e-7;
		break;
	}
	case 3: {
		vAlgae = 1.904e-6;
		vMineral = 5.429e-7;
		break;
	}
	case 4: {
		vAlgae = 1.880e-7;
		vMineral = 2.477e-10;
		break;
	}
	}

	std::cout << "Mineral Particle Count: " << ((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u << std::endl;
	std::cout << "Algae Particle Count: " << ((rMax_algae - rMin_algae) / algae_stepSize) + 1u << std::endl;

	double volume = 0.0;
	int counter = 0;
	std::valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);
	for (double r = rMin_mineral + mineral_stepSize * 0.5; r < rMax_mineral; r += mineral_stepSize) {
		N_mineral[counter] = (vMineral / 10.44) * pow(2.0 * r, -3.4);
		volume += pow(r, 3.0) * N_mineral[counter] * mineral_stepSize;
		++counter;
	}
	double mineralVolume = vMineral / (((4.0 * M_PI) / 3.0) * volume);
	N_mineral *= mineralVolume;

	volume = 0.0;
	counter = 0;
	std::valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);
	for (double r = rMin_algae + algae_stepSize * 0.5; r < rMax_algae; r += algae_stepSize) {
		N_algae[counter] = (vAlgae / 4.97) * pow(2.0 * r, -3.6);
		volume += pow(r, 3.0) * N_algae[counter] * algae_stepSize;
		++counter;
	}
	double algaeVolume = vAlgae / (((4.0 * M_PI) / 3.0) * volume);
	N_algae *= algaeVolume;

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
#else
	const int aerosolType = 1;
	double mean = 8.0e-7;
	double standardDeviation = mean * 1.5;

	const int particleTypes = 1;
	double rMax_water = mean * 10.0 * 2.0;
	double rMin_water = mean / 10.0;
	double water_stepSize = rMin_water * 1.0;
	switch (aerosolType) {
		case 0: {
			break;
		}
		case 1: {
			rMax_water = 4e-4;
			rMin_water = 1e-4;
			water_stepSize = rMin_water * 1.0;
			mean = rMax_water - rMin_water;
			standardDeviation = mean * 20.0;
			break;
		}
		case 2: {
			rMax_water = 5e-8;
			rMin_water = 5e-9;
			water_stepSize = rMin_water * 1.0;
			mean = rMax_water - rMin_water;
			standardDeviation = mean * 20.0;
			break;
		}
	}

	std::valarray<double> N_water;

	double beta_sqr = log(((standardDeviation * standardDeviation) / (mean * mean)) + 1.0);
	double alpha = log(mean) - 0.5 * beta_sqr;
	double a = 4.9757e6;
	double beta = sqrt(beta_sqr);

	double volume = 0.0;

	std::cout << "Particle Count: " << ((rMax_water - rMin_water) / water_stepSize) + 1u << std::endl;

	int counter = 0;
	N_water.resize(static_cast<int>((rMax_water - rMin_water) / water_stepSize) + 1);
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		double x = r * 1e6;
		double logNormal = (1.0 / (r * beta * sqrt(M_PI * 2.0))) * exp(-0.5 * pow((log(r) - alpha) / beta, 2.0));
		double normal = (1.0 / (standardDeviation * sqrt(M_PI * 2.0))) * exp(-0.5 * pow((r - mean) / standardDeviation, 2.0));
		double rainbow = pow(2.0 * x, -3.5);
		double cumulus = 2.373 * pow(x, 6.0) * exp(-1.5 * x);
		double haze = 5.33e4 * pow(x, 4.0) * exp(-8.9 * pow(x, 0.5));
		N_water[counter] = cumulus;
		volume += pow(r, 3.0) * N_water[counter] * water_stepSize;
		++counter;
	}

	double waterWeight = 0.0005;
	double waterDensity = 1000.0;
	double waterVolume = waterWeight < 1e-12 ? 0.0 : waterWeight / waterDensity;
	double airVolume = 1.0 - waterVolume;

	double N = ((waterWeight > 1e-12 ? waterVolume : 1.0) / airVolume) / (4.1887902 * volume);

	counter = 0;
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		N_water[counter] = N * N_water[counter];
		++counter;
	}

	MediumParticles particles{
		std::complex<double>{ 0.0, 0.0 },
		rMin_water,
		rMax_water,
		water_stepSize,
		N_water
	};
#endif

#ifdef CALCULATE_PHASE_FUNCTION
		float phaseFunction[angles * wavelengths];

	#ifdef OCEAN
		MediumParticles particle0 = particles[0];
		MediumParticles particle1 = particles[1];
	#else
		MediumParticles particle0 = particles;
		MediumParticles particle1 = particles;
	#endif
		auto start = std::chrono::high_resolution_clock::now();

		std::valarray<std::thread> th;
		th.resize(threadCount);
		for (int i = 0; i < threadCount; ++i) {
			th[i] = std::thread(CalculatePhaseFunction, i, phaseFunction, particle0, particle1, particleTypes);
		}

		for (int i = 0; i < threadCount; ++i) {
			th[i].join();
		}

		std::valarray<double> anglesRadians;
		anglesRadians.resize(angles);
		for (int n = 0; n < angles; ++n) {
			double dtheta{ M_PI / (angles - 1) };
			double theta = n * dtheta;

			anglesRadians[n] = theta;
		}

		float CDF[angles * wavelengths];
		for (int i = 0; i < wavelengths; ++i) {
			CDF[0 + i * angles] = 0.0f;
			for (int n = 1; n < angles; ++n) {
				float angleStart = anglesRadians[n - 1];
				float angleEnd = anglesRadians[n];
				float phaseStart = phaseFunction[(n - 1) + i * angles];
				float phaseEnd = phaseFunction[n + i * angles];

				float p1 = phaseStart * (cos(angleStart) - cos(angleEnd));
				float p2 = (phaseEnd - phaseStart) / (angleEnd - angleStart);
				float p3 = (sin(angleEnd) - angleEnd * cos(angleEnd)) - (sin(angleStart) - angleStart * cos(angleStart));
				float p4 = angleStart * (cos(angleStart) - cos(angleEnd));

				float integral = (float)2.0f * (float)M_PI * (p1 + p2 * (p3 - p4));
				CDF[n + i * angles] = CDF[(n - 1) + i * angles] + integral;
			}
		}

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::nanoseconds time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

		for (int i = 0; i < wavelengths; ++i) {
			float maxCDF = CDF[(angles - 1) + i * angles];

			for (int n = 0; n < angles; ++n) {
				phaseFunction[n + i * angles] /= maxCDF;
				CDF[n + i * angles] /= maxCDF;
			}
		}

		std::cout << time.count() << std::endl;

		std::ofstream phaseLut("phase.dat", std::ios::binary);
		phaseLut.write(reinterpret_cast<char*>(phaseFunction), sizeof(float) * angles * wavelengths);
		phaseLut.close();
		std::cout << "Finished writing phase!" << std::endl;
		std::ofstream cdfLut("cdf.dat", std::ios::binary);
		cdfLut.write(reinterpret_cast<char*>(CDF), sizeof(float) * angles * wavelengths);
		cdfLut.close();
		std::cout << "Finished writing CDF!" << std::endl;
	#else
		#ifdef OCEAN
			MediumParticles particle0 = particles[0];
			MediumParticles particle1 = particles[1];
		#else
			MediumParticles particle0 = particles;
			MediumParticles particle1 = particles;
		#endif
		auto start = std::chrono::high_resolution_clock::now();

		glm::vec3 xyzCoefficients[7];
		for (int i = 0; i < 7; ++i) {
			xyzCoefficients[i] = glm::vec3(0.0f);
		}

		std::valarray<std::thread> th;
		th.resize(threadCount);
		for (int i = 0; i < threadCount; ++i) {
			th[i] = std::thread(CalculateMediumCoefficients, i, xyzCoefficients, particle0, particle1, particleTypes);
		}

		for (int i = 0; i < threadCount; ++i) {
			th[i].join();
		}
		
		glm::vec3 D65 = xyzCoefficients[6] * xyzToRGBMatrix;

		glm::vec3 rgbScattering = glm::vec3(xyzCoefficients[0] * xyzToRGBMatrix) / D65;
		glm::vec3 rgbExtinction = glm::vec3(xyzCoefficients[1] * xyzToRGBMatrix) / D65;
		glm::vec3 rgbAbsorption = glm::vec3(xyzCoefficients[2] * xyzToRGBMatrix) / D65;
		glm::vec3 rgbHostExtinction = glm::vec3(xyzCoefficients[3] * xyzToRGBMatrix) / D65;
		glm::vec3 rgbParticlesExtinction = glm::vec3(xyzCoefficients[4] * xyzToRGBMatrix) / D65;
		glm::vec3 rgbParticlesAbsorption = glm::vec3(xyzCoefficients[5] * xyzToRGBMatrix) / D65;

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::nanoseconds time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

		std::cout << "scattering = " << "vec3(" << rgbScattering.r << ", " << rgbScattering.g << ", " << rgbScattering.b << ")," << std::endl;
		std::cout << "absorption = " << "vec3(" << rgbAbsorption.r << ", " << rgbAbsorption.g << ", " << rgbAbsorption.b << ")," << std::endl;
		std::cout << "extinction = " << "vec3(" << rgbExtinction.r << ", " << rgbExtinction.g << ", " << rgbExtinction.b << ")," << std::endl;
		std::cout << " " << std::endl;
		std::cout << "hostExtinction = " << "vec3(" << rgbHostExtinction.r << ", " << rgbHostExtinction.g << ", " << rgbHostExtinction.b << ")," << std::endl;
		std::cout << "particleExtinction = " << "vec3(" << rgbParticlesExtinction.r << ", " << rgbParticlesExtinction.g << ", " << rgbParticlesExtinction.b << ")," << std::endl;
		std::cout << "particleAbsorption = " << "vec3(" << rgbParticlesAbsorption.r << ", " << rgbParticlesAbsorption.g << ", " << rgbParticlesAbsorption.b << ")," << std::endl;
		std::cout << std::endl;
		std::cout << time.count() << std::endl;
	#endif
}