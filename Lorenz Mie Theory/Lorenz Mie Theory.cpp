#include "Utility.h"
#include "MieTheory.h"
#include "IndexOfRefraction.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 1800;
const int wavelengths = 16;

void GeneratePhase(ParticleDistribution* particles, bool singleParticle) {
	glm::vec3 mieData[angles];

	for (int i = 0; i < angles; ++i) {
		mieData[i] = glm::vec3(0.0f);
	}

	glm::vec3 scatteringCoefficientRGB = glm::vec3(0.0f);
	glm::vec3 extinctionCoefficientRGB = glm::vec3(0.0f);

	double scatteringAlbedo = 0.0;
	double scatteringCoefficient = 0.0;
	double extinctionCoefficient = 0.0;
	double absorptionCoefficient = 0.0;
	valarray<double> anglesRadians;
	anglesRadians.resize(angles);
	for (int n = 0; n < angles; ++n) {
		double dtheta{ pi / (angles - 1) };
		double theta = n * dtheta;

		anglesRadians[n] = theta;

		for (int i = 0; i < wavelengths; ++i) {
			double r = (double)i / wavelengths;
			double lambda = (390.0 + (441.0 * r)) * 1e-9;

			complex<double> hostIOR = BrineIOR(lambda * 1e9);

			double absorptionMedium = 4.0 * pi * imag(hostIOR) / lambda;

			particles[0].ior = complex<double>(particles[0].ior.real(), InterpAlgaeK(lambda * 1e9));
			particles[1].ior = complex<double>(particles[1].ior.real(), InterpMineralK(lambda * 1e9));

			double phase = 0.0;
			if (!singleParticle) {
				double phaseAsymmetry = 0.0;
				double scattering = 0.0;
				double extinction = 0.0;
				double absorption = 0.0;
				for (int i = 0; i < 2; ++i) {
					BulkMedium Bulk;
					ComputeBulkOpticalProperties(hostIOR, theta, lambda, particles[i], Bulk);
					extinction += absorptionMedium + (Bulk.extinction);
					scattering += Bulk.scattering;
					phase += Bulk.phase * Bulk.scattering;
					phaseAsymmetry += Bulk.phaseAsymmetry * Bulk.scattering;
				}
				phase *= 1.0 / scattering;
				phaseAsymmetry *= 1.0 / scattering;

				scatteringCoefficient = scattering;
				extinctionCoefficient = extinction;
				absorptionCoefficient = extinction - scattering;
				scatteringAlbedo = abs(scatteringCoefficient / extinctionCoefficient);
			}
			else {
				complex<double> S1;
				complex<double> S2;
				double Qsca;
				double Qabs;
				double Qext;
				double phaseAsymmetry;
				ComputeParticleProperties(hostIOR, particles[0].ior, theta, particles[0].rMax, lambda, S1, S2, Qabs, Qsca, Qext, phase, phaseAsymmetry);

				scatteringCoefficient = Qsca;
				extinctionCoefficient = Qext;
				absorptionCoefficient = Qabs;
				scatteringAlbedo = abs(scatteringCoefficient / extinctionCoefficient);
			}

			mieData[n] += SpectrumToXYZ(phase, lambda * 1e9) / (float)wavelengths;
			scatteringCoefficientRGB += SpectrumToXYZ(scatteringCoefficient, lambda * 1e9) / (float)wavelengths;
			extinctionCoefficientRGB += SpectrumToXYZ(extinctionCoefficient, lambda * 1e9) / (float)wavelengths;
		}
		if (n < 1) {
			std::cout << "vec3(" << scatteringCoefficientRGB.r << ", " << scatteringCoefficientRGB.g << ", " << scatteringCoefficientRGB.b << ")," << std::endl;
			std::cout << "vec3(" << extinctionCoefficientRGB.r << ", " << extinctionCoefficientRGB.g << ", " << extinctionCoefficientRGB.b << ")," << std::endl;
		}

		std::cout << "vec3(" << mieData[n].r << ", " << mieData[n].g << ", " << mieData[n].b << ")," << std::endl;
	}

	glm::vec3 CDF[angles];
	CDF[0] = glm::vec3(0.0);
	for (int n = 1; n < angles; ++n) {
		float angleStart = anglesRadians[n - 1];
		float angleEnd = anglesRadians[n];
		glm::vec3 phaseStart = mieData[n - 1];
		glm::vec3 phaseEnd = mieData[n];

		glm::vec3 p1 = phaseStart * (cos(angleStart) - cos(angleEnd));
		glm::vec3 p2 = (phaseEnd - phaseStart) / (angleEnd - angleStart);
		float p3 = (sin(angleEnd) - angleEnd * cos(angleEnd)) - (sin(angleStart) - angleStart * cos(angleStart));
		float p4 = angleStart * (cos(angleStart) - cos(angleEnd));

		glm::vec3 integral = (float)2.0 * (float)pi * (p1 + p2 * (p3 - p4));
		CDF[n] = CDF[n - 1] + integral;
	}

	glm::vec3 maxCDF = CDF[angles - 1];

	for (int n = 0; n < angles; ++n) {
		mieData[n] /= maxCDF;
		CDF[n] /= maxCDF;
	}

	ofstream phaseLut("phase.dat", ios::binary);
	phaseLut.write(reinterpret_cast<char*>(mieData), sizeof(glm::vec3) * angles * 1);
	phaseLut.close();
	std::cout << "Finished writing phase!" << endl;
}

int main() {
	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.01e-6;
	double mineral_stepSize = 0.05e-6;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;
	double algae_stepSize = 0.2e-6;

	double vAlgae = 1.904e-6;
	double vMineral = 5.429e-7;

	int counter = 0;
	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);
	for (double r = rMin_mineral + mineral_stepSize * 0.5; r < rMax_mineral; r += mineral_stepSize) {
		N_mineral[counter] = (vMineral / 10.44) * pow(2.0 * r, -3.4);
		++counter;
	}

	counter = 0;
	valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);
	for (double r = rMin_algae + algae_stepSize * 0.5; r < rMax_algae; r += algae_stepSize) {
		N_algae[counter] = (vAlgae / 4.97) * pow(2.0 * r, -3.6);
		++counter;
	}

	ParticleDistribution mineral{
		complex<double>{ 1.58, 3.61e-4 },
		rMin_mineral,
		rMax_mineral,
		mineral_stepSize,
		N_mineral,

		LogNormalParticleDistribution{
			5.429e-7,
			0.25
		}
	};

	ParticleDistribution algae{
		complex<double>{ 1.41, 8.19e-5 },
		rMin_algae,
		rMax_algae,
		algae_stepSize,
		N_algae,

		LogNormalParticleDistribution{
			1.904e-6,
			0.15
		}
	};

	ParticleDistribution Particles[2];
	Particles[0] = algae;
	Particles[1] = mineral;

	GeneratePhase(Particles, false);
}