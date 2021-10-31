#include "Utility.h"
#include "MieTheory.h"

using namespace std;

int main() {
	const int iterations = 1800;
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const vec3 mieData[] = vec3[](" << endl;

	std::cout << "Start generating Mie Phase:" << endl;

	int angle = 0;
	for (int n = 0; n < iterations; ++n) {
		angle++;

		double dtheta { pi / iterations };
		double theta = n * dtheta;

		double radius = 10.0e-6;

		complex<double> iorHost = { 1.34, 2.42e-9 };

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

		glm::dvec3 phaseAlgae = glm::dvec3(0.0);
		glm::dvec3 scatteringCoefficientAlgae;
		glm::dvec3 extinctionCoefficientAlgae;
		glm::dvec3 absorptionCoefficientAlgae;

		glm::dvec3 phaseMineral = glm::dvec3(0.0);
		glm::dvec3 scatteringCoefficientMineral;
		glm::dvec3 extinctionCoefficientMineral;
		glm::dvec3 absorptionCoefficientMineral;
		double wavelengthIncrement = 20.0;
		for (double wl = 390.0; wl <= 831.0; wl += wavelengthIncrement) {
			//ComputeMediumPhase(complex<double> iorHost, double theta, double lambda, ParticleDistribution& particle)
			double lambda = double(wl * 1e-9);
			double phaseMineralSpectral = ComputeMediumPhase(iorHost, theta, lambda, mineral);
			double phaseAlgaeSpectral   = ComputeMediumPhase(iorHost, theta, lambda, algae);
			phaseMineral += spectrumToXYZ(glm::dvec3(phaseMineralSpectral), wl);
			phaseAlgae   += spectrumToXYZ(glm::dvec3(phaseAlgaeSpectral), wl);
			scatteringCoefficientAlgae   += algae.scattering;
			scatteringCoefficientMineral += mineral.scattering;
		}
		scatteringCoefficientAlgae /= 441.0 / wavelengthIncrement;
		scatteringCoefficientMineral /= 441.0 / scatteringCoefficientMineral;
		phaseMineral /= 441.0 / wavelengthIncrement;
		phaseAlgae /= 441.0 / wavelengthIncrement;
		phaseMineral = xyzToRGB(phaseMineral);
		phaseAlgae   = xyzToRGB(phaseAlgae);
		scatteringCoefficientAlgae   = xyzToRGB(scatteringCoefficientAlgae);
		scatteringCoefficientMineral = xyzToRGB(scatteringCoefficientMineral);

		glm::dvec3 phase = phaseMineral;

		std::cout << cos(theta) << "	vec3(" << phase.r << ", " << phase.g << ", " << phase.b << ")," << endl;

		mieOutput << "	vec3(" << phase.r << ", " << phase.g << ", " << phase.b << ")," << endl;
	}
	std::cout << "Finished generating Mie Phase" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}