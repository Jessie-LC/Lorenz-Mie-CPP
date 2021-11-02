#include "Utility.h"
#include "MieTheory.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 1800;

int main() {
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const vec3 mieData[] = vec3[](" << endl;

	std::cout << "Start generating Mie Phase:" << endl;

	double radius = 10.0e-6;

	complex<double> iorHost = complex<double>(1.33257, 1.67E-08);

	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.01e-6;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;

	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / 1e-8) + 1u);
	int counter = 0;
	for (double r = rMin_mineral + 1e-8 * 0.5; r < rMax_mineral; r += 1e-8) {
		N_mineral[counter] = 5.20019157e-8 * pow(2.0 * r, -3.4);
		++counter;
	}

	valarray<double> N_algae;
	counter = 0;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / 1e-7) + 1u);
	for (double r = rMin_algae + 1e-7 * 0.5; r < rMax_algae; r += 1e-7) {
		N_algae[counter] = 3.82828994e-7 * pow(2.0 * r, -3.6);
		++counter;
	}

	ParticleDistribution mineral{
		complex<double>{ 1.58, 3.61e-4 },
		rMin_mineral,
		rMax_mineral,
		1e-8,
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
		1e-7,
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

	double lambda = 0.550e-6;
	for (int n = 0; n < angles; ++n) {
		double dtheta{ pi / angles };
		double theta = n * dtheta;
		BulkMedium bulk_A;
		BulkMedium bulk_M;
		ComputeBulkOpticalProperties(iorHost, theta, lambda, algae, bulk_A);
		ComputeBulkOpticalProperties(iorHost, theta, lambda, mineral, bulk_M);

		double absorptionMedium = 4.0 * pi * imag(iorHost) / lambda;
		double extinction = absorptionMedium + (bulk_A.extinction + bulk_M.extinction);
		double scattering = bulk_A.scattering + bulk_M.scattering;
		double absorption = bulk_A.absorption + bulk_M.absorption;
		double phase = (1.0 / scattering) * ((bulk_A.phase * bulk_A.scattering) + (bulk_M.phase * bulk_M.scattering));

		mieOutput << phase << "," << endl;
		std::cout << phase << endl;
	}

	std::cout << "Finished generating Mie Phase" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}