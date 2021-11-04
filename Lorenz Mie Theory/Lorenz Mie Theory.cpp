#include "Utility.h"
#include "MieTheory.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 2048;

int main() {
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const float mieData[] = float[](" << endl;

	std::cout << "Start generating Mie Phase:" << endl;

	double radius = 10.0e-6;

	complex<double> iorHost = complex<double>(1.00029, 0.0);

	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.02e-6;
	double mineral_stepSize = 1e-8;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;
	double algae_stepSize = 1e-8;

	int counter = 0;
	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);

	valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);

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

	double rMax_water = 50e-6;
	double rMin_water = 10e-6;
	double water_stepSize = 1e-6;

	LogNormalParticleDistribution CloudLogNormal{
		25.0e-6,
		0.15
	};

	valarray<double> N_cloud;

	double x_vs = CloudLogNormal.v;
	double beta_sqr = log(((CloudLogNormal.c * CloudLogNormal.c) / (x_vs * x_vs)) + 1.0);
	double alpha = log(x_vs) - 0.5 * beta_sqr;
	double beta = sqrt(beta_sqr);

	counter = 0;
	N_cloud.resize(static_cast<unsigned int>((rMax_water - rMin_water) / water_stepSize) + 1u);
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		double x = r;
		double tmp = (log(x) - alpha) / beta;
		N_cloud[counter] = 1.0 / (x * beta * sqrt(tau)) * exp(-0.5 * tmp * tmp);
		++counter;
	}

	counter = 0;
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		N_cloud[counter] *= pow(r, 3.0) / water_stepSize;
		++counter;
	}

	ParticleDistribution clouds{
		complex<double>{ 1.335601, 2.46E-09 },
		rMin_water,
		rMax_water,
		water_stepSize,
		N_cloud,

		CloudLogNormal
	};

	//double phaseTexture[angles * 1];

	double lambda = 0.550e-6;
	int counter1 = 0;
	for (int n = 0; n < angles; ++n) {
		counter1++;
		double dtheta{ pi / angles };
		double theta = n * dtheta;
		BulkMedium bulk_C;
		//BulkMedium bulk_M;
		ComputeBulkOpticalProperties(iorHost, theta, lambda, clouds, bulk_C);
		//ComputeBulkOpticalProperties(iorHost, theta, lambda, mineral, bulk_M);

		/*
		complex<double> S1;
		complex<double> S2;
		double Qsca;
		double Qabs;
		double Qext;
		double phase;
		ComputeParticleProperties(iorHost, clouds.ior, theta, 10e-6, lambda, S1, S2, Qabs, Qsca, Qext, phase);
		*/
	
		double absorptionMedium = 4.0 * pi * imag(iorHost) / lambda;
		double extinction = absorptionMedium + (bulk_C.extinction);
		double scattering = bulk_C.scattering;
		double absorption = bulk_C.absorption;
		double phase = (1.0 / scattering) * bulk_C.phase * bulk_C.scattering;

		mieOutput << "	" << phase << ", " << endl;
		std::cout << phase << endl;
	}

	cout << "Finished generating Mie Phase!" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}