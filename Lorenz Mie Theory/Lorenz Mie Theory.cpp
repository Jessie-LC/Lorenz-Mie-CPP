#include "Utility.h"
#include "MieTheory.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 1800;

int main() {
	complex<double> iorHost = complex<double>(1.335601, 2.46E-09);

	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.2e-6;
	double mineral_stepSize = 2e-7;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;
	double algae_stepSize = 5e-7;

	int counter = 0;
	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);
	for (double r = rMin_mineral + mineral_stepSize * 0.5; r < rMax_mineral; r += mineral_stepSize) {
		N_mineral[counter] = (3.07e-7 / 10.44) * pow(2.0 * r, -3.4);
		//N_mineral[counter] *= 3.0 / (4.0 * pi * r * r * r) * 1.0 / mineral_stepSize;
		std::cout << "Number Density Mineral: " << N_mineral[counter] << endl;
		++counter;
	}

	counter = 0;
	valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);
	for (double r = rMin_algae + algae_stepSize * 0.5; r < rMax_algae; r += algae_stepSize) {
		N_algae[counter] = (3.87e-7 / 4.97) * pow(2.0 * r, -3.6);
		//N_algae[counter] *= 3.0 / (4.0 * pi * r * r * r) * 1.0 / algae_stepSize;
		std::cout << "Number Density Algae: " << N_algae[counter] << endl;
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
	/*
	double rMax_water = 20.0e-6;
	double rMin_water = 1.0e-6; //If problems occur, make this smaller.
	double water_stepSize = 0.1e-6;

	LogNormalParticleDistribution CloudLogNormal{
		6.0e-6,
		CloudLogNormal.mean * 0.25
	};

	valarray<double> N_cloud;

	double x_vs = CloudLogNormal.mean;
	double beta_sqr = log(((CloudLogNormal.standardDeviation * CloudLogNormal.standardDeviation) / (x_vs * x_vs)) + 1.0);
	double alpha = log(x_vs) - 0.5 * beta_sqr;
	double a = 5.3e4;
	double beta =  sqrt(beta_sqr);

	int counter = 0;
	N_cloud.resize(static_cast<unsigned int>((rMax_water - rMin_water) / water_stepSize) + 1u);
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		double x = r;
		double tmp = (log(x) - alpha) / beta;
		double y = 2.1;
		N_cloud[counter] = 1.0 / (x * beta * sqrt(tau)) * exp(-0.5 * tmp * tmp);
		//N_cloud[counter] = a * pow(x, alpha) * exp(-beta * pow(x, y));
		std::cout << "Radius UM:" << r * 1e6 << "	Number Density:" << N_cloud[counter] << endl;
		++counter;
	}

	counter = 0;
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		N_cloud[counter] *= 3.0 / (4.0 * pi * r * r * r) * 1.0 / water_stepSize; //Not sure if this is required.
		++counter;
	}

	std::cout << "Number Distribution finished generating" << endl;

	ParticleDistribution clouds{
		complex<double>{ 1.335601, 2.46E-09 },
		rMin_water,
		rMax_water,
		water_stepSize,
		N_cloud,

		CloudLogNormal
	};
	*/

	ParticleDistribution Particles[2];
	Particles[0] = algae;
	Particles[1] = mineral;

	double lambda = 0.650e-6;
	
	valarray<double> mieData;
	mieData.resize(angles);

	double scatteringAlbedo = 0.0;
	double scatteringCoefficient = 0.0;
	double extinctionCoefficient = 0.0;
	double absorptionCoefficient = 0.0;
	std::cout << "Start generating volume properties via Mie Theory" << endl;
	std::cout << "Angle" << "  " << "Scattering" << "  " << "Extinction" << "  " << "Bulk Phase" << endl;
	for (int n = 0; n < angles; ++n) {
		double dtheta{ pi / angles };
		double theta = n * dtheta;
	
		double phase = 0.0;
		double scattering = 0.0;
		double extinction = 0.0;
		double absorption = 0.0;
		for (int i = 0; i < 2; ++i) {
			BulkMedium Bulk;
			ComputeBulkOpticalProperties(iorHost, theta, lambda, Particles[i], Bulk);

			double absorptionMedium = 4.0 * pi * imag(iorHost) / lambda;
			extinction += absorptionMedium + (Bulk.extinction);
			scattering += Bulk.scattering;
			absorption += Bulk.absorption;
			phase += Bulk.phase * Bulk.scattering;
		}
		phase = (1.0 / scattering) * phase;

		scatteringCoefficient = scattering;
		extinctionCoefficient = extinction;
		absorptionCoefficient = absorption;
		scatteringAlbedo = glm::clamp(scattering / extinction, 0.0, 10.0);

		mieData[n] = phase;

		std::cout << ((double)n / (double)angles) * 180.0 << "  " << scattering << "  " << extinction << "  " << phase << endl;
	}

	std::cout << "Finished generating volume properties!" << endl;

	glm::vec3 phaseTexture[angles];
	for (int x = 0; x < angles; ++x) {
		phaseTexture[x] = glm::vec3(mieData[x]);
	}

	glm::vec3 CDF[angles];
	for (int n = 0; n < angles; ++n) {
		float integral = 0.0;
		for (int m = 0; m < n; ++m) {
			float dTheta = pi / angles;
			integral += (mieData[m] * sin(m * dTheta)) * dTheta;
		}
		CDF[n] = glm::vec3(integral) / (float)(pi / 2.0);
	}

	glm::vec3 maxCDF = CDF[angles - 1];

	for (int n = 0; n < angles; ++n) {
		phaseTexture[n] /= maxCDF;
		CDF[n] /= maxCDF;
		std::cout << "CDF:	" << CDF[n].x << "	Phase:	" << phaseTexture[n].x << endl;
	}

	ofstream phaseLut("phase.dat", ios::binary);
	phaseLut.write(reinterpret_cast<char*>(phaseTexture), sizeof(glm::vec3) * angles * 1);
	phaseLut.close();
	std::cout << "Finished writing phase!" << endl;

	ofstream cdfLut("cdf.dat", ios::binary);
	cdfLut.write(reinterpret_cast<char*>(CDF), sizeof(glm::vec3) * angles * 1);
	cdfLut.close();
	std::cout << "Finished writing CDF lut!" << endl;

	std::cout << "Scattering Albedo:	" << scatteringAlbedo << endl;
	std::cout << "Scattering Coefficient:	" << scatteringCoefficient / maxCDF.x << endl;
	std::cout << "Extinction Coefficient:	" << extinctionCoefficient / maxCDF.x << endl;
	std::cout << "Absorption Coefficient:	" << absorptionCoefficient / maxCDF.x << endl;
}