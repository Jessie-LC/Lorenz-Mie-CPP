#include "Utility.h"
#include "MieTheory.h"

#include <tuple>
#include <limits>

using namespace std;

const int angles = 3600;
const int wavelengths = 441;

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

double HenyeyGreensteinPhase(double cosTheta, double g) {
	const double norm = 0.25 / pi;

	double gg = g * g;
	return norm * ((1.0 - gg) / pow(1.0 + gg - 2.0 * g * cosTheta, 3.0 / 2.0));
}

double HenyeyGreensteinLegendre(double cosTheta, double g) {
	double gg = g * g;
	double P = 0.0;
	for (int i = 0; i < 500; ++i) {
		P += sgn((2.0 * i + 1.0) * pow(g, i)) * Pn(cosTheta, i);
	}
	return P / 12.5663706;
}

double CornetteShanks(double cosTheta, double g) {
	double gg = g * g;
	double p1 = 1.5 * ((1.0 - gg) / (2.0 + gg));
	double p2 = (1.0 + sqr(cosTheta)) / (pow((1.0 + gg - 2.0 * g * cosTheta), 3.0 / 2.0));
	double phase = (p1 * p2);
	phase /= 12.5663706;
	return phase;
}

double CornetteShanksLegendre(double cosTheta, double g) {
	double gg = g * g;
	double P = 0.0;
	for (int i = 0; i < 500; ++i) {
		double tmp1 = i * (i - 1) / (2.0 * i - 1.0);
		double tmp2 = (((5.0 * sqr((double)i) - 1.0) / (2.0 * i - 1.0)) + (sqr(i + 1.0) / (2.0 * i + 3.0))) * pow(g, i);
		double tmp3 = (((i + 1.0) * (i + 2.0)) / (2.0 * i + 3.0)) * pow(g, i + 2.0);
		P += sgn(tmp1 + tmp2 + tmp3) * Pn(cosTheta, i);
	}
	return 1.5 * (1.0 / (2.0 + gg)) * P;
}

/*

float LambdaSample() {
	float r = RandNextF(), g;
	int i;
	for (i = 0; i < 441 && r >= 0.0; i++) {
		g = dot(cie[i], 1.0 / vec3(113.042)) / 3.0;
		r -= g * 1.0;
	}

	return float(390 + i) + r / g;
}

float LambdaPdf(float w) {
	float n = (w - 390.0);
	int i = int(n);
	if(i < 0 || i >= (830-390)) {
		return 0.0;
	}
	return (dot(cie[i], 1.0 / vec3(113.042)) * 441.0);
}

*/

int main() {
	ofstream mieOutput;
	mieOutput.open("Mie Outputs.txt");

	complex<double> iorHost = complex<double>{ 1.00028, 0.0 };

	double rMax_mineral = 100e-6;
	double rMin_mineral = 0.05e-6;
	double mineral_stepSize = rMin_mineral;

	double rMax_algae = 100e-6;
	double rMin_algae = 0.225e-6;
	double algae_stepSize = 2e-7;

	mieOutput << "Mineral N" << endl;

	int counter = 0;
	valarray<double> N_mineral;
	N_mineral.resize(static_cast<unsigned int>((rMax_mineral - rMin_mineral) / mineral_stepSize) + 1u);
	for (double r = rMin_mineral + mineral_stepSize * 0.5; r < rMax_mineral; r += mineral_stepSize) {
		N_mineral[counter] = (5.429e-7 / 10.44) * pow(2.0 * r, -3.4);
		mieOutput << N_mineral[counter] << endl;
		++counter;
	}
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;

	mieOutput << "Algae N" << endl;
	counter = 0;
	valarray<double> N_algae;
	N_algae.resize(static_cast<unsigned int>((rMax_algae - rMin_algae) / algae_stepSize) + 1u);
	for (double r = rMin_algae + algae_stepSize * 0.5; r < rMax_algae; r += algae_stepSize) {
		N_algae[counter] = (1.904e-6 / 4.97) * pow(2.0 * r, -3.6);
		mieOutput << N_algae[counter] << endl;
		++counter;
	}
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;

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

	double rMax_water = 2e-5;
	double rMin_water = 2e-7;
	double water_stepSize = rMin_water * 2.0;

	LogNormalParticleDistribution CloudLogNormal{
		8.0e-7,
		CloudLogNormal.mean * 1.5
	};

	valarray<double> N_cloud;

	double beta_sqr = log(((CloudLogNormal.standardDeviation * CloudLogNormal.standardDeviation) / (CloudLogNormal.mean * CloudLogNormal.mean)) + 1.0);
	double alpha = log(CloudLogNormal.mean) - 0.5 * beta_sqr;
	double a = 4.9757e6;
	double beta = sqrt(beta_sqr);

	double volume = 0.0;

	counter = 0;
	N_cloud.resize(static_cast<unsigned int>((rMax_water - rMin_water) / water_stepSize) + 1u);
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		double x = r * 1e6;
		double tmp = (log(x) - alpha) / beta;
		double y = 0.5;
		double cumulus = 2.373 * pow(x, 6.0) * exp(-1.5 * x);
		double normal = (1.0 / (CloudLogNormal.standardDeviation * sqrt(tau))) * exp(-0.5 * pow((r - CloudLogNormal.mean) / CloudLogNormal.standardDeviation, 2.0));
		double logNormal = (1.0 / (r * beta * sqrt(tau))) * exp(-0.5 * pow((log(r) - alpha) / beta, 2.0));
		double haze = 5.33e4 * pow(x, 4.0) * exp(-8.9 * pow(x, 0.5));
		N_cloud[counter] = cumulus;
		volume += pow(r, 3.0) * N_cloud[counter] * water_stepSize;
		++counter;
	}

	double waterWeight = 0.0005;
	double waterDensity = 1000.0;
	double waterVolume = waterWeight < 1e-12 ? 0.0 : waterWeight / waterDensity;
	double airVolume = 1.0 - waterVolume;

	double N = ((waterWeight > 1e-12 ? waterVolume : 1.0) / airVolume) / (4.1887902 * volume);

	counter = 0;
	for (double r = rMin_water + water_stepSize * 0.5; r < rMax_water; r += water_stepSize) {
		N_cloud[counter] = N * N_cloud[counter];
		std::cout << "Radius: " << r << " DSD: " << N_cloud[counter] << endl;
		++counter;
	}

	ParticleDistribution clouds{
		complex<double>{ 1.3330, 1.9600e-9 },
		rMin_water,
		rMax_water,
		water_stepSize,
		N_cloud,

		CloudLogNormal
	};

	ParticleDistribution Particles[2];
	Particles[0] = clouds;
	Particles[1] = mineral;

	double lambda = 0.550e-6;
	
	glm::vec3 mieData[angles];

	double absorptionMedium = 4.0 * pi * imag(iorHost) / lambda;

	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;

	double scatteringAlbedo = 0.0;
	double scatteringCoefficient = 0.0;
	double extinctionCoefficient = 0.0;
	double absorptionCoefficient = 0.0;
	valarray<double> anglesRadians;
	anglesRadians.resize(angles);

	auto t0 = std::chrono::high_resolution_clock::now();
	for (int n = 0; n < angles; ++n) {
		double dtheta{ pi / (angles - 1) };
		double theta = n * dtheta;

		anglesRadians[n] = theta;

		//*
		double phase = 0.0;
		double phaseAsymmetry = 0.0;
		double scattering = 0.0;
		double extinction = 0.0;
		double absorption = 0.0;
		for (int i = 0; i < 1; ++i) {
			BulkMedium Bulk;
			ComputeBulkOpticalProperties(iorHost, theta, lambda, Particles[i], Bulk);
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
		scatteringAlbedo = glm::clamp(abs(scatteringCoefficient / extinctionCoefficient), 0.0, 10.0);
		//*/

		/*
		complex<double> S1;
		complex<double> S2;
		double Qsca;
		double Qabs;
		double Qext;
		double phase;
		double phaseAsymmetry;
		ComputeParticleProperties(iorHost, complex<double>{ 1.3330, 1.9600e-9 }, theta, rMax_water, lambda, S1, S2, Qabs, Qsca, Qext, phase, phaseAsymmetry);

		scatteringCoefficient = Qsca;
		extinctionCoefficient = Qext;
		absorptionCoefficient = Qabs;
		scatteringAlbedo = glm::clamp(abs(scatteringCoefficient / extinctionCoefficient), 0.0, 10.0);
		//*/

		mieData[n] = glm::vec3(phase);

		cout << phase << endl;

		//fprintf(stderr, "\b\b\b\b\%3d%c", (100 * n / angles), '%');
	}
	cout << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;

	//*
	glm::vec3 phaseTexture[angles];
	for (int x = 0; x < angles; ++x) {
		phaseTexture[x] = mieData[x];
	}

	glm::vec3 CDF[angles];
	CDF[0] = glm::vec3(0.0);
	for (int n = 1; n < angles; ++n) {
		float angleStart = anglesRadians[n - 1];
		float angleEnd   = anglesRadians[n];
		glm::vec3 phaseStart = mieData[n - 1];
		glm::vec3 phaseEnd   = mieData[n];

		glm::vec3 p1 = phaseStart * (cos(angleStart) - cos(angleEnd));
		glm::vec3 p2 = (phaseEnd - phaseStart) / (angleEnd - angleStart);
		float p3 = (sin(angleEnd) - angleEnd * cos(angleEnd)) - (sin(angleStart) - angleStart * cos(angleStart));
		float p4 = angleStart * (cos(angleStart) - cos(angleEnd));

		glm::vec3 integral = (float)2.0 * (float)pi * (p1 + p2 * (p3 - p4));
		CDF[n] = CDF[n - 1] + integral;
	}

	glm::vec3 maxCDF = CDF[angles - 1];

	mieOutput << "Phase" << endl;
	for (int n = 0; n < angles; ++n) {
		phaseTexture[n] /= maxCDF;
		CDF[n] /= maxCDF;
		mieOutput << phaseTexture[n].x << endl;
		std::cout << "CDF:	" << CDF[n].x << "	Phase:	" << phaseTexture[n].x << endl;
	}

	std::cout << maxCDF.x << endl;

	mieOutput << " " << endl;
	mieOutput << " " << endl;
	mieOutput << " " << endl;

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

	mieOutput << "Coefficients" << endl;
	mieOutput << scatteringCoefficient / maxCDF.x << endl;
	mieOutput << extinctionCoefficient / maxCDF.x << endl;
	mieOutput << absorptionCoefficient / maxCDF.x << endl;
	//*/
}