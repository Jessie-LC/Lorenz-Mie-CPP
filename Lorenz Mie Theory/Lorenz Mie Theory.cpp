#include "Utility.h"
#include "MieTheory.h"

using namespace std;

int main() {
	const int iterations = 1800;
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const vec3 mieData[] = vec3[](" << endl;

	cout << "Start generating Mie Phase:" << endl;

	valarray<glm::dvec3> phase;
	phase.resize(iterations + 1);
	for (double wl = 390.0; wl <= 830.0; wl += 21.0) {
		int n = 0;
		int angle = 0;
		double lambda = double(wl * 1e-9);

		glm::dvec3 scatteringCoefficient;
		glm::dvec3 extinctionCoefficient;
		glm::dvec3 absorptionCoefficient;

		valarray<complex<double>> S1;
		valarray<complex<double>> S2;
		S1.resize(iterations + 1);
		S2.resize(iterations + 1);

		double sum;
		double Qsca;
		double Qabs;
		double Qext;

		for (int iteration = 0; iteration < iterations; ++iteration) {
			n++;
			angle++;

			double dtheta{ pi / iterations };
			double theta = iteration * dtheta;

			double radius = 10.0e-6;

			complex<double> iorHost = { 1.00029, 0.0 };

			complex<double> iorParticle = { 1.3330, 1.26e-4 };

			double phaseSpectral = ComputeMiePhase(iorHost, iorParticle, theta, radius, lambda, S1, S2, Qabs, Qsca, Qext, sum, angle);

			//double phase = sqr(abs(S1)) + sqr(abs(S2));
			//phase = phase / 4.0 * pi * sum;
			phase[angle] += ((sqr(abs(S1[angle])) + sqr(abs(S2[angle]))) / (4.0 * pi * sum)) * spectrumToXYZ(wl);
			scatteringCoefficient += Qsca * spectrumToXYZ(wl);
			extinctionCoefficient += Qext * spectrumToXYZ(wl);
			absorptionCoefficient += Qabs * spectrumToXYZ(wl);
		}
	}

	for (int n = 0; n < iterations; ++n) {
		phase[n] = xyzToRGB(phase[n]);
		cout << n << "    " << phase[n].r << ", " << phase[n].g << ", " << phase[n].b << endl;

		mieOutput << "    vec3(" << phase[n].r << ", " << phase[n].g << ", " << phase[n].b << ")," << endl;
	}
	cout << "Finished generating Mie Phase" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}