#include "Utility.h"
#include "MieTheory.h"

using namespace std;

int main() {
	const int iterations = 1800;
	ofstream mieOutput;
	mieOutput.open("Generate Mie Phase.txt");
	mieOutput << "const vec3 mieData[] = vec3[](" << endl;

	cout << "Start generating Mie Phase:" << endl;

	int n = 0;
	int angle = 0;
	for (int iteration = 0; iteration < iterations; ++iteration) {
		n++;
		angle++;

		double dtheta { pi / iterations };
		double theta = iteration * dtheta;

		double radius = 1.0e-6;

		complex<double> iorHost = { 1.00029, 0.0 };

		complex<double> iorParticle = { 1.3330, 1.26e-9 };

		glm::dvec3 phase = glm::dvec3(0.0);
		glm::dvec3 scatteringCoefficient;
		glm::dvec3 extinctionCoefficient;
		glm::dvec3 absorptionCoefficient;
		
		complex<double> S1;
		complex<double> S2;
		double Qsca;
		double Qabs;
		double Qext;
		for (double wl = 390.0; wl <= 830.0; wl += 5.0) {
			double lambda = double(wl * 1e-9);
			double phaseSpectral = ComputeMiePhase(iorHost, iorParticle, theta, radius, lambda, S1, S2, Qabs, Qsca, Qext);
			phase += phaseSpectral * spectrumToXYZ(wl);
			scatteringCoefficient += Qsca * spectrumToXYZ(wl);
			extinctionCoefficient += Qext * spectrumToXYZ(wl);
			absorptionCoefficient += Qabs * spectrumToXYZ(wl);
		}
		phase = xyzToRGB(phase);
		scatteringCoefficient = xyzToRGB(scatteringCoefficient);
		extinctionCoefficient = xyzToRGB(extinctionCoefficient);
		absorptionCoefficient = xyzToRGB(absorptionCoefficient);

		cout << angle << "	" << phase.r << ", " << phase.g << ", " << phase.b << endl;

		mieOutput << "	vec3(" << phase.r << ", " << phase.g << ", " << phase.b << ")," << endl;
	}
	cout << "Finished generating Mie Phase" << endl;
	mieOutput << ");" << endl;
	mieOutput.close();
}