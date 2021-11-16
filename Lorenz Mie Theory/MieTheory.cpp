#include "Utility.h"
#include "MieTheory.h"

unsigned int M = 0u;
valarray<complex<double>> particleA;
valarray<complex<double>> hostA;
complex<double> a, b;

void AallN(valarray<complex<double>>& A, const complex<double>& z, unsigned int userM = 0) {
	if (userM > M) M = userM;

	A.resize(M + 2);

	A[M + 1] = complex<double>(0.0, 0.0);

	if (M == 0) return;

	for (unsigned int n = M; n <= M; --n) {
		complex<double> tmp = (n + 1.0) / z;
		A[n] = tmp - pow(tmp + A[n + 1], -1.0);
	}
}

complex<double> PsiZeta(unsigned int n, const complex<double>& z, const complex<double>& oldB, const valarray<complex<double>>& A) {
	static unsigned int old_n = 0;
	static complex<double> oldPsiZeta;

	if (old_n == 0) {
		oldPsiZeta = 0.5 * (1.0 - exp(complex<double>(-2.0 * z.imag(), 2.0 * z.real())));
	}
	if (n == old_n + 1) {
		complex<double> n_z = static_cast<double>(n) / z;
		oldPsiZeta = oldPsiZeta * (n_z - A[n - 1])*(n_z - oldB);
		if (n == M) {
			old_n = 0;
		}
		else {
			old_n = n;
		}
	} else if (n != old_n && n != M) {
		cerr << "Error in computation of PSIn*ZETAn" << endl;
	}

	return oldPsiZeta;
}

complex<double> B(unsigned int n, const complex<double>& z, valarray<complex<double>>& A) {
	static unsigned int old_n = 0;
	static complex<double> oldB;

	if (old_n == 0) {
		oldB = complex<double>(0.0, 1.0);
	}

	if (n > A.size()) {
		AallN(A, z, n);
	}

	if (n == old_n + 1) {
		oldB = A[n] + (complex<double>(0.0, 1.0) / PsiZeta(n, z, oldB, A));
		if (n == M) {
			old_n = 0;
		}
		else {
			old_n = n;
		}
	}
	else if (n != old_n && n != M) {
		cerr << "Error in computation of Bn" << endl;
	}

	return oldB;
}

complex<double> R(unsigned int n, const complex<double>& z, valarray<complex<double>>& A) {
	static unsigned int old_n = 0;
	static complex<double> oldR;

	if (old_n == 0) {
		oldR = 0.5 * (1.0 - exp(complex<double>(2.0 * z.imag(), -2.0 * z.real())));
	}

	if (n > A.size()) {
		AallN(A, z, n);
	}

	if (n == old_n + 1) {
		complex<double> n_z = static_cast<double>(n) / z;
		oldR = oldR * (B(n, z, A) + n_z) / (A[n] + n_z);
		if (n == M) {
			old_n = 0;
		}
		else {
			old_n = n;
		}
	}
	else if (n != old_n && n != M) {
		cerr << "Error in computation of Rn" << endl;
	}

	return oldR;
}

complex<double> SphJn_Derivative(int n, complex<double> z) {
	return SphJn(n - 1, z) - ((double)(n + 1) / z) * SphJn(n, z);
}

complex<double> SphYn_Derivative(int n, complex<double> z) {
	return SphYn(n - 1, z) - ((double)(n + 1) / z) * SphYn(n, z);
}

complex<double> Psi(int n, complex<double> z) {
	return z * SphJn(n, z);
}

complex<double> PsiPrime(int n, complex<double> z) {
	return 1.0 * SphJn(n, z) + z * SphJn_Derivative(n, z);
}

complex<double> Zeta(int n, complex<double> z) {
	return z * (SphJn(n, z) - complex<double>(0.0, 1.0) * SphYn(n, z));
}

complex<double> ZetaPrime(int n, complex<double> z) {
	complex<double> g = SphJn(n, z) - complex<double>(0.0, 1.0) * SphYn(n, z);
	complex<double> gPrime = SphJn_Derivative(n, z) - (complex<double>(0.0, 1.0) * SphYn_Derivative(n, z));
	return 1.0 * g + z * gPrime;
}

void LorenzMie_ab(unsigned int n, double size, const complex<double>& iorHost, const complex<double>& iorParticle) {
	static unsigned int old_n = 0;

	if (n == old_n) {
		return;
	}
	if (n != old_n + 1) {
		cerr << "Error: The Lorenz-Mie coefficients must be computed consecutively starting from n=1 and counting upwards to n=M" << endl;
		exit(0);
	}
	if (n == M) {
		old_n = 0;
	}
	else {
		old_n = n;
	}

	complex<double> hostZ = iorHost * size;
	complex<double> particleZ = iorParticle * size;

	/*
	a = (iorHost * PsiPrime(n, particleZ) * Psi(n, hostZ) - iorParticle * Psi(n, particleZ) * PsiPrime(n, hostZ)) / 
		(iorHost * PsiPrime(n, particleZ) * Zeta(n, hostZ) - iorParticle * Psi(n, particleZ) * ZetaPrime(n, hostZ));
	b = (iorParticle * PsiPrime(n, particleZ) * Psi(n, hostZ) - iorHost * Psi(n, particleZ) * PsiPrime(n, hostZ)) /
		(iorParticle * PsiPrime(n, particleZ) * Zeta(n, hostZ) - iorHost * Psi(n, particleZ) * ZetaPrime(n, hostZ));
	*/

	complex<double> B_n = B(n, hostZ, hostA);
	complex<double> R_n = R(n, hostZ, hostA);

	a = R_n * (iorHost * particleA[n] - iorParticle * hostA[n]) / (iorHost * particleA[n] - iorParticle * B_n);
	b = R_n * (iorParticle * particleA[n] - iorHost * hostA[n]) / (iorParticle * particleA[n] - iorHost * B_n);
}

unsigned int TermsToSum(const complex<double> z) {
	double size = abs(z);
	return static_cast<unsigned int>(ceil(size + 4.3 * cbrt(size) + 1.0));
}

void ComputeParticleProperties(complex<double> iorHost, complex<double> iorParticle, double theta, double radius, double lambda, complex<double>& S1, complex<double>& S2, double& Qabs, double& Qsca, double& Qext, double& phase) {
	/*
		https://cseweb.ucsd.edu//~henrik/papers/lorenz_mie_theory/computing_scattering_properties_using_lorenz_mie_theory.pdf

		Secontion 2: Scattering in Participating Media
	*/
	double size = 2.0 * pi * radius / lambda;
	M = TermsToSum(iorHost * size);
	AallN(hostA, iorHost * size);
	AallN(particleA, iorParticle * size);

	double sum = 0.0;

	LorenzMie_ab(1, size, iorHost, iorParticle);

	double cosTheta = cos(theta);

	double crossSectionEXT = 0.0;
	for (unsigned int n = 1; n < M; ++n) {
		double PiN = computePi(cosTheta, n);
		double TauN = computeTau(cosTheta, n);

		complex<double> a_n = a;
		complex<double> b_n = b;

		double tmp = (2.0 * n + 1.0) / (n * (n + 1.0));
		S1 += tmp * (a_n * PiN + b_n * TauN);
		S2 += tmp * (a_n * TauN + b_n * PiN);
		sum += (2.0 * n + 1.0) * (sqr(abs(a_n)) + sqr(abs(b_n)));

		crossSectionEXT += (2.0 * n + 1.0) * real((a_n + b_n) / (iorHost * iorHost));

		LorenzMie_ab(n + 1, size, iorHost, iorParticle);
	}

	double alpha = 4.0 * pi * radius * imag(iorHost) / lambda;
	double y = alpha < 10e-6 ? 1.0 : (2.0 * (1.0 + (alpha - 1.0) * exp(alpha))) / pow(alpha, 2.0);
	double term1 = pow(lambda, 2.0) * exp(-alpha);
	double term2 = 2.0 * pi * y * sqr(abs(iorHost));

	Qsca = term1 / term2 * sum;
	Qext = (sqr(lambda) / tau) * crossSectionEXT;

	Qabs = Qext - Qsca;

	complex<double> k = 2.0 * pi * iorHost / lambda;

	phase = (sqr(abs(S1)) + sqr(abs(S2))) / (2.0 * sqr(abs(k)) * Qsca);
}

void ComputeBulkOpticalProperties(complex<double> iorHost, double theta, double lambda, ParticleDistribution& particle, BulkMedium& bulk) {
	/*
		https://cseweb.ucsd.edu//~henrik/papers/lorenz_mie_theory/computing_scattering_properties_using_lorenz_mie_theory.pdf

		Secontion 2.2: Bulk Optical Properties
	*/
	double phase = 0.0;
	double scatteringSigma = 0.0;
	double extinctionSigma = 0.0;
	int counter = 0;
	for (double r = particle.rMin + particle.stepSize * 0.5; r < particle.rMax; r += particle.stepSize) {
		complex<double> S1;
		complex<double> S2;
		double Qsca;
		double Qabs;
		double Qext;
		double particlePhase;
		ComputeParticleProperties(iorHost, particle.ior, theta, r, lambda, S1, S2, Qabs, Qsca, Qext, particlePhase);

		double sigmaS = Qsca * particle.N[counter] * particle.stepSize;

		scatteringSigma += sigmaS;
		extinctionSigma += Qext * particle.N[counter] * particle.stepSize;
		phase += Qsca * particlePhase * particle.N[counter] * particle.stepSize;

		++counter;
	}

	bulk.phase = (1.0 / scatteringSigma) * phase;
	bulk.scattering = scatteringSigma;
	bulk.extinction = extinctionSigma;
	bulk.absorption = extinctionSigma - scatteringSigma;
}