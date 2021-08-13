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
		A[n] = tmp - 1.0 / (tmp + A[n + 1]);
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
		oldPsiZeta = oldPsiZeta * (n_z - A[n - 1])*(n_z * oldB);
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
		oldB = A[n] + complex<double>(0.0, 1.0) / PsiZeta(n, z, oldB, A);
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

	complex<double> B_n = B(n, hostZ, hostA);
	complex<double> R_n = R(n, hostZ, hostA);

	a = R_n * (iorHost * particleA[n] - iorParticle * hostA[n]) / (iorHost * particleA[n] - iorParticle * B_n);
	b = R_n * (iorParticle * particleA[n] - iorHost * hostA[n]) / (iorParticle * particleA[n] - iorHost * B_n);
}

unsigned int TermsToSum(const complex<double> z) {
	double size = abs(z);
	return static_cast<unsigned int>(ceil(size + 4.3 * cbrt(size) + 1.0));
}

double ComputeMiePhase(complex<double> iorHost, complex<double> iorParticle, double cosTheta, double radius, double lambda, complex<double>& S1, complex<double>& S2, double& Qabs, double& Qsca, double& Qext) {
	double size = 2.0 * pi * radius / lambda;
	M = TermsToSum(iorHost * size);
	AallN(hostA, iorHost * size);
	AallN(particleA, iorParticle * size);

	double sum = 0.0;

	LorenzMie_ab(1, size, iorHost, iorParticle);

	valarray<double> Pi;
	valarray<double> Tau;
	Pi.resize(M + 1);
	Tau.resize(M + 1);
	Pi[0] = 1.0;
	Pi[1] = 3.0 * cosTheta;
	Tau[0] = cosTheta;
	Tau[1] = 2.0 * cosTheta * Pi[1] - 3.0;
	for (int n = 2; n < M; ++n) {
		Pi[n] = ((2.0 * (n + 1.0) - 1.0) / n) * cosTheta * Pi[n - 1] - ((n + 1.0) / n) * Pi[n - 2];
		Tau[n] = (n + 1.0) * cosTheta * Pi[n] - (n + 2.0) * Pi[n - 1];
	}

	for (unsigned int n = 1; n < M; ++n) {
		complex<double> a_n = a;
		complex<double> b_n = b;

		double tmp = (2.0 * n + 1.0) / (n * (n + 1.0));
		S1 += tmp * (a_n * Pi[n] + b_n * Tau[n]);
		S1 += tmp * (a_n * Tau[n] + b_n * Pi[n]);
		sum += (2.0 * n + 1.0) * (sqr(abs(a_n)) + sqr(abs(b_n)));

		LorenzMie_ab(n + 1, size, iorHost, iorParticle);
	}
	double phase = sqr(abs(S1)) + sqr(abs(S2));
		   phase = phase / 4.0 * pi * sum;

	return phase;
}