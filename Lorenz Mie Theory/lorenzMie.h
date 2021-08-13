#include "constants.h"
#include "Utility.h"
#include "MieTheory.h"

double computeMieTheory(complex<double> iorMedium, complex<double> iorParticle, double cosTheta, double r, double lambda) {
	double k0 = 2.0 * pi / lambda;
	complex<double> k = iorMedium * k0;
	complex<double> x = k * r;
	complex<double> m_r = iorParticle / iorMedium;
	complex<double> y = x * m_r;

	const double p = 4.3;
	int N = floor(abs(x) + 4 * cbrt(abs(x)) + 2);
	int M = max(N, static_cast<int>(ceil(abs(y)))) + 15;

	//---------------Logarithmic derivative DD calculated by downward recurrence-------------------
	// beginning with initial value (0.,0.) at M
	// y:=x*m_r   argument of DD
	vector<complex<double>> DD(M, { 0.0 + 0.0i });
	for (int n = M - 1; n > 1; --n)
		DD[n - 2] = (1.0 * n) / y - 1.0 / (DD[n - 1] + (1.0 * n) / y);

	//-----------Reccati-Bessel function PSI(0:N) calculated by downward recurrence----------
	//Reference: Mischenko et al. 2002, Scattering, Absorption and Emission of Light by Small Particles 3rd, pp.167-169
	// x:=k*r_p argument of PSI
	//RR=0.0d0; ! R(n):=PSI(n)/PSI(n-1)
	vector<complex<double>> RR(M + 1, { 0.0 + 0.0i });
	RR[M] = x / (2.0 * M + 1.0); //starting value of downward recurrence
	for (int n = M - 1; n > -1; --n)
		RR[n] = 1.0 / ((2.0 * n + 1.0) / x - RR[n + 1]); // R(n) := Rn

	vector<complex<double>> psi(M + 1, { 0.0 + 0.0i });
	psi[0] = RR[0] * cos(x);
	for (int n = 1; n <= N; ++n)
		psi[n] = RR[n] * psi[n - 1]; // PSI(n) := PSIn

	// -------Reccati-Bessel function chi(0:N) calculated by upward recurrence------
	// beginning with initial values chi(-1)=sin(x), chi(0)=-cos(x)
	// Reference: Bphren and Huffman 1983 (Appendix A, p.478)
	// CHI(x):= x*y(x) where y(x) is spherical bessel function of second kind
	// This is contrast to the definition of BH83: CHI(x):= -x*y(x)
	// x:=k*r_p argument of CHI
	vector<complex<double>> chi(N + 1, { 0.0 + 0.0i });
	chi[0] = -cos(x);
	chi[1] = (1.0 / x) * chi[0] - sin(x);
	for (int n = 2; n <= N; ++n)
		chi[n] = ((2.0 * n - 1.0) / x) * chi[n - 1] - chi[n - 2];

	vector<complex<double>> qsi(N + 1, { 0.0 + 0.0i });
	for (int n = 0; n <= N; ++n)
		qsi[n] = psi[n] + chi[n] * complex<double>(0.0 + 1.0i);

	//-----------Reccati-Bessel function PSIM(0:N) calculated by downward recurrence----------
	// Reference: Mischenko et al. 2002, Scattering, Absorption and Emission of Light by Small Particles 3rd, pp.167-169
	// y:=x*m_r argument of PSI
	fill(RR.begin(), RR.end(), 0.0 + 0.0i);
	RR[M] = y / (2.0 * M + 1.0);
	for (int n = M - 1; n > -1; --n)
		RR[n] = 1.0 / ((2.0 * n + 1.0) / y - RR[n + 1]); // R(n) := Rn

	vector<complex<double>> psim(N + 1, { 0.0 + 0.0i });
	psim[0] = RR[0] * cos(y);
	for (int n = 1; n <= N; ++n)
		psim[n] = RR[n] * psim[n - 1];

	vector<complex<double>> a(N, { 0.0 + 0.0i });
	vector<complex<double>> b(N, { 0.0 + 0.0i });
	for (int n = 0; n < N; ++n) {
		a[n] = ((iorParticle / (iorMedium * m_r) * DD[n] + (n + 1.0) / x) * psi[n + 1] - psi[n]) / ((iorParticle / (iorMedium * m_r) * DD[n] + (n + 1.0) / x) * qsi[n + 1] - qsi[n]); // BH83, Eq.4.88
		b[n] = ((iorMedium * m_r / iorParticle * DD[n] + (n + 1.0) / x) * psi[n + 1] - psi[n]) / ((iorMedium * m_r / iorParticle * DD[n] + (n + 1.0) / x) * qsi[n + 1] - qsi[n]); // BH83, Eq.4.88
	}

	vector<double> fn1(N, { 0.0 }), fn2(N, { 0.0 });
	for (int n = 1; n < N; ++n) {
		fn1[n] = 2.0 * (n + 1) + 1.0;
		fn2[n] = fn1[n] / ((n + 1.0) * (n + 2.0));
	}

	vector<double> Pi(N, { 0.0 });
	vector<double> Tau(N, { 0.0 });
	//Pi[0] = 1.0;
	//Pi[1] = 3 * cosTheta * Pi[0];
	//Tau[0] = cosTheta * Pi[0];
	//Tau[1] = 2 * cosTheta * Pi[1] - 3.0 * Pi[0];
	for (int n = 0; n < N; ++n) {
		Pi[n] = computePi(cosTheta, n);// ((2.0 * (n + 1.0) - 1.0) / n)* cosTheta* (Pi[n - 1] - (((n + 1.0) / n), Pi[n - 2]));
		Tau[n] = computeTau(cosTheta, n);// (n + 1.0)* cosTheta* (Pi[n] - ((n + 2.0) * Pi[n - 1]));
	}

	complex<double> S1 = complex<double>{ 0.0, 0.0 };
	complex<double> S2 = complex<double>{ 0.0, 0.0 };
	double sum = 0.0;
	for (int n = 0; n < N; ++n) {
		S1 += fn2[n] * (a[n] * Pi[n] + b[n] + Tau[n]);
		S2 += fn2[n] * (b[n] * Pi[n] + a[n] + Tau[n]);

		sum += (2.0 * n + 1.0) * (pow(abs(a[n]), 2.0) + pow(abs(b[n]), 2.0));
	}
	double phase = (pow(abs(S1), 2.0) + pow(abs(S2), 2.0)) / ((4.0 * pi) * sum);

	return phase;
}