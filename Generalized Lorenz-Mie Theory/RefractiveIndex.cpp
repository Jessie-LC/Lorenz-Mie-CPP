#include "RefractiveIndex.h"
#include <glm.hpp>
#include "Constants.h"

int BinarySearch(int lowIndex, int highIndex, double toFind, const double arr[]) {
	while (lowIndex < highIndex) {
		int midIndex = (lowIndex + highIndex) >> 1;
		if (arr[midIndex] < toFind) {
			lowIndex = midIndex + 1;
		}
		else if (arr[midIndex] > toFind) {
			highIndex = midIndex;
		}
		else {
			return midIndex;
		}
	}
	return -1;
}

std::complex<double> IceIOR(double wavelength) {
	return std::complex<double>(sqrt(1 + 0.496 / (1 - pow(0.071 / wavelength, 2)) + 0.190 / (1 - pow(0.134 / wavelength, 2))), 0.0);
}

std::complex<double> BrineIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, OceanWavelengths[0]);
	double end = std::min(lambdaMax, OceanWavelengths[16]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 16, start, OceanWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double brineK = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 17 && end >= OceanWavelengths[entry]; ++entry) {
			double a = OceanWavelengths[entry],
				b = OceanWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = BrineK[entry],
				fb = BrineK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		brineK = result / (lambdaMax - lambdaMin);
	}

	brineK = glm::clamp(brineK, 1e-15, brineK);

	const double n1 = 1.779e-4;
	const double n2 = -1.05e-6;
	const double n3 = 1.6e-8;
	const double n4 = -2.02e-6;
	const double n5 = 15.868;
	const double n6 = 0.01155;
	const double n7 = -0.00423;
	const double n8 = -4382.0;
	const double n9 = 1.1455e6;
	const double T = 20.0;
	const double S = 0.0;
	double brineN = 1.31405 + (n1 + n2 * T + n3 * pow(T, 2.0)) * S + n4 * pow(T, 2.0) + ((n5 + n6 * S + n7 * T) / wavelength) + (n8 / pow(wavelength, 2.0)) + (n9 / pow(wavelength, 3.0));

	return std::complex<double>(brineN, brineK);
}

double InterpAlgaeK(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, OceanWavelengths[0]);
	double end = std::min(lambdaMax, OceanWavelengths[16]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 16, start, OceanWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double algaeK = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 17 && end >= OceanWavelengths[entry]; ++entry) {
			double a = OceanWavelengths[entry],
				b = OceanWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = AlgaeK[entry],
				fb = AlgaeK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		algaeK = result / (lambdaMax - lambdaMin);
	}

	return algaeK;
}

double InterpMineralK(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, OceanWavelengths[0]);
	double end = std::min(lambdaMax, OceanWavelengths[16]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 16, start, OceanWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double mineralK = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 17 && end >= OceanWavelengths[entry]; ++entry) {
			double a = OceanWavelengths[entry],
				b = OceanWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = MineralK[entry],
				fb = MineralK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		mineralK = result / (lambdaMax - lambdaMin);
	}

	return mineralK;
}

std::complex<double> WaterIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, WaterWavelengths[0]);
	double end = std::min(lambdaMax, WaterWavelengths[31]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 31, start, WaterWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double waterN = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 31 && end >= WaterWavelengths[entry]; ++entry) {
			double a = WaterWavelengths[entry],
				b = WaterWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = WaterN[entry],
				fb = WaterN[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		waterN = result / (lambdaMax - lambdaMin);
	}

	double waterK = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 31 && end >= WaterWavelengths[entry]; ++entry) {
			double a = WaterWavelengths[entry],
				b = WaterWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = WaterK[entry],
				fb = WaterK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		waterK = result / (lambdaMax - lambdaMin);
	}

	return std::complex<double>(waterN, waterK);
}

std::complex<double> IronIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, IronWavelengths[0]);
	double end = std::min(lambdaMax, IronWavelengths[43]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 43, start, IronWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double N = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 43 && end >= IronWavelengths[entry]; ++entry) {
			double a = IronWavelengths[entry],
				b = IronWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = IronN[entry],
				fb = IronN[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		N = result / (lambdaMax - lambdaMin);
	}

	double K = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 43 && end >= IronWavelengths[entry]; ++entry) {
			double a = IronWavelengths[entry],
				b = IronWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = IronK[entry],
				fb = IronK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		K = result / (lambdaMax - lambdaMin);
	}
	return std::complex<double>(N, K);
}

std::complex<double> CopperIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = std::max(lambdaMin, CopperWavelengths[0]);
	double end = std::min(lambdaMax, CopperWavelengths[76]);

	int idx = int(std::max(glm::distance(0.0, double(BinarySearch(0, 76, start, CopperWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double N = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 76 && end >= CopperWavelengths[entry]; ++entry) {
			double a = CopperWavelengths[entry],
				b = CopperWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = CopperN[entry],
				fb = CopperN[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		N = result / (lambdaMax - lambdaMin);
	}

	double K = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 76 && end >= CopperWavelengths[entry]; ++entry) {
			double a = CopperWavelengths[entry],
				b = CopperWavelengths[entry + 1],
				ca = std::max(a, start),
				cb = std::min(b, end),
				fa = CopperK[entry],
				fb = CopperK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		K = result / (lambdaMax - lambdaMin);
	}
	return std::complex<double>(N, K);
}

std::complex<double> Silica(double wavelength) {
    int index = int(((wavelength - 390.0) / 441.0) * 77.0);
    
    return std::complex<double>(SilicaN[index], SilicaK[index]);
}

std::complex<double> AirIOR(double wavelength) {
	double real = 1.0 + 8.06051E-5 + 2.480990E-2 / (132.274 - pow(wavelength, -2.0)) + 1.74557E-4 / (39.32957 - pow(wavelength, -2.0));
	double imaginary = 0.0;
	return std::complex<double>(real, imaginary);
}