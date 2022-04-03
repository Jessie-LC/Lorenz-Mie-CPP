#include "IndexOfRefraction.h"

#include "constants.h"

using namespace std;

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

complex<double> WaterIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = glm::max(lambdaMin, WaterWavelengths[0]);
	double end = glm::min(lambdaMax, WaterWavelengths[31]);

	int idx = int(glm::max(glm::distance(0.0, double(BinarySearch(0, 31, start, WaterWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double waterN = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 32 && end >= WaterWavelengths[entry]; ++entry) {
			double a = WaterWavelengths[entry],
				b = WaterWavelengths[entry + 1],
				ca = max(a, start),
				cb = min(b, end),
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
		for (; entry < 32 && end >= WaterWavelengths[entry]; ++entry) {
			double a = WaterWavelengths[entry],
				b = WaterWavelengths[entry + 1],
				ca = max(a, start),
				cb = min(b, end),
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
	return complex<double>(waterN, waterK);
}


complex<double> GoldIOR(double wavelength) {
	double lambdaMin = wavelength - 0.5, lambdaMax = wavelength + 0.5;
	double start = glm::max(lambdaMin, GoldWavelengths[0]);
	double end = glm::min(lambdaMax, GoldWavelengths[86]);

	int idx = int(glm::max(glm::distance(0.0, double(BinarySearch(0, 86, start, GoldWavelengths))), 1.0) - 1.0);

	if (end <= start) {
		return 0.0;
	}

	double goldN = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 86 && end >= GoldWavelengths[entry]; ++entry) {
			double a = GoldWavelengths[entry],
				b = GoldWavelengths[entry + 1],
				ca = max(a, start),
				cb = min(b, end),
				fa = GoldN[entry],
				fb = GoldN[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		goldN = result / (lambdaMax - lambdaMin);
	}

	double goldK = 0.0;
	{
		double result = 0.0;
		int entry = idx;
		for (; entry < 86 && end >= GoldWavelengths[entry]; ++entry) {
			double a = GoldWavelengths[entry],
				b = GoldWavelengths[entry + 1],
				ca = max(a, start),
				cb = min(b, end),
				fa = GoldK[entry],
				fb = GoldK[entry + 1],
				invAB = 1.0 / (b - a);

			if (ca >= cb) {
				continue;
			}

			double interpA = glm::mix(fa, fb, (ca - a) * invAB);
			double interpB = glm::mix(fa, fb, (cb - a) * invAB);

			result += 0.5 * (interpA + interpB) * (cb - ca);
		}

		goldK = result / (lambdaMax - lambdaMin);
	}
	return complex<double>(goldN, goldK);
}