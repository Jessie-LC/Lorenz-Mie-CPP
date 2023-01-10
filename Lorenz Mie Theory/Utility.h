#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include <valarray>
#include <chrono>
#include "glm.hpp"
#include "constants.h"

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const double tau = pi * 2.0;

const double epsilon = 1e-5;

struct LegendreDerivatives {
	double value;
	double derivative1;
	double derivative2;
};

void LegendreBulk(double x, int n, double out0[], double out1[], double out2[]);
LegendreDerivatives LegendreAll(double x, int n);

double Factorial(int n);
double T_Gamma(double z);
std::complex<double> SphJn(int n, std::complex<double> z);
std::complex<double> SphYn(int n, std::complex<double> z);
glm::vec3 SpectrumToXYZ(float spectrum, float w);