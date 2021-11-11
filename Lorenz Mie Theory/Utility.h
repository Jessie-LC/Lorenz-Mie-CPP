#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include <valarray>
#include "glm.hpp"
#include "constants.h"

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const double tau = pi * 2.0;

const double epsilon = 1e-5;

double Pn(double x, int n);

double derivativeLegendre(int n, double x);

double computePi(double mu, int n);

double derivativePi(double x, int n);

double computeTau(double mu, int n);

std::complex<double> SphJn(int n, std::complex<double> z);
std::complex<double> SphYn(int n, std::complex<double> z);

double sqr(double x);

glm::dvec3 xyzToRGB(glm::dvec3 xyz);

glm::dvec3 SpectrumToXYZ(glm::dvec3 spectrum, float w);