#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include "glm.hpp"

using namespace std;

double ComputeMiePhase(complex<double> iorHost, complex<double> iorParticle, double theta, double r, double lambda, valarray<complex<double>>& S1, valarray<complex<double>>& S2, double& Qabs, double& Qsca, double& Qext, double& sum, int angle);