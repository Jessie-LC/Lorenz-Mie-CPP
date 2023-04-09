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

double InterpAlgaeK(double wavelength);
double InterpMineralK(double wavelength);

std::complex<double> BrineIOR(double wavelength);
std::complex<double> WaterIOR(double wavelength);