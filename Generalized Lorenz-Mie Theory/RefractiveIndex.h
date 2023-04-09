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

const std::complex<double> VacuumIOR = std::complex<double>(1.0, 0.0);

std::complex<double> BrineIOR(double wavelength);
std::complex<double> WaterIOR(double wavelength);
std::complex<double> IronIOR(double wavelength);
std::complex<double> CopperIOR(double wavelength);
std::complex<double> AirIOR(double wavelength);