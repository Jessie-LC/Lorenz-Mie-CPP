#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <future>
#include <thread>
#include <complex>
#include "glm.hpp"

std::complex<double> WaterIOR(double wavelength);
std::complex<double> GoldIOR(double wavelength);