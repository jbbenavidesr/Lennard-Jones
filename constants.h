#include <iostream>
#include <cmath>
#include <fstream>

// Geometry constants
const double Lx = 20.0, Ly = 20.0, Lz = 50.0;
const double V = Lx * Ly * Lz;

// PEFRL
const double Zi = 0.1786178958448091e0;
const double Lambda = 0.2123418310626054 * (-1);
const double Xi = 0.06626458266981849 * (-1);

const double coef1 = (1 - 2 * Lambda) / 2;
const double coef2 = (1 - 2 * (Xi + Zi));

// Implementation
const int Nx = 4, Ny = 4, Nz = 4;
const int N = Nx * Ny * Nz;
const double m0 = 1.0;
const double tMax = 100.0;
const double tEq = 25.0;
const double dt = 0.001;
const double rCut = 3.0;
const double TB = 1;

// Random number generator
const unsigned long long seed = 1;
