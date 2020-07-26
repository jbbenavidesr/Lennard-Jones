/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include "molecular_dynamics.h"

int main(int argc, char *argv[])
{
    Body Molecule[N];
    MolecularDynamics LennardJones;
    CRandom ran64(seed);
    int tSteps = (int)(tMax / dt);
    double vMax = 0.01;

    // Initial configuration
    LennardJones.init(ran64, vMax, Molecule);
    LennardJones.calculate_all_forces(Molecule);

    std::cout << "t,x,y,z,vx,vy,vz,u\n";
    for (int n = 0; n < tSteps; n++)
    {
        for (int i = 0; i < N; i++)
        {
            Molecule[i].printState(n * dt);
            //std::cout << dt * n << ',' << LennardJones.T(Molecule) << std::endl;
        }

        LennardJones.velocityVerletStep(Molecule, dt);
    }

    return 0;
}