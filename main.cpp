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
    double vMax = 0.5;
    double L, xi0, U, m2, m4, auxM2, auxM4;
    double Q = 0.2;

    std::cout << "T,L\n";

    for (double T = 1; T < 1.2; T += 1)
    {
        //m2 = 0;
        //m4 = 0;
        L = 0;
        vMax = T;
        xi0 = T / Q;
        HeatBath NoseHoover(xi0, T, Q);

        // Initial configuration
        LennardJones.init(ran64, vMax, Molecule);
        LennardJones.calculate_all_forces(Molecule);

        for (int n = 0; n < tSteps; n++)
        {
            //std::cout << dt * n << ',' << LennardJones.getL() << std::endl;
            /*
            for (int i = 0; i < N; i++)
            {
                Molecule[i].printState(n * dt);
                //std::cout << dt * n << ',' << LennardJones.T(Molecule) << std::endl;
            }
            */
            LennardJones.velocityVerletThermalStep(Molecule, NoseHoover, dt);
            if (n * dt > tEq)
            {
                L += LennardJones.getL();
                //std::tie(auxM2, auxM4) = LennardJones.Binder(Molecule, 16);
                //m2 += auxM2;
                //m4 += auxM4;
            }
        }
        L /= (tSteps - (tEq / dt));
        //m2 /= (tSteps - (tEq / dt));
        //m4 /= (tSteps - (tEq / dt));
        // U = m4 / (m2 * m2);
        std::cout << T << ',' << L << std::endl;
    }

    // std::cout << "t,x,y,z,vx,vy,vz,u\n";

    return 0;
}