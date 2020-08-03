/**
 * Animation module.
 * For calculating and saving data please use main.cpp
 */
#include "molecular_dynamics.h"
#include "animation.hpp"

int main(int argc, char *argv[])
{
    double T = 0.1, Q = 0.1, tdibujo;
    Body Molecule[N];
    MolecularDynamics LennardJones;
    CRandom ran64(seed);
    int tSteps = (int)(tMax / dt);

    start_animation(argc);

    double vMax = T;
    HeatBath NoseHoover(0, 0, T, Q);

    // Initial configuration
    LennardJones.init(ran64, vMax, Molecule);
    LennardJones.calculate_all_forces(Molecule);

    tdibujo = 0;
    for (int q = 0; q < tSteps; q++)
    {

        if (tdibujo > 0.1)
        {
            begin_frame(argc);
            for (int i = 0; i < N; i++)
                Molecule[i].print();
            end_frame(argc);
            tdibujo = 0.0;
        }
        tdibujo += dt;
        LennardJones.thermalRESPA(Molecule, NoseHoover, dt);
    }

    return 0;
}