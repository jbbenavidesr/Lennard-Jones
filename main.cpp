/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include <string>
#include <sstream>
#include <fstream>

#include "molecular_dynamics.h"

std::string filename(int n);

int main(int argc, char *argv[])
{
    int max = 8 * 20; // number of cores * number of temperatures per core. Change this accordingly
    double T_0 = 0.1;
    double Tmax = 1.1, dT = Tmax / (double)max;

#pragma omp parallel for shared(T_0, Tmax, dT)
    for (int n = 0; n < max; n++)
    {
        double T = T_0 + dT * n;

        Body Molecule[N];

        MolecularDynamics LennardJones;
        CRandom ran64(seed);
        int tSteps = (int)(tMax / dt);
        double vMax;
        double L, xi0, U, m2, m4, auxM2, auxM4;
        double Q = 0.1;

        std::ofstream File(filename(n));

        File << "t,L\n";
        //File << "t,x,y,z,vx,vy,vz,u\n";
        //m2 = 0;
        //m4 = 0;
        L = 0;
        vMax = T;
        HeatBath NoseHoover(0, 0, T, Q);

        // Initial configuration
        LennardJones.init(ran64, vMax, Molecule);
        LennardJones.calculate_all_forces(Molecule);

        for (int q = 0; q < tSteps; q++)
        {
            File << dt * q << ',' << LennardJones.getL() << std::endl;
            /*
            for (int i = 0; i < N; i++)
            {
                Molecule[i].printState(q * dt, File);
                //std::cout << dt * n << ',' << LennardJones.T(Molecule) << std::endl;
            }
            */
            //if (q * dt > tEq)
            //{
            //L += LennardJones.getL();
            //LennardJones.velocityVerletStep(Molecule, dt);
            //std::tie(auxM2, auxM4) = LennardJones.Binder(Molecule, 16);
            //m2 += auxM2;
            //m4 += auxM4;
            //}
            //else
            //{
            LennardJones.thermalRESPA(Molecule, NoseHoover, dt);
            //}
        }
        //L /= (tSteps - (tEq / dt));
        //m2 /= (tSteps - (tEq / dt));
        //m4 /= (tSteps - (tEq / dt));
        // U = m4 / (m2 * m2);
        //File << T << ',' << L << std::endl;
        File.close();
    }

    return 0;
}

std::string filename(int n)
{
    std::string name;
    std::stringstream n_s;
    n_s << n;
    name = "./" + n_s.str() + ".csv";
    return name;
}
