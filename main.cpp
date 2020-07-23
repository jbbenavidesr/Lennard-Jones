/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include "molecular_dynamics.h"
#include "random64.h"

int main(int argc, char *argv[])
{
    Body Molecule[N];
    CRandom ran64(1);
    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);
    double dz = Lz / (Nz + 1);
    int tSteps = (int)(tMax / dt);

    double x, y, z, vx, vy, vz, vmax;
    vmax = 2.7;

    // Initial configuration
    for (int k = 0; k < N; k++) // Run through every molecule
    {
        // Initial positions in cubic lattice
        x = dx + (k % Nx) * dx;
        y = dy + ((k / Nx) % Ny) * dy;
        z = dz + (k / (Nx * Ny)) * dz;

        // Initial random velocities
        vx = vmax * (2 * ran64.r() - 1);
        vy = vmax * (2 * ran64.r() - 1);
        vz = vmax * (2 * ran64.r() - 1);

        Molecule[k].init(x, y, z, vx, vy, vz, m0);
    }

    // Reset center of mass velocity
    Vector3D vCM;
    for (int k = 0; k < N; k++)
    {
        vCM += Molecule[k].getV();
    }
    vCM = vCM / N;
    Vector3D resetV;
    for (int k = 0; k < N; k++)
    {
        resetV = Molecule[k].getV() - vCM;
        Molecule[k].setV(resetV);
    }

    std::cout << "t,x,y,z,vx,vy,vz,u\n";
    for (int n = 0; n < tSteps; n++)
    {
        // ----- Get forces -----
        // Every force starts in 0
        for (int k = 0; k < N; k++)
        {
            Molecule[k].resetForce();

            // Print state
            Molecule[k].printState(n * dt);
        }

        // calculate forces
        for (int i = 0; i < (N - 1); i++)
            for (int j = i + 1; j < N; j++)
            {
                // Calculate rMin
                Vector3D dr;
                dr = rMin(Molecule[i], Molecule[j]);

                //Calculate forces
                double normR = norm(dr);
                if (normR < rCut)
                {
                    double f = 24 * (2 * std::pow(normR, -14) - std::pow(normR, -8));
                    Molecule[i].addForce(f * dr);
                    Molecule[j].addForce((-1 * f) * dr);
                }
            }

        for (int k = 0; k < N; k++)
        {
            // Update velocities
            Molecule[k].moveV(dt);

            // Update Positions
            Molecule[k].moveR(dt);
        }
    }

    return 0;
}