#include <iostream>
#include <cmath>
#include <fstream>

#include "random64.h"
#include "vector3D.h"

const int Nline = 10;
const int N = Nline * Nline * Nline;
const double L = 10.0; // V= L*L*L
const double T = 10.0;
const double dt = 0.001;
const double rCut = 3.0;
const double m = 1.0;
const double TB = 4.0;

double lennard_jones(double r)
{
    return 24 * (2 * std::pow(r, -14) - std::pow(r, -8));
}

double potential(double r)
{
    return 4 * (std::pow(r, -12) - std::pow(r, -6));
}

Vector3D r_min(Vector3D r1, Vector3D r2)
{
    Vector3D r;
    r = r1 - r2;
    double x, y, z;
    x = r.x() - L * std::round(r.x() / L);
    y = r.y() - L * std::round(r.y() / L);
    z = r.z() - L * std::round(r.z() / L);
    Vector3D rmin(x, y, z);
    return rmin;
}

/* 
 * Calculates all forces between molecules in a given time step.
*/
void force(Vector3D *F, Vector3D *r, double *U)
{
    // Every force starts in 0
    for (int k = 0; k < N; k++)
    {
        F[k].load(0, 0, 0);
        U[k] = 0;
    }

    // calculate forces
    for (int i = 0; i < (N - 1); i++)
        for (int j = i + 1; j < N; j++)
        {
            Vector3D dr;
            dr = r_min(r[i], r[j]);
            if (norm(dr) < rCut)
            {
                double f = lennard_jones(norm(dr));
                F[i] += f * dr;
                F[j] -= f * dr;
                // std::cout << norm(dr) << " " << f * norm(dr) << std::endl;
                U[i] += potential(norm(dr));
                U[j] += potential(norm(dr));
            }
        }
}

void moveV(Vector3D *v, Vector3D *vnew, Vector3D *F)
{
    for (int k = 0; k < N; k++)
    {
        vnew[k] = v[k] + (dt / m) * F[k];
    }
}

void moveR(Vector3D *r, Vector3D *rnew, Vector3D *vnew)
{
    for (int k = 0; k < N; k++)
    {
        rnew[k] = r[k] + dt * vnew[k];
    }
}

double Temp(Vector3D *V)
{
    double t = 0;
    for (int k = 0; k < N; k++)
    {
        t += norm2(V[k]);
    }
    t /= 3 * N;
    return t;
}

void thermostat(Vector3D *V, double tau)
{
    double l;
    double T = Temp(V);
    l = 1 + (dt / tau) * ((TB / T) - 1);
    for (int k = 0; k < N; k++)
        V[k] *= l;
}

int main(void)
{
    CRandom ran64(1);
    double l = L / (Nline + 1);
    int Nsteps = (int)(T / dt);
    double x, y, z, vx, vy, vz;
    double Teq = 10.0;
    double tau = 5.0;

    /* Arrays for positions and velocities and force:
     * in positions r is for t and rnew is t+dt
     * in velocities v is for t-(dt/2) and vnew for t+(dt/2)
     * force is only for t and F[k] is the total force on particle k    
    */
    Vector3D r[N], rnew[N], v[N], vnew[N], F[N];
    double U[N];

    // Initial configuration
    for (int k = 0; k < N; k++) // Run through every molecule
    {
        // Initial positions
        x = l + (k % Nline) * l;
        y = l + ((k / Nline) % Nline) * l;
        z = l + (k / (Nline * Nline)) * l;
        r[k].load(x, y, z);

        // Initial random velocities
        double vmax = 2.7;
        vx = vmax * (2 * ran64.r() - 1);
        vy = vmax * (2 * ran64.r() - 1);
        vz = vmax * (2 * ran64.r() - 1);
        v[k].load(vx, vy, vz);
    }

    // Reset center of mass velocity
    Vector3D vCM;
    for (int k = 0; k < N; k++)
    {
        vCM += v[k];
    }
    vCM = vCM / N;
    for (int k = 0; k < N; k++)
    {
        v[k] -= vCM;
    }
    std::cout << "t,x,y,z,vx,vy,vz,u\n";
    for (int n = 0; n < Nsteps; n++)
    {
        force(F, r, U);
        moveV(v, vnew, F);
        moveR(r, rnew, vnew);

        for (int k = 0; k < N; k++)
        {
            std::cout << n * dt << "," << r[k].x() << "," << r[k].y() << "," << r[k].z() << "," << v[k].x() << "," << v[k].y() << "," << v[k].z() << "," << U[k] << std::endl;
            v[k] = vnew[k];
            r[k] = rnew[k];
        }
        /*
        if (n * dt > Teq)
            thermostat(v, tau);

        if (n * dt > Teq + 10)
            tau = 1000;
        */
    }

    return 0;
}