#include "constants.h"
#include "vector3D.h"
#include "random64.h"

// ----- Class Declarations -----
class HeatBath
{
};

class Body
{
private:
    // It has a position r, velocity v and a total force F acts on it.
    Vector3D r, v, F;
    // It has a mass and a total potential energy U.
    double m, U;
    // For this first implementation I'll ignore the geometry and assume it as a point mass.

public:
    void init(double x, double y, double z, double vx, double vy, double vz, double mass = 1.0);
    void printState(double t);

    // --Inline Functions --

    // Makes the total force equal to 0 for recalculation.
    inline void resetForce(void) { F.load(0, 0, 0); };

    // Adds a force to the total force on the body
    inline void addForce(Vector3D dF) { F += dF; };

    // Makes the total potential equal to 0 for recalculation.
    inline void resetU(void) { U = 0; };

    // Adds a potential to the total  Potential energy of the body
    inline void addU(double dU) { U += dU; };

    // update of velocity
    inline void moveV(double dt, double coef = 1.0) { v += (dt * coef / m) * F; };

    // update of position
    inline void moveR(double dt, double coef = 1.0) { r += (dt * coef) * v; };

    // Get position
    inline Vector3D getR(void) { return r; };

    // Get velocity
    inline Vector3D getV(void) { return v; };

    // Set position
    inline void setR(Vector3D newR) { r = newR; };

    // Set velocity
    inline void setV(Vector3D newV) { v = newV; };

    //-- Friends --

    friend Vector3D rMin(Body b1, Body b2);
    friend class MolecularDynamics;
    friend class HeatBath;
};

class MolecularDynamics
{
public:
    void init(CRandom &ran64, double vMax, Body *Molecule);
    void calculate_force_pair(Body &molecule1, Body &molecule2);
    void calculate_all_forces(Body *molecule);
    void leapFrogStep(Body *molecule, double dt);
    void leapFrogThermalStep(Body *molecule, HeatBath &thermo, double dt);

    void velocityVerletStep(Body *molecule, double dt);
    void velocityVerletThermalStep(Body *molecule, HeatBath &thermo, double dt);

    //void move_with_pefrl(Body *molecule, double dt);
    Vector3D forceLJ(Vector3D dr);
    double T(Body *molecule);
};

// ----- Function implementations -----

/*
 * Initializes the molecule in a given initial position, with a given velocity and a given mass
 * Force and potential start being 0.
 * 
 * @param double x, y, z: Initial position of the body.
 * @param double vx, vy, vz: Initial velocity of the body.
 * @param double m: Mass of the body. Defaults to 1.
 */
void Body::init(double x, double y, double z, double vx, double vy, double vz, double mass)
{
    r.load(x, y, z);
    v.load(vx, vy, vz);
    m = mass;
    resetForce();
    resetU();
}

/*
 * Prints the current state of the body in csv format. For the moment it uses cout
 * a future implementation is set to make it write in a given file directly.
 * 
 * @param double t: The current timestep in the simulation
 */
void Body::printState(double t)
{
    std::cout << t << "," << r.x() << "," << r.y() << "," << r.z() << "," << v.x() << "," << v.y() << "," << v.z() << "," << U << std::endl;
}

/*
 * Initializes a MD system with lattice arrangement and random velocities within a maximum value per component.
 * 
 * @param CRandom ran64: Random number generator for initial velocities.
 * @param double vmax: Each velocity component starts somewhere in (-vmax, vmax)
 */
void MolecularDynamics::init(CRandom &ran64, double vmax, Body *Molecule)
{
    double x, y, z, vx, vy, vz;

    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);
    double dz = Lz / (Nz + 1);

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
            Molecule[k].v -= vCM;
        }
    }
}

/*
 * Calculates force between a pair of molecules. 
 * Force between molecules is calculated using a Lennard-Jones force.  
 * 
 * @params molecule1, molecule2 : Pair of molecules for which the force is to be calculated.
 */
void MolecularDynamics::calculate_force_pair(Body &molecule1, Body &molecule2)
{
    // Calculate rMin
    Vector3D dr;
    dr = rMin(molecule1, molecule2);

    if (norm(dr) < rCut)
    {
        Vector3D f = forceLJ(dr);
        molecule1.addForce(f);
        molecule2.addForce((-1) * f);
    }
}

/**
 * Calculates force between all molecules in the system. 
 * 
 * @params molecule: Array of Molecules in the system.
 */
void MolecularDynamics::calculate_all_forces(Body *molecule)
{
    for (int i = 0; i < N; i++)
        molecule[i].resetForce();

    for (int i = 0; i < (N - 1); i++)
        for (int j = i + 1; j < N; j++)
            calculate_force_pair(molecule[i], molecule[j]);
}

/*
 * Evolves one timestep of the system using a second order LeapFrog algorithm.
 * 
 * @param Body molecule: Array with all the particles of the system.
 * @param double dt: Time step of the simulation
 */
void MolecularDynamics::leapFrogStep(Body *molecule, double dt)
{
    calculate_all_forces(molecule);
    for (int k = 0; k < N; k++)
    {
        // Update velocities
        molecule[k].moveV(dt);

        // Update Positions
        molecule[k].moveR(dt);
    }
}

/*
 * Evolves one timestep of the system using a second order velocity verlet algorithm.
 * 
 * @param Body molecule: Array with all the particles of the system.
 * @param double dt: Time step of the simulation
 */
void MolecularDynamics::velocityVerletStep(Body *molecule, double dt)
{

    for (int k = 0; k < N; k++)
    {
        // Update velocities to t + dt/2
        molecule[k].moveV(dt, 0.5);

        // Update Positions to t +dt
        molecule[k].moveR(dt);
    }
    // get force in t + dt
    calculate_all_forces(molecule);

    // Update velocities to t + dt
    for (int k = 0; k < N; k++)
        molecule[k].moveV(dt, 0.5);
}

/*
 * Calculates force using the Lennard-Jones potential for a couple of molecules.
 * 
 * @param Vector3D dr: vector from one molecule to the other.
 * 
 * @return Vector3D Force vector between these molecules.
*/
Vector3D MolecularDynamics::forceLJ(Vector3D dr)
{
    double r = norm(dr);
    return (24 * (2 * std::pow(r, -14) - std::pow(r, -8))) * dr;
}

/*
 * Calculates the Temperature of the system in a given moment
 * 
 * @param Body molecule: Array with all the particles of the system.
 */
double MolecularDynamics::T(Body *molecule)
{
    double T = 0;
    for (int i = 0; i < N; i++)
    {
        T += norm2(molecule[i].getV());
    }
    T *= m0 / (3 * N);
    return T;
}

/*
 * Calculates the minimum image distance of 2 bodies. This is based in the image
 * strategy for implementing periodic boundaries.
 * 
 * @param Body b1, b2: the 2 bodies whose minimum distance we want to calculate.
 * 
 * @return Vector3D: Returns the minimum distance between the bodies.
 */
Vector3D rMin(Body b1, Body b2)
{
    Vector3D dr;
    dr = b1.r - b2.r;
    double x, y, z;
    x = dr.x() - Lx * std::round(dr.x() / Lx);
    y = dr.y() - Ly * std::round(dr.y() / Ly);
    z = dr.z() - Lz * std::round(dr.z() / Lz);
    Vector3D rmin(x, y, z);
    return rmin;
}
