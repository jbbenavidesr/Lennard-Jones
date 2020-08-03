#include <tuple>

#include "constants.h"
#include "vector3D.h"
#include "random64.h"

// ----- Class Declarations -----
class HeatBath
{
private:
    double p, eta, Q, T;

public:
    HeatBath(double eta0, double p0, double Tb = 1.0, double Q0 = 1.0);
    void moveP(double K, double dt);

    inline void moveEta(double dt) { eta += dt * p / Q; };
    inline double getEta(void) { return eta; };
    inline double getP(void) { return p; };
    inline double getQ(void) { return Q; };
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
    void printState(double t, std::ofstream &File);

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

    //Get kinetic energy
    inline double getK() { return 0.5 * m * norm2(v); };

    // print gnuplot animation
    inline void print() { std::cout << "," << r.x() << "+0.2*cos(t)," << r.z() << "+0.2*sin(t)"; };

    //-- Friends --

    friend Vector3D rMin(Body b1, Body b2);
    friend class MolecularDynamics;
    friend class HeatBath;
};

class MolecularDynamics
{
private:
    double L = 0;

public:
    void init(CRandom &ran64, double vMax, Body *Molecule);
    void calculate_force_pair(Body &molecule1, Body &molecule2);
    void calculate_all_forces(Body *molecule);
    void leapFrogStep(Body *molecule, double dt);
    void leapFrogThermalStep(Body *molecule, HeatBath &thermo, double dt);

    void velocityVerletStep(Body *molecule, double dt);
    void thermalRESPA(Body *molecule, HeatBath &thermo, double dt);

    //void move_with_pefrl(Body *molecule, double dt);
    double forceLJ(double r2);
    double potentialLJ(double r2);
    double T(Body *molecule);
    double K(Body *molecule);
    std::tuple<double, double> Binder(Body *molecule, int Mb);

    inline double getL() { return L; }
};

// ----- Function implementations -----

/*
 * Initializes Heat bath. 
 */
HeatBath::HeatBath(double eta0, double p0, double Tb, double Q0)
{
    eta = eta0;
    p = p0;
    Q = Q0;
    T = Tb;
}

/*
 * Updates the value of Ps
 */
void HeatBath::moveP(double K, double dt)
{
    p += dt * (2 * K - (3 * N * T));
}

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
void Body::printState(double t, std::ofstream &File)
{
    File << t << "," << r.x() << "," << r.y() << "," << r.z() << "," << v.x() << "," << v.y() << "," << v.z() << "," << U << '\n';
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
    double dz = 0.5 * Lz / (Nz + 1);

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

    double norm2_dr = norm2(dr);
    L += std::sqrt(norm2_dr);

    if (norm2_dr < (rCut * rCut))
    {
        Vector3D f = forceLJ(norm2_dr) * dr;
        double u = potentialLJ(norm2_dr);
        molecule1.addForce(f);
        molecule1.addU(u);
        molecule2.addForce((-1) * f);
        molecule2.addU(u);
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
    {
        molecule[i].resetForce();
        molecule[i].resetU();
    }

    L = 0;

    for (int i = 0; i < (N - 1); i++)
        for (int j = i + 1; j < N; j++)
            calculate_force_pair(molecule[i], molecule[j]);

    L /= 0.5 * N * (N - 1);
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
 * Evolves one timestep of the system using a second order velocity verlet algorithm.
 * 
 * @param Body molecule: Array with all the particles of the system.
 * @param HeatBath thermo: NoseHoover thermal bath.
 * @param double dt: Time step of the simulation
 */
void MolecularDynamics::thermalRESPA(Body *molecule, HeatBath &thermo, double dt)
{
    double kin, dtHalf, factor;
    kin = K(molecule);
    dtHalf = 0.5 * dt;
    thermo.moveP(kin, dtHalf); // 1
    thermo.moveEta(dtHalf);    // 2
    factor = std::exp((-1) * thermo.getP() * dtHalf / thermo.getQ());
    for (int i = 0; i < N; i++)
    {
        molecule[i].setV(factor * molecule[i].getV()); // 3
        molecule[i].moveV(dtHalf);                     // 4
        molecule[i].moveR(dt);                         // 5
    }

    calculate_all_forces(molecule); // 6

    for (int i = 0; i < N; i++)
    {
        molecule[i].moveV(dtHalf);                     // 7
        molecule[i].setV(factor * molecule[i].getV()); // 8
    }
    kin = K(molecule);
    thermo.moveEta(dtHalf);    // 9
    thermo.moveP(kin, dtHalf); // 10
}

/*
 * Calculates force using the Lennard-Jones potential for a couple of molecules.
 * 
 * @param double r2: distance squared from one molecule to the other.
 * 
 * @return double Force between these molecules.
*/
double MolecularDynamics::forceLJ(double r2)
{
    double dr6 = r2 * r2 * r2;
    double dr8 = dr6 * r2;
    dr6 = 1.0 / dr6;
    dr8 = 1.0 / dr8;

    return 24 * dr8 * (2 * dr6 - 1.0) - 8 * c2;
}

/*
 * Calculates potential using the Lennard-Jones potential for a couple of molecules.
 * 
 * @param double r2: distance squared from one molecule to the other.
 * 
 * @return double Potential Energy between these molecules.
*/
double MolecularDynamics::potentialLJ(double r2)
{
    double dr6 = r2 * r2 * r2;
    dr6 = 1.0 / dr6;

    return 4 * (dr6 * (dr6 - 1.0) + c2 * r2 + c0);
}

/*
 * Calculates the Kinetic Energy of the system in a given moment
 * 
 * @param Body molecule: Array with all the particles of the system.
 */
double MolecularDynamics::K(Body *molecule)
{
    double K = 0;
    for (int i = 0; i < N; i++)
    {
        K += molecule[i].getK();
    }
    return K;
}

/*
 * Calculates the Temperature of the system in a given moment
 * 
 * @param Body molecule: Array with all the particles of the system.
 */
double MolecularDynamics::T(Body *molecule)
{
    double T = K(molecule);
    T *= 2 / (3 * N);
    return T;
}

std::tuple<double, double> MolecularDynamics::Binder(Body *molecule, int Mb)
{
    Vector3D origin, r;
    double lx, ly, lz, rohMean, m2, m4, cVol;
    int nx, ny, nz, Mb3;
    lx = Lx / Mb;
    ly = Ly / Mb;
    lz = Lz / Mb;
    cVol = lx * ly * lz;
    m2 = 0;
    m4 = 0;
    Mb3 = Mb * Mb * Mb;
    rohMean = N / (cVol * Mb3);

    int nMol[Mb3];
    for (int i = 0; i < Mb3; i++)
        nMol[i] = 0;

    for (int k = 0; k < N; k++)
    {
        r = molecule[k].getR();
        nx = ((int)(r.x() / lx)) % Mb;
        ny = ((int)(r.y() / ly)) % Mb;
        nz = ((int)(r.z() / lz)) % Mb;
        nMol[nx + ny * Mb + nz * Mb * Mb] += 1;
    }

    for (int i = 0; i < Mb3; i++)
    {
        m2 += ((nMol[i] / cVol) - rohMean) * ((nMol[i] / cVol) - rohMean);
        m4 += std::pow(((nMol[i] / cVol) - rohMean), 4);
    }
    m2 /= Mb3;
    m4 /= Mb3;
    return std::make_tuple(m2, m4);
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
