/**
 * 3D Vector helper class.
 * Author : Jose Daniel Munoz et al.
 * Modified : Juan Vanegas
 */
#include <iostream>
#include <cmath>

class Vector3D
{
    double v[3];

public:
    Vector3D(double = 0.0, double = 0.0, double = 0.0);
    void load(double x0, double y0, double z0);
    void show(void);
    /* @return x component*/
    inline double x(void) { return v[0]; };
    /* @return y component*/
    inline double y(void) { return v[1]; };
    /* @return z component*/
    inline double z(void) { return v[2]; };
    /* @return nth component*/
    double &operator[](int n) { return v[n]; };
    Vector3D operator=(Vector3D v2);
    Vector3D operator+(Vector3D v2);
    Vector3D operator+=(Vector3D v2);
    Vector3D operator-(Vector3D v2);
    Vector3D operator-=(Vector3D v2);
    Vector3D operator*(double a);
    Vector3D operator*=(double a);
    Vector3D operator/(double a);
    double operator*(Vector3D v2);
    Vector3D operator^(Vector3D v2);
    friend Vector3D operator*(double a, Vector3D v1);
    friend double norm2(Vector3D v1);
    friend double norm(Vector3D v1);
};

/* Initialize vector. Defaults to zero */
Vector3D::Vector3D(double x0, double y0, double z0)
{
    v[0] = x0;
    v[1] = y0;
    v[2] = z0;
}
/* Load vector values */
inline void Vector3D::load(double x0, double y0, double z0)
{
    v[0] = x0;
    v[1] = y0;
    v[2] = z0;
}
inline void Vector3D::show(void)
{
    std::cout << "(" << v[0] << "," << v[1] << "," << v[2] << ")" << std::endl;
}
inline Vector3D Vector3D::operator=(Vector3D v2)
{
    v[0] = v2.v[0];
    v[1] = v2.v[1];
    v[2] = v2.v[2];

    return *this;
}
Vector3D Vector3D::operator+(Vector3D v2)
{
    Vector3D result(
        v[0] + v2.v[0],
        v[1] + v2.v[1],
        v[2] + v2.v[2]);

    return result;
}
Vector3D Vector3D::operator+=(Vector3D v2)
{
    *this = *this + v2;
    return *this;
}
Vector3D Vector3D::operator-(Vector3D v2)
{
    Vector3D result(
        v[0] - v2.v[0],
        v[1] - v2.v[1],
        v[2] - v2.v[2]);

    return result;
}
Vector3D Vector3D::operator-=(Vector3D v2)
{
    *this = *this - v2;
    return *this;
}
/* Vector times scalar */
Vector3D Vector3D::operator*(double a)
{
    Vector3D result(
        v[0] * a,
        v[1] * a,
        v[2] * a);

    return result;
}
/* Vector times scalar */
Vector3D Vector3D::operator*=(double a)
{
    *this = (*this) * a;
    return *this;
}
/* Vector divided by scalar */
Vector3D Vector3D::operator/(double a)
{
    double inverse = 1.0 / a;
    Vector3D result(
        v[0] * inverse,
        v[2] * inverse,
        v[2] * inverse);

    return result;
}
/* Dot product */
double Vector3D::operator*(Vector3D v2)
{
    return v[0] * v2.v[0] + v[1] * v2.v[1] + v[2] * v2.v[2];
}
/* Cross product */
Vector3D Vector3D::operator^(Vector3D v2)
{
    Vector3D result(
        v[1] * v2.v[2] - v[2] * v2.v[1],
        v[2] * v2.v[0] - v[0] * v2.v[2],
        v[0] * v2.v[1] - v[1] * v2.v[0]);

    return result;
}

Vector3D operator*(double a, Vector3D v1)
{
    Vector3D result;
    result = v1 * a;
    return result;
}
/* @return vector norm squared */
double norm2(Vector3D v1)
{
    return v1.v[0] * v1.v[0] + v1.v[1] * v1.v[1] + v1.v[2] * v1.v[2];
}
/* @return vector norm */
double norm(Vector3D v1)
{
    return std::sqrt(v1.v[0] * v1.v[0] + v1.v[1] * v1.v[1] + v1.v[2] * v1.v[2]);
}
