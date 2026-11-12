/*
Definitions of 3D vectors, 3x3 matrices, and quaternions with basic operations used in rigid 
body dynamics and control. These definitions are used in this project.

KEY CONVENTIONS USED:
1. Torque about the CoG of a force F applied at point p: T = cross(p-cg, F)
*/

#pragma once
#include <stdexcept>
#include <cmath>

// 3D vector definition
struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    double operator[](int i)      const { return i == 0 ? x : (i == 1 ? y : z); }
    double& operator[](int i)           { return i == 0 ? x : (i == 1 ? y : z); }
    Vec3 operator+(const Vec3& b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s)      const { return {x / s, y / s, z / s}; }
    Vec3 operator-()             const { return {-x, -y, -z}; }
    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
};

inline Vec3 operator*(double s, const Vec3& v) { return {v.x * s, v.y * s, v.z * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) { 
    double n = norm(v);
    return n > 0 ? v / n : v; 
}

// 3x3 matrix definition
struct Mat3 {
    double matrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    static Mat3 identity() {
        Mat3 I;
        I.matrix[0][0] = 1;
        I.matrix[1][1] = 1;
        I.matrix[2][2] = 1;
        return I;
    }
    static Mat3 diagonal(double a, double b, double c) {
        Mat3 D;
        D.matrix[0][0] = a;
        D.matrix[1][1] = b;
        D.matrix[2][2] = c;
        return D;
    }
    Vec3 operator*(const Vec3& v) const {
        return {
            matrix[0][0] * v.x + matrix[0][1] * v.y + matrix[0][2] * v.z,
            matrix[1][0] * v.x + matrix[1][1] * v.y + matrix[1][2] * v.z,
            matrix[2][0] * v.x + matrix[2][1] * v.y + matrix[2][2] * v.z
        };
    }
    Mat3 operator*(const Mat3& b) const {
        Mat3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.matrix[i][j] = matrix[i][0] * b.matrix[0][j] +
                                      matrix[i][1] * b.matrix[1][j] +
                                      matrix[i][2] * b.matrix[2][j];
            }
        }
        return result;
    }
    Mat3 transpose() const {
        Mat3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.matrix[i][j] = matrix[j][i];
            }
        }
        return result;
    }
    double determinant() const {
        return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
               matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
               matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    }
    Mat3 inverse() const {
        double det = determinant();
        if (det == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }
        Mat3 inv;
        inv.matrix[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
        inv.matrix[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det;
        inv.matrix[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;
        inv.matrix[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det;
        inv.matrix[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
        inv.matrix[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det;
        inv.matrix[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
        inv.matrix[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det;
        inv.matrix[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;
        return inv;
    }
};
