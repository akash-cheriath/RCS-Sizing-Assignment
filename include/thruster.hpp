// Defines the RCS thruster layouts, configuration (B) matrix, control authority analysis and thruster allocation

#pragma once
#include <array>
#include <string>
#include <vector>
#include <stdexcept>
#include "config_map.hpp"
#include "const_def.hpp"
#include "simplex_solver.hpp"
#include "math_structs.hpp"

struct Thruster {
    std::string name;           // E.g. "P1_FWD_a"
    std::string type;           // E.g. AFT (+X-force), FWD (-X-force), ROLLP, ROLLN, ...
    Vec3 pos;                   // Position in body frame
    Vec3 dir;                   // Unit force direction in body frame
    double force_value;         // [N] Nominal vacuum thrust  
    int pod;                    // Pod index 
    int string_id;              // Redundancy string (0 = A, 1 = B)            
};

struct Layout {
    std::string name;
    std::vector<Thruster> thr;
    int size() const{ return static_cast<int>(thr.size()); }
};

struct MassProps {
    std::string name;
    double mass;                // [kg]
    Vec3 cog;                   // [m] centre of gravity in body frame 
    Mat3 inertia;               // [kg m^2] Intertia tensor about the CoG 
};

struct PodGeometry {
    double base_radius;         // [m] Heat shield max radius 
    double half_angle;          // [rad] Sidewall half-cone angle of the capsule
    double pod_x;               // [m] Axial position of the pod
    double clock0;              // [rad] Clocking of the first pod from +Y towards +Z
    double cant_aft;            // [rad] Outward cant of the aft-pointing nozzles w.r.t axial firing direction
    double cant_fwd;            // [rad] Outward cant of the forward-pointing nozzles w.r.t axial firing direction (also used by thrusters at nose cone)
    double cant_roll;           // [rad] Outward cant of the tangential nozzles w.r.t ideal firing direction
    double force;               // [N] Thrust per thruster
    double nose_x = 0.0;        // [m] Axial position of the nose
    double nose_r = 0.0;        // [m] Radius of the thruster arrangement on the nose
    double sideWallRadius(double x) const { return base_radius - (x * std::tan(half_angle)); }
};

// Layout Builders
// Thrusters sit on the conical sidewall of the capsule, arranged in pods, with local unit vectors at clock angle phi : radial r, tangential t

namespace layouts {

    inline void addThruster(Layout& Lay, const PodGeometry& podGeom, int pod, double phi,
                            const std::string& type, int string_id, const std::string& tag ="") {
        double sidewall_radius = podGeom.sideWallRadius(podGeom.pod_x);
        Vec3 pos(podGeom.pod_x, sidewall_radius * std::cos(phi), sidewall_radius * std::sin(phi));
        Vec3 x_hat(1, 0, 0);
        Vec3 r_hat(0, std::cos(phi), std::sin(phi));
        Vec3 t_hat(0, -std::sin(phi), std::cos(phi));
        Vec3 dir;

        if (type == "AFT")          // Exhaust aftwards and outwards -> force in +X-direction accelerating the capsule
            dir = (x_hat * std::cos(podGeom.cant_aft)) - (r_hat * std::sin(podGeom.cant_aft)); 
        else if (type == "FWD")     // Exhaust forwards and outwards -> force in -X-direction accelerating the capsule
            dir = (-x_hat * std::cos(podGeom.cant_fwd)) - (r_hat * std::sin(podGeom.cant_fwd));
        else if (type == "ROLLP")   // Exhaust tangentially (-ve) and outwards -> + Roll-torque
            dir = (t_hat * std::cos(podGeom.cant_roll)) - (r_hat * std::sin(podGeom.cant_roll));
        else if (type == "ROLLN")   // Exhaust tangentially (+ve) and outwards -> - Roll-torque
            dir = (-t_hat * std::cos(podGeom.cant_roll)) - (r_hat * std::sin(podGeom.cant_roll));
        
        Thruster thruster;
        thruster.name = "P" + std::to_string(pod + 1) + "_" + type + tag;
        thruster.type = type;
        thruster.pos = pos;
        thruster.dir = normalize(dir);
        thruster.force_value = podGeom.force;
        thruster.pod = pod;
        thruster.string_id = string_id;
        Lay.thr.push_back(thruster);
    }

    // POD LAYOUTS
    // Each pod carries nAft AFT, nFwd FWD and nRoll pairs of ROLLP/ROLLN thrusters
    // Redundancy strings : If a pod has 2 thrusters of same type, one goes into each string ; if it has only 1 thruster of any type,
    //                      opposite pods (e.g. 0,2) form string A and the other (e.g. 1,3) forms the string B - Each string holds
    //                      balanced opposite pairs

    // GENERIC Ring Pods at any clock angles [rad] (e.g., 3, 4, or 6 pods, also uneven spacing) :
    // Every pod carries nAft AFT, nFwd FWD and nROLLpairs of ROLLP/ROLLN
    // Redundancy allocated as mentioned above
    // alternateRoll = true (with nRoll = 1) : ONE roll jet per pod -> ROLLP on even and ROLLN on odd pods
    //                                       - Opposite pods still form pure roll couples (need even number of pods)

    inline Layout Ring(PodGeometry podGeom, const std::vector<double>& clocks, int nAft, int nFwd, int nRollPairs,
                       const std::string& name, bool alternateRoll = false) {
        
        Layout Lay;
        Lay.name = name;
        auto stringFor = [](int n, int k, int p) {
            return (n >=2) ? (k % 2) : (p % 2);
        };
        for (int p = 0; p < static_cast<int>(clocks.size()); ++p) {
            const double phi = clocks[p];
            for (int k = 0; k < nAft; ++k)
                addThruster(Lay, podGeom, p, phi, "AFT", stringFor(nAft, k, p), k ? "_b" : "_a");
            for (int k = 0; k < nFwd; ++k)
                addThruster(Lay, podGeom, p, phi, "FWD", stringFor(nFwd, k, p), k ? "_b" : "_a");
            if (alternateRoll) {
                addThruster(Lay, podGeom, p, phi, (p % 2 == 0) ? "ROLLP" : "ROLLN", (p / 2) % 2, "_a");
                continue;
            }
            for (int k = 0; k < nRollPairs; ++k) {
                addThruster(Lay, podGeom, p, phi, "ROLLP", stringFor(nRollPairs, k, p), k ? "_b" : "_a");
                addThruster(Lay, podGeom, p, phi, "ROLLN", stringFor(nRollPairs, k, p), k ? "_b" : "_a");
            }
        }
        return Lay;
    }
    
    // NOSE Ring (Dragon-like) : Forward firing thrusters on a ring of radius r around the docking port
    // Exhaust foward and cant_fwd outwards -> -X force (braking, deorbiting, CAM)
    // nPerPos indicates number of thrusters at a particular position. Pod positions on the nose numbered from podOffset

    inline void addNoseRing(Layout& Lay, const PodGeometry& podGeom, double nose_thr_x, double nose_thr_r, const std::vector<double>& clocks,
                            int nPerPos, int podOffset) {
        for (int p = 0; p < static_cast<int>(clocks.size()); ++p) {
            const double phi = clocks[p];
            const Vec3 x_hat(1, 0, 0);
            Vec3 r_hat(0, std::cos(phi), std::sin(phi));
            for (int k = 0; k < nPerPos; ++k) {
                Thruster thruster;
                thruster.name = "N" + std::to_string(p + 1) + "_NOSE" + (nPerPos > 1 ? (k ? "_b" : "_a") : "");
                thruster.type = "NOSE";
                thruster.pos = Vec3(nose_thr_x, nose_thr_r * std::cos(phi), nose_thr_r * std::sin(phi));
                thruster.dir = normalize((-x_hat * std::cos(podGeom.cant_fwd)) - (r_hat * std::sin(podGeom.cant_fwd)));
                thruster.force_value = podGeom.force;
                thruster.pod = podOffset + p;
                thruster.string_id = nPerPos >= 2 ? k % 2 : p % 2;
                Lay.thr.push_back(thruster);
            }
        }
    }
    
} // namespace layouts

// ------------------------------ CONFIGURATION MATRIX -----------------------------------------
// Configuration matrix B (6 x N) : column j = [F_j d_j ; (pos_j - cog) x F_j d_j]
// Rows 0-2 = force [N] in body-fixed frame, Rows 3-5 = Torque about the CoG [N m]
using Matrix = std::vector<std::vector<double>>;
using Vec6 = std::array<double, 6>;

inline Matrix buildB(const Layout& Lay, const Vec3& cog) {
    Matrix B(6, std::vector<double>(Lay.size(), 0.0));
    for (int j = 0; j < Lay.size(); ++j) {
        Vec3 force = Lay.thr[j].dir * Lay.thr[j].force_value;
        Vec3 torque = cross(Lay.thr[j].pos - cog, force);
        B[0][j] = force.x;
        B[1][j] = force.y;
        B[2][j] = force.z;
        B[3][j] = torque.x;
        B[4][j] = torque.y;
        B[5][j] = torque.z;
    }
    return B;
}

// -------------------------- CONTROL AUTHORITY ----------------------------------
// To calculate the LARGEST pure force or torque that a layout can make in one direction
// Calculated using the simplex solver
// 
// Linear Program unknowns for N thrusters :    x = [u_1,...,u_N, s, w_1,....,w_N]
// u_j      duty cycle of thruster j (0 = off and 1 = fully on)
// s        size of the output along d [N for a force, Nm for a torque]
// w_j      slack = unused part of the thruster j (the solver only takes equations)
// 
// Equations (authorityRows) :
// B u - s d = 0    on every axis k with rows[k] = true; s along d, nothing on other axis
//                  An axis set to false is left free (e.g. small side forces during entry roll do not matter)
// u_j + w_j = 1    healthy thruster, so u_j <= 1 (= 0 : failed thruster, held off)
//
// [A] [u_1; u_2; ... ; u_N ; s ; w_1; w_2; ...; w_N] = [b] = [0; 0;...; 0; 1; 1;...; 1]

inline void authorityRows(const Matrix& B, const Vec6& d, const std::array<bool, 6>& constrained,
                          const std::vector<bool>& active, Matrix& A, std::vector<double>& b) {
    
    const int num_thr = static_cast<int>(B[0].size());
    const int num_var = (2 * num_thr) + 1;
    for (int k = 0; k < 6; ++k) {
        if(!constrained[k])
            continue;
        std::vector<double> row(num_var, 0.0);
        for (int j = 0; j < num_thr; ++j)
            row[j] = B[k][j];
        row[num_thr] = -d[k];
        A.push_back(row);
        b.push_back(0.0);
    }
    for (int j =0; j < num_thr; ++j) {
        std::vector<double> row(num_var, 0.0);
        row[j] = 1.0;
        row[num_thr + 1 + j] = 1.0;
        A.push_back(row);
        b.push_back(active[j] ? 1.0 : 0.0);
    }
} 

// Stage 1 : Obtain the largest pure output s_max along d (0 if the layout cannot make it) :

inline double maxOutput(const Matrix& B, const Vec6& d, const std::array<bool, 6>& constrained, 
                        const std::vector<bool>& active) {
    const int num_thr = static_cast<int>(B[0].size());
    Matrix A;
    std::vector<double> b;
    authorityRows(B, d, constrained, active, A, b);
    std::vector<double> c(2 * num_thr + 1, 0.0);
    c[num_thr] = -1.0;                          // -ve because Simplex solver minimizes the cost
    SimplexResult result = solveSimplex(A, b, c);
    return (result.ok && result.x[num_thr] > 1e-9) ? result.x[num_thr] : 0.0;
}

// Stage 2 : Calculating efficiency of making the output s along d, all thrusters healthy :
// eff = s / (least total thrust fired) [1 = no thrust wasted]
// But the same output (s_max) that is passed from Stage 1 optimization can be made in several ways, 
// including some wasteful thrusters. 
// (E.g. extra jets firing to produce more force, but against each other to cancel the torque).
// The budget calculates the cheapest option. 
//
// Hence, Same equations plus "s_fixed" is passed to solver to minimize sum F_j u_j

inline double efficientOutput(const Matrix& B, const std::vector<double>& thrust_mag, const Vec6& d,
                              const std::array<bool, 6>& constrained, double s) {
    if (s <= 0)
        return 0.0;
    const int num_thr = static_cast<int>(B[0].size());
    Matrix A;
    std::vector<double> b;
    authorityRows(B, d, constrained, std::vector<bool>(num_thr, true), A, b);
    s *= 1.0 - 1e-7;                        // to avoid any round-off issue
    std::vector<double> fix(2 * num_thr + 1, 0.0);
    fix[num_thr] = 1.0;
    A.push_back(fix);
    b.push_back(s);
    std::vector<double> c(2 * num_thr +1, 0.0);
    for (int j = 0; j < num_thr; ++j)
        c[j] = thrust_mag[j];

    SimplexResult result = solveSimplex(A, b, c);
    if (!result.ok)
        throw std::runtime_error("Efficiency linear program did not solve");
    double used = 0;
    for (int j = 0; j < num_thr; ++j)
        used += thrust_mag[j] * result.x[j];
    return s / used;
}

// Counter of allocation LPs that did not solve (diagnostic : should stay 0).
inline int& allocFailures() { static int n = 0; return n; }

// ------------------------ THRUSTER ALLOCATION --------------------------------
// Used to assess which thruster to fire, and by how much, so that the vehicle gets the commanded
// force and torque (6 values) as closely as possible
//
// Linear Program Unknowns : x = [u_1,...,u_N, p_1,...,p_6, q_1,...,q_6, w_1,...,w_N]
//      u_j         duty cycle of thruster j (0...1)
//      p_k, q_k    shortfall / overshoot on axis k : Part of the commanded thrust that thrusters cannot make
//      w_j         slack, so that u_j + w_j = 1
//      
//      Equations           : b u + p - q = F
//                              u_j + w_j = 1
//
//      Cost                : sum u_j + (100 * torque_errors [N m]) + (forceWeight * force_errors [N]) 
//
// forceWeight = 1 for a burn as it matters, 0.01 during attitude control only

inline std::vector<double> allocate(const Matrix& B, const Vec6& force_torque, double forceWeight,
                                    const std::vector<bool>& active) {
    const int num_thr = static_cast<int>(B[0].size());
    const int num_var = num_thr + 12 + num_thr;
    Matrix A;
    std::vector<double> b;
    for (int k = 0; k < 6; ++k) {
        std::vector<double> row(num_var, 0.0);
        for (int j= 0; j < num_thr; ++j)
            row[j] = B[k][j];
        row[num_thr + k] = 1.0;
        row[num_thr + 6 + k] = -1.0;
        A.push_back(row);
        b.push_back(force_torque[k]);
    }

    for (int j = 0; j < num_thr; ++j) {
        std::vector<double> row(num_var, 0.0);
        row[j] = 1.0;
        row[num_thr + 12 + j] = 1.0;
        A.push_back(row);
        b.push_back(active[j] ? 1.0 : 0.0);
    }

    std::vector<double> c(num_var, 0.0);
    for (int j = 0; j < num_thr; ++j)
        c[j] = 1.0;
    for (int k = 0; k < 6; ++k)
        c[num_thr + k] = c[num_thr + 6 + k] = (k < 3) ? forceWeight : 100.0;     

    SimplexResult result = solveSimplex(A, b, c);

    std::vector<double> u(num_thr, 0.0);
    if (result.ok){
        for (int j = 0; j < num_thr; ++j)
            u[j] = clamp(result.x[j], 0.0, 1.0);
    }
    else
        ++allocFailures();
    
    return u;
}

inline std::vector<double> thrustVector(const Layout& Lay){
    std::vector<double> thrustVec;
    for (auto& t : Lay.thr)
        thrustVec.push_back(t.force_value);
    return thrustVec;
}


