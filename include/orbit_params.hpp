/*
Defines the orbit class, which is used to represent the orbital parameters of a spacecraft. This helps to size
every orbital maneuvers with the required delta-v, that is required for sizing the RCS thrusters. 
*/

#pragma once
#include "math_structs.hpp"
#include "const_def.hpp"
#include "config_map.hpp"

namespace orbit {

    inline double orbitalVelocity(double r)     { return std::sqrt(cnst::MU_EARTH / r); }
    inline double meanMotion(double r)          { return std::sqrt(cnst::MU_EARTH / (r * r * r)); }

    // Two-impulse maneuver delta-v calculation using the Hohmann transfer approximation between two circular orbits 
    // r1: initial orbit radius, r2: final orbit radius
    struct HohmannTransfer { double deltaV1, deltaV2; };
    inline HohmannTransfer hohmannTransfer(double r1, double r2) {
        double v1 = orbitalVelocity(r1);
        double v2 = orbitalVelocity(r2);
        double a_transfer = (r1 + r2) / 2.0;
        double v_transfer1 = std::sqrt(cnst::MU_EARTH * (2.0 / r1 - 1.0 / a_transfer));
        double v_transfer2 = std::sqrt(cnst::MU_EARTH * (2.0 / r2 - 1.0 / a_transfer));
        double deltaV1 = std::fabs(v_transfer1 - v1);
        double deltaV2 = std::fabs(v2 - v_transfer2);
        return {deltaV1, deltaV2};
    }

    // Single retrograde burn from a circular orbit to reduce the perigee to r_perigee during a deorbit maneuver.
    inline double deorbitBurn(double r_circular, double r_perigee) {
        double v_circular = orbitalVelocity(r_circular);
        double a_transfer = (r_circular + r_perigee) / 2.0;
        double v_transfer = std::sqrt(cnst::MU_EARTH * (2.0 / r_circular - 1.0 / a_transfer));
        return std::fabs(v_circular - v_transfer);
    }

    // Perigee the COMBINED vehicle (capsule + SM) must aim at in the deorbit burn.
    // After SM separation the capsule fires a small posigrade burn dv at apogee, which raises its perigee to
    // r_perigee_capsule. So the combined vehicle aims a little lower: its apogee speed is the capsule's minus dv.
    
    inline double combinedVehiclePerigee(double r_apogee, double r_perigee_capsule, double dv) {
        double v_apo_capsule = std::sqrt(2.0 * cnst::MU_EARTH * r_perigee_capsule / (r_apogee * (r_apogee + r_perigee_capsule)));
        double v_apo_combined = v_apo_capsule - dv;
        return v_apo_combined * v_apo_combined * r_apogee * r_apogee / (2.0 * cnst::MU_EARTH - v_apo_combined * v_apo_combined * r_apogee);
    }

    // Simple exponential atmospheric density model for low Earth orbit (LEO) altitudes. This can be used to estimate drag effects on 
    // the spacecraft. 
    // rho = rho0 * exp(-h / H), where rho0 is the reference density at sea level, h is the altitude above Earth's surface, and 
    //                                 H is the scale height.
    // Source : Vallado, "Fundamentals of Astrodynamics and Applications", 4th Edition, 2013, Table 8.4 "Exponential Atmosphere Model"

    inline double densityAtAltitude(double altitude) {
        struct AtmosphereLayer {
            double h0;      // Base altitude [m]
            double rho0;    // Base density [kg/m^3]
            double H;       // Scale height [m]
        };

        static const AtmosphereLayer layers[] = {
            {0.0, 1.225, 7.249e3},         {25e3, 3.899e-2, 6.349e3},    {30e3, 1.774e-2, 6.682e3},
            {40e3, 3.972e-3, 7.554e3},     {50e3, 1.057e-3, 8.382e3},    {60e3, 3.206e-4, 7.714e3},
            {70e3, 8.770e-5, 6.549e3},     {80e3, 1.905e-5, 5.799e3},    {90e3, 3.396e-6, 5.382e3},
            {100e3, 5.297e-7, 5.877e3},    {110e3, 9.661e-8, 7.263e3},   {120e3, 2.438e-8, 9.473e3},
            {130e3, 8.484e-9, 12.636e3},   {140e3, 3.845e-9, 16.149e3},  {150e3, 2.070e-9, 22.523e3},   
            {180e3, 5.464e-10, 29.740e3},  {200e3, 2.789e-10, 37.105e3}, {250e3, 7.248e-11, 45.546e3},
            {300e3, 2.418e-11, 53.628e3},  {350e3, 9.518e-12, 53.298e3}, {400e3, 3.725e-12, 58.515e3},  
            {450e3, 1.585e-12, 60.828e3},  {500e3, 6.967e-13, 63.822e3}
        };

        const int numLayers = sizeof(layers) / sizeof(layers[0]);
        int i = 0;
        while (i < numLayers - 1 && altitude >= layers[i + 1].h0) {
            ++i;
        }
        return layers[i].rho0 * std::exp((layers[i].h0 - altitude) / layers[i].H);
    }

    // Quantifying the finite-burn loss of a long tangential burn as per ESA guideline (eta_gravity_loss), as the dv budget 
    // is impulsive. The deorbit burn may last ~ 13 min, and the capusle may cover ~50 deg around orbit during this burn.
    // Real burn needs a higher dv than the ideal dv given by the thrust applied at an ideal point.
    inline double finiteBurnEfficiency (double arc){
        if (arc < 1e-6) return 1.0;
        return 0.5 * (1.0 + std::sin(0.5 * arc) / (0.5 * arc));
    }


    // Clohessy-Wiltshire (CW) model:
    // x'' = 2 n z'             + ax
    // y'' = -n^2 y             + ay
    // z'' = -2 n x' + 3 n^2 z  + az,
    // defined in LVLH frame where x-axis = +V-bar (along target velocity), y-axis = +H-bar (opp. to orbit ang momentum)
    // and z-axis = +R-bar (towards centre of Earth)

    // State Transition sub-matrices of above equation sourced from "Terminal Guidance System for Satellite Rendezvous"
    // [r(t) ; v(t)] = [phi_rr(t) phi_rv(t); phi_vr(t) phi_vv(t)] [r0 ; v0], where r and v are 3*1 vectors and phi_xx are
    // each 3x3 matrices

    struct CWStm { Mat3 phi_rr, phi_rv, phi_vr, phi_vv; };
    inline CWStm cwStm(double n, double t) {
        double s = std::sin(n * t);
        double c = std::cos(n * t);

        CWStm phi;
        
        // The below assignments are done for non-zero entries of the matrices. Rest of the entries are 0 assigned automatically
        phi.phi_rr.matrix[0][0] = 1.0;
        phi.phi_rr.matrix[0][2] = 6.0 * ((n * t) - s);
        phi.phi_rr.matrix[1][1] = c;
        phi.phi_rr.matrix[2][2] = 4.0 - (3.0 * c);

        phi.phi_rv.matrix[0][0] = (1.0 / n) * ((4.0 * s) - (3.0 * n * t));
        phi.phi_rv.matrix[0][2] = (2.0 / n) * ( 1.0 - c );
        phi.phi_rv.matrix[1][1] = (1.0 / n) * s;
        phi.phi_rv.matrix[2][0] = (2.0 / n) * ( c - 1.0 );
        phi.phi_rv.matrix[2][2] = (1.0 / n) * s;

        phi.phi_vr.matrix[0][2] = 6.0 * n * (1.0 - c);
        phi.phi_vr.matrix[1][1] = -n * s;
        phi.phi_vr.matrix[2][2] = 3.0 * n * s;

        phi.phi_vv.matrix[0][0] = (4.0 * c) - 3.0;
        phi.phi_vv.matrix[0][2] = 2.0 * s;
        phi.phi_vv.matrix[1][1] = c;
        phi.phi_vv.matrix[2][0] = - 2.0 * s;
        phi.phi_vv.matrix[2][2] = c;

        return phi;
    }

    // CW two-impulse transfer : from (r0, v0) reach rf in time T, then stop : Eg WP0 -> WP1
    // Position after T : rf = phi_rr r0 + phi_rv v0+
    // First burn dv1   : dv1 = v0+ - v0
    // Arrival velocity vf : vf = phi_vr r0 + phi_vv v0+
    // Second burn dv2  : dv2 = -vf (To stop the vehicle)
    // phi_rv singular for a transfer of whole number of orbits - T should avoid it!!
    inline void cwTwoImpulse(double n, const Vec3& r0, const Vec3& v0, const Vec3& rf, double T, Vec3& dv1, Vec3& vf){
        CWStm phi = cwStm(n, T);
        Vec3 v0_plus = phi.phi_rv.inverse() * (rf - (phi.phi_rr * r0));
        dv1 = v0_plus - v0;
        vf = (phi.phi_vr * r0) + (phi.phi_vv * v0_plus);
    }

    // Natural (no-thrust) CW accelerations
    inline Vec3 cwAccel(double n, const Vec3& r, const Vec3& v) {
        auto cw_accel_x = 2.0 * n * v.z;
        auto cw_accel_y = -n * n * r.y;
        auto cw_accel_z = (-2.0 * n * v.x) + (3.0 * n * n * r.z);
        return {cw_accel_x, cw_accel_y, cw_accel_z};
    }

} // namespace orbit