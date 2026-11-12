// Translates the config files into C++ structures. 

#pragma once
#include "thruster.hpp"
#include "config_map.hpp"
#include "const_def.hpp"
#include "math_structs.hpp"

struct ThrusterPerf {
    double F;           // [N] Nominal thrust
    double isp;         // [s] Steady-state Isp
    double isp_pulse;   // [s] Effective Isp in short-pulse operation
    double t_min;       // [s] Minimum on-time
};

struct Vehicle {
    MassProps combined_vehicle;     // Combined capsule + service module
    MassProps capsule;
    PodGeometry podGeom;
    ThrusterPerf thr;
};

inline MassProps loadMassProps(const Config& conf, const std::string& name) {
    MassProps mp;
    mp.name = name;
    mp.mass = conf.numValue(name + "_mass");
    mp.cog = Vec3(conf.numValue(name + "_cog_x"), conf.numValue(name + "_cog_y"), conf.numValue(name + "_cog_z"));
    mp.inertia = Mat3::diagonal(conf.numValue(name + "_Ixx"), conf.numValue(name + "_Iyy"), conf.numValue(name + "_Izz"));
    return mp;
}

inline Vehicle loadVehicle(const Config& conf) {
    Vehicle veh;
    veh.combined_vehicle = loadMassProps(conf, "combined");
    veh.capsule = loadMassProps(conf, "capsule");
    veh.podGeom.base_radius = 0.5 * conf.numValue("capsule_base_diameter");
    veh.podGeom.cant_aft = conf.numValue("cant_aft") * cnst::DEG;
    veh.podGeom.cant_fwd = conf.numValue("cant_fwd") * cnst::DEG;
    veh.podGeom.cant_roll = conf.numValue("cant_roll") * cnst::DEG;
    veh.podGeom.clock0 = conf.numValue("pod_clock0") * cnst::DEG;
    veh.podGeom.force = conf.numValue("thrust");
    veh.podGeom.half_angle = conf.numValue("sidewall_half_angle") * cnst::DEG;
    veh.podGeom.nose_r = conf.numValue("nose_thr_radius");
    veh.podGeom.nose_x = conf.numValue("nose_thr_x");
    veh.podGeom.pod_x = conf.numValue("pod_x_pos");
    veh.thr.F = veh.podGeom.force;
    veh.thr.isp = conf.numValue("isp_steady");
    veh.thr.isp_pulse = conf.numValue("isp_pulse");
    veh.thr.t_min = conf.numValue("min_pulse_time") * 1e-3;
    return veh;
}
