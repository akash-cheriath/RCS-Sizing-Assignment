// Checks which thrust level F per thruster satisfies ALL sizing drivers of one architecture? 

#include <algorithm>
#include <cmath>
#include "analyses.hpp"
#include "orbit_params.hpp"

namespace {

// Angular acceleration needed to rotate theta in time t with the rate limit omega_max (rest to rest).
double accelForTime(double theta, double omega_max, double t) {
    if (t <= theta / omega_max) return 1e9;                          // impossible under the rate limit
    if (t >= 2.0 * theta / omega_max) return 4.0 * theta / (t * t);  // bang-bang, rate limit not reached
    return omega_max / (t - theta / omega_max);                           // accelerate - coast - brake
}

}  // namespace

ThrustWindow runThrustWindow(const Config& cfg, const Vehicle& veh, const LayoutResult& LayRes, const BudgetResult& Budget) {
    const double deg2Rad             = cnst::DEG;
    const double inertiaMarginNom    = 1.0 + cfg.numValue("mass_prop_margin");       // healthy cases
    const double inertiaMarginFail   = 1.0 + cfg.numValue("failure_case_margin");   // failure cases

    // ---- mission quantities
    const double rTarget    = cnst::RADIUS_EARTH + cfg.numValue("h_target") * 1e3;
    const double rPerigee   = orbit::combinedVehiclePerigee(rTarget, cnst::RADIUS_EARTH + cfg.numValue("perigee_target") * 1e3,
                                                    cfg.numValue("srvc_mod_sep_dv"));
    const double dvDeorbit  = orbit::deorbitBurn(rTarget, rPerigee) * (1.0 + cfg.numValue("dv_margin"));
    const double massD      = Budget.mass_at_deorbit;
    const double massProx   = Budget.mass_at_proxops;
    const double throttle   = cfg.numValue("main_burn_throttle");         // main burns use this share of the largest -X force
    const double tD0        = cfg.numValue("deorbit_burn_max") * 60.0;
    const double tD1        = cfg.numValue("deorbit_burn_max_1f") * 60.0;
    const double tD2        = cfg.numValue("deorbit_burn_max_2f") * 60.0;

    const Vec3 vReentry     = reentryVelocityDir(cfg.numValue("alpha_trim") * deg2Rad);
    const double Is         = dot(vReentry, veh.capsule.inertia * vReentry);                    // inertia about the velocity vector
    const double thRev      = cfg.numValue("bank_reversal") * deg2Rad;
    const double omegaMax   = cfg.numValue("roll_rate_max") * deg2Rad;
    const double aBank0     = accelForTime(thRev, omegaMax, cfg.numValue("bank_reversal_time"));
    const double aBank2     = accelForTime(thRev, omegaMax, cfg.numValue("bank_reversal_time_2f"));
    const double omega0     = cfg.numValue("tipoff_rate") * deg2Rad;
    const double tDetumble  = cfg.numValue("detumble_time");
    const double camDv      = cfg.numValue("cam_dv"); 
    const double camT       = cfg.numValue("cam_time");
    const double tolOmega   = cfg.numValue("dock_ang_rate_max") * deg2Rad;
    const double frac       = cfg.numValue("mib_tol");

    ThrustWindow Window;
    auto perN = [](double a) { 
        return a > 1e-9 ? a : 1e-9; 
    };                                       
                                                       // lost authority -> huge thrust
    auto add = [&](const std::string& d, const std::string& fc, const std::string& basis, double F, bool isMin,
                   double Frel = -1.0) { 
        Window.drivers.push_back({d, fc, basis, F, isMin, Frel < 0 ? F : Frel}); 
    };

    // 1-3. Deorbit burn duration (combined vehicle, pure -X force)
    add("deorbit burn <= " + cfg.strValue("deorbit_burn_max", "") + " min", "none", "m dv / t / a(-Fx)",
        massD * dvDeorbit / tD0 / perN(throttle * LayRes.nominal[0][1]), true);
    add("deorbit burn <= " + cfg.strValue("deorbit_burn_max_1f", "") + " min", "1 fail", "m dv / t / a1(-Fx)",
        massD * dvDeorbit / tD1 / perN(throttle * LayRes.fail1[0][1]), true);
    add("deorbit burn <= " + cfg.strValue("deorbit_burn_max_2f", "") + " min", "2 fail", "m dv / t / a2(-Fx)",
        massD * dvDeorbit / tD2 / perN(throttle * LayRes.fail2[0][1]), true, massD * dvDeorbit / tD2 / perN(throttle * LayRes.axialBest[2]));

    // 4-5. Reentry bank reversal (capsule, torque about the velocity vector)
    add("bank reversal in " + cfg.strValue("bank_reversal_time", "") + " s", "none", "I_vel (1+m) alpha / a_bank",
        Is * inertiaMarginNom * aBank0 / perN(LayRes.bank[0]), true);
    add("bank reversal in " + cfg.strValue("bank_reversal_time_2f", "") + " s", "2 fail", "I_vel (1+m_f) alpha / a2_bank",
        Is * inertiaMarginFail * aBank2 / perN(LayRes.bank[2]), true);


    // 6. Collision avoidance manoeuvre (nose at the station: -X), 2 failures
    add("CAM " + cfg.strValue("cam_dv", "") + " m/s in " + cfg.strValue("cam_time", "") + " s", "2 fail",
        "m dv / t / a2(-Fx)", massProx * camDv / camT / perN(LayRes.fail2[0][1]), true,
        massProx * camDv / camT / perN(LayRes.axialBest[2]));

    // 7. Detumble (combined vehicle), worst axis, healthy and 1-failure cases
    {
        double Fd = 0;
        for (int ax = 0; ax < 3; ++ax) {
            const double I = veh.combined_vehicle.inertia.matrix[ax][ax];
            const double a0 = std::min(LayRes.torqueFree[0][2 * ax], LayRes.torqueFree[0][2 * ax + 1]);
            const double a1 = std::min(LayRes.torqueFree1[0][2 * ax], LayRes.torqueFree1[0][2 * ax + 1]);
            Fd = std::max({Fd, I * inertiaMarginNom * omega0 / tDetumble / perN(a0), I * inertiaMarginFail * omega0 / tDetumble / perN(a1)});
        }
        add("detumble " + cfg.strValue("tipoff_rate", "") + " deg/s in " + cfg.strValue("detumble_time", "") + " s",
            "none / 1 fail", "I (1+m) w0 / t / a", Fd, true);
    }

    // 8. Every pure axis after one failure (prox-ops must stay controllable)
    if (!LayRes.oneFT) 
        add("6-DOF pure control", "1 fail", "an axis is lost", 1e9, true);

    // 9. Upper bound: one minimum pulse of one thruster (docking fine control, rotation). IDSS Table 3.3.1.1-2
    //    limits the vector sum of the pitch/yaw rates and, separately, the roll rate (both to tolOmega)
    
    const std::string pct = std::to_string(static_cast<int>(std::round(100 * frac)));
    {
        const double FmaxRoll     = frac * tolOmega * veh.combined_vehicle.inertia.matrix[0][0] /
                                    (veh.thr.t_min * std::max(LayRes.maxArm[0], 1e-9));
        const double FmaxPitchYaw = frac * tolOmega / (veh.thr.t_min * std::max(LayRes.maxRatePitchYaw, 1e-12));
        add("single-jet rotation pulse <= " + pct + " % of IDSS rate tol (pitch/yaw vector sum)", "-",
            "mib_tol tol / (tmin dw_py); roll: mib_tol tol Ixx / (tmin arm_x)", std::min(FmaxRoll, FmaxPitchYaw), false);
    }

    // 10. Upper bound: one minimum pulse of the largest pure sideways force must not change the lateral
    //     velocity by more than mib_tol of the IDSS lateral-rate tolerance (fine control inside the
    //     approach corridor; conservative: all sideways jets fire for one minimum pulse)
    {
        const double aLat = std::max({LayRes.nominal[0][2], LayRes.nominal[0][3], LayRes.nominal[0][4], LayRes.nominal[0][5]});
        const double tolVel = cfg.numValue("dock_lat_rate_max");
        add("single lateral pulse <= " + pct + " % of IDSS lateral-rate tol", "-", "mib_tol tol_v m / (tmin a_lat)",
            frac * tolVel * massProx / (veh.thr.t_min * std::max(aLat, 1e-9)), false);
    }

    Window.F_max = 1e300;
    for (const auto& d : Window.drivers) {
        if (d.isMin && d.thrust > Window.F_min) { 
            Window.F_min = d.thrust; 
            Window.driver_min = d.name + " (" + d.failure_case + ")"; 
        }
        if (d.isMin && d.thrust_relax > Window.F_min_relaxed) { 
            Window.F_min_relaxed = d.thrust_relax; 
            Window.driver_min_relaxed = d.name + " (" + d.failure_case + ")"; 
        }
        if (!d.isMin && d.thrust < Window.F_max) { 
            Window.F_max = d.thrust; 
            Window.driver_max = d.name; }
    }
    return Window;
}