// Calculates the delta-v and propellant budget

#include <algorithm>
#include <cmath>
#include <iostream>
#include "analyses.hpp"
#include "csv_export.hpp"
#include "orbit_params.hpp"

BudgetResult runDvBudget(const Config& cfg, const Vehicle& veh, const LayoutResult& LayRes, const std::string& outDir) {
    const bool verbose      = !outDir.empty();
    const double g0         = cnst::G0; 
    const double deg2Rad    = cnst::DEG;
    const double isp        = veh.thr.isp;
    const double ispPulse   = veh.thr.isp_pulse;
    const double r1         = cnst::RADIUS_EARTH + cfg.numValue("h_insertion") * 1e3;
    const double r2         = cnst::RADIUS_EARTH + cfg.numValue("h_target") * 1e3;
    const double inc        = cfg.numValue("inclination") * deg2Rad;
    const double n2         = orbit::meanMotion(r2);
    const double margin     = 1.0 + cfg.numValue("dv_margin");

    const double etaMain    = LayRes.effBurn;                   // main burns on the FWD (-X) thrusters
    const double etaRcs     = LayRes.effAvg[0];
    const double etaRcsCaps = LayRes.effAvg[1];
    const double arm        = (LayRes.maxArm[0] + LayRes.maxArm[1] + LayRes.maxArm[2]) / 3.0;   // torque arm per newton

    double mass = veh.combined_vehicle.mass;
    int cfgIdx = 0;                                     // 0 combined vehicle, 1 capsule
    double propTotal = 0, dvTotal = 0, attTotal = 0;
    BudgetResult res;

    Csv out(outPath(outDir, "dv_budget.csv"));
    out.header({"phase", "category", "config", "dv_m_s", "isp_s", "efficiency", "mass_start_kg", "propellant_kg",
                "mass_end_kg", "note"});
    auto trans = [&](const std::string& name, const std::string& cat, double dv, double ispSteady, double eta,
                     const std::string& note) {
        const double m0 = mass; 
        const double mp = m0 * (1.0 - std::exp(-dv / (eta * ispSteady * g0)));
        mass -= mp; 
        propTotal += mp; 
        dvTotal += dv;
        out << name << cat << (cfgIdx ? "capsule" : "combined") << dv << ispSteady << eta << m0 << mp << mass << note;
        out.endl();
    };

    auto attitude = [&](const std::string& name, double mp, const std::string& note) {
        const double m0 = mass;
        mass -= mp; 
        propTotal += mp; 
        attTotal += mp;
        out << name << "attitude" << (cfgIdx ? "capsule" : "combined") << 0.0 << ispPulse << 1.0 << m0 << mp << mass << note;
        out.endl();
    };

    // --------------------------------- LAUNCH AND RENDEZVOUS ----------------------------------------
    {   
        // Detumble: angular momentum h = |I w0| with the tip-off rate on all axes
        const double omega0 = cfg.numValue("tipoff_rate") * deg2Rad;
        const Vec3 h = veh.combined_vehicle.inertia * Vec3(omega0, omega0, omega0);
        attitude("Detumble after launcher separation", cfg.numValue("detumble_pulse_factor") * norm(h) / (arm * ispPulse * g0),
                 "factor x |I w0| / (arm Isp g0)");
    }

    {   
        // Launcher dispersions (3 sigma): altitude = Hohmann from the low orbit instead of the nominal one;
        // Plane (inclination + RAAN) combined with the circularisation burn at 400 km (vector addition)
        const auto hNom = orbit::hohmannTransfer(r1, r2);
        const auto hLow = orbit::hohmannTransfer(r1 - cfg.numValue("dispersion_h") * 1e3, r2);
        const double dvAlt = (hLow.deltaV1 + hLow.deltaV2) - (hNom.deltaV1 + hNom.deltaV2);
        const double dInc = cfg.numValue("dispersion_inc") * deg2Rad; 
        const double dRaan = cfg.numValue("dispersion_raan") * deg2Rad;
        const double dvPlane = orbit::orbitalVelocity(r2) * std::sqrt(dInc * dInc + std::pow(std::sin(inc) * dRaan, 2));
        const double dvComb = std::sqrt(hNom.deltaV2 * hNom.deltaV2 + dvPlane * dvPlane) - hNom.deltaV2;
        trans("Injection dispersion correction (3 sigma)", "rendezvous", dvAlt + dvComb, isp, etaMain,
              "altitude + plane (combined with circularisation)");
    }

    {   
        // Phasing at the insertion altitude: duration from the phase angle; drag make-up (high solar activity)
        const double tPh = cfg.numValue("phase_angle") * deg2Rad / (orbit::meanMotion(r1) - n2);
        const double v = orbit::orbitalVelocity(r1);
        const double aDrag = 0.5 * cfg.numValue("density_solar_factor") * orbit::densityAtAltitude(r1 - cnst::RADIUS_EARTH) *
                          v * v * cfg.numValue("drag_cd") * cfg.numValue("drag_area") / mass;
        trans("Drag make-up during phasing", "rendezvous", aDrag * tPh, isp, etaMain,
              "phasing " + std::to_string(tPh / 3600.0).substr(0, 4) + " h");
    }

    {   
        const auto h = orbit::hohmannTransfer(r1, r2);
        trans("Orbit raise 200 -> 400 km", "rendezvous", (h.deltaV1 + h.deltaV2) * margin, isp, etaMain,
              "Hohmann x (1 + dv margin)");
    }

    trans("Terminal phase (TPI + midcourse + V-null)", "rendezvous", cfg.numValue("tpi_mcc_v_null_dv"), isp, etaMain, "Orion");
    trans("Rendezvous dispersions (3 sigma)", "rendezvous", cfg.numValue("dispersion_terminal_dv"), isp, etaMain, "Orion");

    // --------------------------------- PROXIMITY OPERATION ----------------------------------------

    res.mass_at_proxops = mass;
    double dvProx = 0;

    {   
        // WP0 hold on R-bar: hovering below the station needs a = 3 n^2 z
        const double z = cfg.numValue("wp0_below_trgt_rbar"); 
        const double T = cfg.numValue("wp0_hold") * 60.0;
        const double dv = 3.0 * n2 * n2 * z * T;
        dvProx += dv;
        trans("WP0 hold below the station (R-bar)", "prox-ops", dv, ispPulse, etaRcs, "3 n^2 z T");
    }

    {   
        // WP0 -> WP1: two-impulse Clohessy-Wiltshire transfer
        Vec3 dv1, vf;
        orbit::cwTwoImpulse(n2, Vec3(0, 0, cfg.numValue("wp0_below_trgt_rbar")), Vec3(),
                            Vec3(cfg.numValue("wp1_vbar"), 0, 0), cfg.numValue("wp0_1_orbit_frac") * 2.0 * cnst::PI / n2,
                            dv1, vf);
        const double dv = norm(dv1) + norm(vf);
        dvProx += dv;
        trans("WP0 -> WP1 transfer", "prox-ops", dv, ispPulse, etaRcs, "CW two-impulse |dv1| + |dv2|");
    }

    dvProx += cfg.numValue("wp1_hold_dv");
    trans("WP1 hold on +V-bar", "prox-ops", cfg.numValue("wp1_hold_dv"), ispPulse, etaRcs, "allowance");
    
    {   
        // V-bar approach WP1 -> WP2: Coriolis compensation 2 n v t + start / stop impulses
        const double x1 = cfg.numValue("wp1_vbar"); 
        const double x2 = cfg.numValue("wp2_vbar");
        const double vf = cfg.numValue("approach_speed_far");
        const double vn = cfg.numValue("approach_speed_near");
        const double xs = std::max(std::min(cfg.numValue("approach_switch_range"), x1), x2);
        const double t1 = (x1 - xs) / vf; 
        const double t2 = (xs - x2) / vn;
        const double dv = 2.0 * n2 * (vf * t1 + vn * t2) + vf + std::fabs(vf - vn) + vn;
        dvProx += dv;
        trans("V-bar approach WP1 -> WP2", "prox-ops", dv, ispPulse, etaRcs, "2 n v t + impulses");
    }

    {   
        // WP2 hold, final approach at the contact speed, push after contact
        const double vc = cfg.numValue("contact_speed");
        const double x2 = cfg.numValue("wp2_vbar");
        const double dv = cfg.numValue("wp2_hold_dv") + 2.0 * n2 * x2 + vc + cfg.numValue("capture_push_dv");
        dvProx += dv;
        trans("WP2 hold + final approach + capture push", "prox-ops", dv, ispPulse, etaRcs, "allowance + 2 n x + impulse");
    }

    trans("Prox-ops control allowance", "prox-ops", cfg.numValue("proxops_allowance") * dvProx, ispPulse, etaRcs,
          "pulsing, navigation noise, cross-coupling");

    // -------------------------------------- CONTINGENCIES ---------------------------------------------

    trans("Collision avoidance manoeuvres", "contingency", cfg.numValue("cam_count") * cfg.numValue("cam_dv"), isp,
          etaMain, "CAM");
    trans("Retreat + second docking attempt", "contingency", cfg.numValue("reapproach_dv"), ispPulse, etaRcs, "CCT-REQ-1130");
    trans("Re-rendezvous after 24 h", "contingency", cfg.numValue("rerendezvous_dv"), isp, etaMain, "CCT-REQ-1130");
    trans("Debris avoidance manoeuvres", "contingency",
          cfg.numValue("debris_avoid_mnvr_count") * cfg.numValue("debris_avoid_mnvr_dv"), isp, etaMain, "free flight");
    
    // ------------------------------- UNDOCKING, DEPARTURE, DEORBIT -------------------------------------

    trans("Undocking push", "departure", cfg.numValue("undock_push_dv"), ispPulse, etaRcs, "");
    trans("Departure burns", "departure", cfg.numValue("departure_dv"), ispPulse, etaRcs, "leave the approach ellipsoid");
    trans("Phasing for the landing opportunity", "departure", cfg.numValue("deorbit_phasing_dv"), isp, etaMain, "");

    res.mass_at_deorbit = mass;

    {   
        // Deorbit burn: the combined vehicle is aimed below the capsule entry perigee (the capsule raise burn
        // after SM separation lifts it onto the target); long burn -> finite-burn loss (ESA)
        const double rp         = orbit::combinedVehiclePerigee(r2, cnst::RADIUS_EARTH + cfg.numValue("perigee_target") * 1e3,
                                                        cfg.numValue("srvc_mod_sep_dv"));
        const double dvImp      = orbit::deorbitBurn(r2, rp);
        const double forceMain  = cfg.numValue("main_burn_throttle") * LayRes.nominal[0][1] * veh.thr.F;   // main-burn force [N]
        const double tBurn      = mass * dvImp / forceMain;
        const double eta        = orbit::finiteBurnEfficiency(tBurn * n2);
        res.deorbit_burn_min    = tBurn / 60.0;

        trans("Deorbit burn", "deorbit", dvImp * margin / eta, isp, etaMain,
              "impulsive " + std::to_string(dvImp).substr(0, 5) + " m/s, burn " + std::to_string(tBurn / 60.0).substr(0, 4) +
              " min, finite-burn eta " + std::to_string(eta).substr(0, 5));
    }

    {   
        const double m0 = mass;
        mass -= veh.combined_vehicle.mass - veh.capsule.mass;
        out << "Service module separation" << "mass event" << "combined" << 0.0 << 0.0 << 0.0 << m0 << 0.0 << mass
            << "SM burns up";
        out.endl();
    }

    cfgIdx = 1;
    trans("Capsule raise burn", "deorbit", cfg.numValue("srvc_mod_sep_dv"), ispPulse, etaRcsCaps, "lifts the capsule perigee");

    // -------------------------------------- REENTRY ---------------------------------------------------

    {
        const Vec3 vReentry = reentryVelocityDir(cfg.numValue("alpha_trim") * deg2Rad);
        const double Is = dot(vReentry, veh.capsule.inertia * vReentry);                          // inertia about the velocity vector
        const double omega = cfg.numValue("roll_rate_max") * deg2Rad;
        const double omega_s = cfg.numValue("slew_rate") * deg2Rad;
        attitude("Pre-entry 180 deg pitch", 2.0 * veh.capsule.inertia.matrix[1][1] * omega_s / (arm * ispPulse * g0),
                 "2 I w / (arm Isp g0)");
        const double mRev = cfg.numValue("bank_reversal_count") * 2.0 * Is * omega / (arm * ispPulse * g0);
        attitude("Bank reversals", mRev, "count x 2 I w / (arm Isp g0)");
        attitude("Entry rate damping", std::max(mRev, cfg.numValue("entry_rate_damping")), "allowance");
    }

    // ----------------------------------- ATTITUDE (Whole Mission) -----------------------------------------

    cfgIdx = 0;
    {   
        const Mat3& I3 = veh.combined_vehicle.inertia;
        const double I = (I3.matrix[0][0] + I3.matrix[1][1] + I3.matrix[2][2]) / 3.0;           // mean inertia

        // Limit cycle with couples (Orion orbit-control paper): (2 F tmin arm)^2 / (4 I deadband arm Isp g0) per axis
        const double h = 2.0 * veh.thr.F * veh.thr.t_min * arm;
        const double mdot = 3.0 * h * h / (4.0 * I * cfg.numValue("deadband") * deg2Rad * arm * ispPulse * g0);

        attitude("Attitude hold limit cycle (free flight)", mdot * cfg.numValue("free_flight_time") * 3600.0,
                 "3 axes, deadband");
        attitude("Attitude slews", cfg.numValue("slew_count") * 2.0 * I * cfg.numValue("slew_rate") * deg2Rad / (arm * ispPulse * g0),
                 "count x 2 I w / (arm Isp g0)");
    }
    attitude("Attitude-control propellant margin (ESA)", cfg.numValue("attitude_prop") * attTotal,
             std::to_string(static_cast<int>(std::round(100 * cfg.numValue("attitude_prop")))) + " % of the attitude lines");

    
    // --------------------------------------- TOTAL DELTA V ------------------------------------------------

    res.dv_total = dvTotal;
    res.propellant_total = propTotal;
    res.propellant_loaded = propTotal / (1.0 - cfg.numValue("prop_residual"));
    out << "TOTAL usable propellant" << "total" << "-" << dvTotal << 0 << 0 << veh.combined_vehicle.mass << propTotal << 0 << "";
    out.endl();
    out << "Loaded propellant (incl. residuals)" << "total" << "-" << 0 << 0 << 0 << 0 << res.propellant_loaded << 0
        << "/(1 - residual)";
    out.endl();
    if (verbose)
        std::cout << "  " << LayRes.name << " with " << veh.thr.F << " N, Isp " << veh.thr.isp << " s: delta-v " << dvTotal
                  << " m/s, usable " << propTotal << " kg, loaded " << res.propellant_loaded << " kg ("
                  << 100.0 * res.propellant_loaded / veh.combined_vehicle.mass << " % of the vehicle), deorbit burn "
                  << res.deorbit_burn_min << " min\n";
    return res;
}


