// Checks the safe operation in the proximity of the target vehicle

#include <cmath>
#include <iostream>
#include <sstream>
#include "analyses.hpp"
#include "csv_export.hpp"
#include "orbit_params.hpp"

namespace {

struct State { Vec3 r, v; };

// What happened along one free-flight track with impulsive +V-bar burns
struct Track {
    double minRange     = 1e9; 
    double tMinRange    = 0;   // closest approach to the station
    double maxCone      = 0;                     // largest angle off +V-bar while still inside the KOS [rad]
    double tExit        = -1;                      // first time outside the KOS [s] (-1 : never)
    double tEnter       = -1;                     // first time inside the KOS after having been outside [s]
    double overflight   = -1;                 // range when the vehicle first passes over / behind the station (x <= 0)
    double dvTotal      = 0;
};

struct Burn { double t, dv; };              // time [s], delta-v along +x (V-bar) [m/s]

Track fly(double n, State state, const std::vector<Burn>& burns, double tEnd, double kos, Csv& log,
          const std::string& name, bool startInsideKos) {
    const double dt = 1.0;
    const auto P = orbit::cwStm(n, dt);
    Track track;
    size_t ib = 0;
    bool outside = !startInsideKos;
    for (double t = 0; t <= tEnd + 1e-9; t += dt) {
        while (ib < burns.size() && burns[ib].t <= t + 1e-9) { state.v.x += burns[ib].dv; track.dvTotal += std::fabs(burns[ib].dv); ++ib; }
        const double rg = norm(state.r);
        if (t > 0 && rg < track.minRange) { track.minRange = rg; track.tMinRange = t; }
        if (!outside) {                                         // still inside the KOS after undocking
            if (rg > 10.0) track.maxCone = std::max(track.maxCone, std::acos(clamp(state.r.x / rg, -1.0, 1.0)));
            if (rg > kos) { outside = true; track.tExit = t; }
        } else if (rg < kos && track.tEnter < 0) {
            track.tEnter = t;                                      // (re-)entered the KOS
        }
        if (track.overflight < 0 && t > 0 && state.r.x <= 0.0 && startInsideKos) track.overflight = rg;
        if (static_cast<long>(t) % 10 == 0)
            log << name << t << state.r.x << state.r.y << state.r.z << rg, log.endl();

        State nx{P.phi_rr * state.r + P.phi_rv * state.v, P.phi_vr * state.r + P.phi_vv * state.v};
        state = nx;
    }
    return track;
}

std::string fmt(double v, int digits = 1) {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(digits);
    o << v;
    return o.str();
}

}  // namespace

void runProximityChecks(const Config& cfg, const std::string& outDir) {
    const double r2         = cnst::RADIUS_EARTH + cfg.numValue("h_target") * 1e3;
    const double n          = orbit::meanMotion(r2); 
    const double period     = 2.0 * cnst::PI / n;
    const double kos        = cfg.numValue("kos_radius"); 
    const double margin     = cfg.numValue("kos_margin_min");
    const double tSafe      = cfg.numValue("passive_safety_time") * 3600.0;
    const double corridor   = cfg.numValue("app_corridor_half_angle") * cnst::DEG;
    const double cone       = cfg.numValue("dep_corridor_angle") * cnst::DEG;

    Csv log(outPath(outDir, "proximity_trajectories.csv"));
    log.header({"case", "t_s", "x_vbar_m", "y_hbar_m", "z_rbar_m", "range_m"});
    Csv out(outPath(outDir, "proximity_checks.csv"));
    out.header({"id", "requirement", "value", "limit", "status", "note"});
    auto row = [&](const std::string& id, const std::string& req, const std::string& val, const std::string& lim,
                   const std::string& status, const std::string& note) {
        out << id << req << val << lim << status << note;
        out.endl();
        std::cout << "  " << id << "  " << status << "  " << req << ": " << val << " (" << lim << ")"
                  << (note.empty() ? "" : " - " + note) << "\n";
    };

    // --------------- Approach: WP0 (R-bar, below) -> WP1 (+V-bar), two-impulse CW transfer

    const Vec3 rWp0(0, 0, cfg.numValue("wp0_below_trgt_rbar")), rWp1(cfg.numValue("wp1_vbar"), 0, 0);
    const double T = cfg.numValue("wp0_1_orbit_frac") * period;
    Vec3 dv1, vf;
    orbit::cwTwoImpulse(n, rWp0, Vec3(), rWp1, T, dv1, vf);
    const Track tf = fly(n, {rWp0, dv1}, {}, T, kos, log, "WP0 -> WP1 transfer", false);

    {
        const std::string st = tf.minRange < kos ? "FAIL" : tf.minRange < kos * (1.0 + margin) ? "FLAG" : "PASS";
        row("MIS-04a", "WP0 -> WP1 transfer outside the KOS", "closest " + fmt(tf.minRange, 0) + " m at " +
            fmt(tf.tMinRange / 60.0) + " min", "> " + fmt(kos, 0) + " m, margin " + fmt(100 * margin, 0) + " %", st,
            tf.tMinRange > T - 2.0 ? "closest point is WP1 itself" : "");
    }

    // ---------------- Braking burn at WP1 missed: free drift from the arrival state

    const Track tm = fly(n, {rWp1, vf}, {}, tSafe, kos, log, "WP1 braking missed", false);
    row("MIS-04c", "WP1 braking burn missed: no KOS entry", tm.tEnter < 0 ? "never enters, closest " + fmt(tm.minRange, 0) + " m"
                                                                          : "enters after " + fmt(tm.tEnter / 60.0) + " min",
        "within " + fmt(tSafe / 3600.0, 0) + " h", tm.tEnter < 0 ? "PASS" : "FLAG",
        "arrival speed " + fmt(norm(vf), 2) + " m/s" +
            (tm.tEnter < 0 ? "" : "; not passively safe: the WP1 braking burn is safety-critical (backed up by the "
                                  "1-fault-tolerant RCS; otherwise a CAM within that time)"));

    // ----------------- Final approach on the docking axis
    {
        const double wp2 = cfg.numValue("wp2_vbar");
        row("MIS-04b", "final approach inside the approach corridor", "on the +V-bar docking axis (0 deg)",
            "+-" + fmt(corridor / cnst::DEG, 0) + " deg: half-width " + fmt(kos * std::tan(corridor)) + " m at the KOS, " +
                fmt(wp2 * std::tan(corridor)) + " m at WP2",
            "PASS", "by design; lateral fine control = thrust-window driver 'single lateral pulse'");
    }

    // ------------------ Departure from the +V-bar port: push, DB0, DB1 along +x (posigrade)

    const double push = cfg.numValue("undock_push_dv"), db0 = cfg.numValue("dep_burn0_dv"), db1 = cfg.numValue("dep_burn1_dv");
    const double t0 = cfg.numValue("dep_burn0_time"), t1 = cfg.numValue("dep_burn1_time");
    const State start{Vec3(cfg.numValue("undock_start_range"), 0, 0), Vec3()};
    const Track dburnNominal = fly(n, start, {{0, push}, {t0, db0}, {t1, db1}}, tSafe, kos, log, "departure nominal", true);
    const Track dburn0_miss = fly(n, start, {{0, push}}, tSafe, kos, log, "departure DB0 missed", true);
    const Track dburn1_miss = fly(n, start, {{0, push}, {t0, db0}}, tSafe, kos, log, "departure DB1 missed", true);

    {
        const bool ok = dburnNominal.maxCone <= cone && dburnNominal.tExit >= 0 && dburnNominal.tEnter < 0;
        row("RCS-15a", "nominal departure: inside the cone, leaves the KOS, no re-entry",
            "max " + fmt(dburnNominal.maxCone / cnst::DEG) + " deg, exit after " + fmt(dburnNominal.tExit / 60.0) + " min, " +
                (dburnNominal.tEnter < 0 ? "no re-entry" : "RE-ENTERS"),
            "<= " + fmt(cone / cnst::DEG, 0) + " deg, " + fmt(tSafe / 3600.0, 0) + " h", ok ? "PASS" : "FAIL",
            "passes over the station at " + fmt(dburnNominal.overflight, 0) + " m");
    }

    {
        const bool ok = dburn0_miss.overflight < 0 || dburn0_miss.overflight > kos;
        row("RCS-15b", "DB0 missed (push only): passes over the station outside the KOS",
            "overflight at " + fmt(dburn0_miss.overflight, 0) + " m", "> " + fmt(kos, 0) + " m", ok ? "PASS" : "FLAG",
            ok ? "" : "not passively safe: DB0 is safety-critical (backed up by the 1-fault-tolerant RCS)");
    }

    row("RCS-15c", "DB1 missed (push + DB0): no KOS re-entry",
        dburn1_miss.tEnter < 0 ? "no re-entry, overflight at " + fmt(dburn1_miss.overflight, 0) + " m" : "re-enters after " + fmt(dburn1_miss.tEnter / 60.0) + " min",
        fmt(tSafe / 3600.0, 0) + " h", dburn1_miss.tEnter < 0 ? "PASS" : "FAIL", "");
    
    {
        const double used = db0 + db1, budget = cfg.numValue("departure_dv");
        row("RCS-15d", "departure burns inside the budget line", fmt(used, 2) + " m/s (DB0 + DB1)",
            "<= " + fmt(budget, 2) + " m/s (departure_dv)", used <= budget ? "PASS" : "FAIL", "");
    }
    
}