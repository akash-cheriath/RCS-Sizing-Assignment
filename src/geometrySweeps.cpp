// Sweeps across different values of axial distance of the pods as well as the cant angles of the thrusters
//
#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <sstream>
#include "analyses.hpp"
#include "csv_export.hpp"

namespace {

struct Variant {
    std::string trade;          // "station", "cant_fwd", "cant_aft", "cant_roll"
    double value;               // [m] or [deg]
    bool admissible;
    bool isConfig;              // the value of the config
};

struct Point {
    Variant var;
    std::string faultTol;             // fault-tolerance level reached
    int faultTolRank;                 // 0 strict, 1 relaxed, 2 1FT only, 3 none
    double fmin, fmax, ratio, kg;
    std::string fits;                 // classes inside the window
    std::string fitsBuyable;          // ids of the classes that fit and can be bought, e.g. "A"
};

int rankOf(const LayoutResult& LayRes) { 
    return LayRes.strict2FT ? 0 : LayRes.relaxed2FT ? 1 : LayRes.oneFT ? 2 : 3; 
}

// lo, lo + step, ... (<= hi), plus hi and the config value; sorted, no duplicates
std::vector<double> grid(double lo, double hi, double step, double cfgValue) {
    std::vector<double> v;
    for (int k = 0; lo + k * step <= hi + 1e-9; ++k) 
        v.push_back(std::round((lo + k * step) * 1e6) / 1e6);
    for (double x : {hi, cfgValue}) {
        if (std::none_of(v.begin(), v.end(), [&](double y) { return std::fabs(x - y) < 1e-6; })) 
            v.push_back(x);
    }
    std::sort(v.begin(), v.end());
    return v;
}

std::string num(double v) { 
    std::ostringstream out; 
    out << v; 
    return out.str(); 
}

// All variants of one architecture (runs in its own thread)
std::vector<Point> sweepArchitecture(const Config& cfg, const Vehicle& veh, const Architecture& arch,
                                     const std::vector<ThrusterClass>& classes) {
    const double xBehind    = cfg.numValue("ring_behind_cog_x"); 
    const double xFront     = cfg.numValue("ring_front_cog_x");
    const bool twoRing      = arch.family == "Double ring";
    const double aftMin     = cfg.numValue("sidewall_half_angle") + 12.0;   // plume-clearance rule (wall, SM)
    const double fwdMin     = cfg.numValue("cant_fwd_min_plume");           // plume away from the station

    std::vector<Variant> vars;
    const double cfgX = twoRing ? xFront : veh.podGeom.pod_x;
    const std::vector<double> xs = twoRing
        ? grid(cfg.numValue("front_ring_x_min"), cfg.numValue("front_ring_x_max"), cfg.numValue("front_ring_x_step"), cfgX)
        : grid(cfg.numValue("pod_x_min"), cfg.numValue("pod_x_max"), cfg.numValue("pod_trade_step"), cfgX);

    for (double x : xs) 
        vars.push_back({"station", x, true, std::fabs(x - cfgX) < 1e-6});

    const struct { const char* name; const char* key; double cfgDeg; } cants[] = {
        {"cant_fwd", "cant_trade_fwd", veh.podGeom.cant_fwd / cnst::DEG},
        {"cant_aft", "cant_trade_aft", veh.podGeom.cant_aft / cnst::DEG},
        {"cant_roll", "cant_trade_roll", veh.podGeom.cant_roll / cnst::DEG}
    };

    for (const auto& c : cants)
        for (double d : grid(cfg.numValue(std::string(c.key) + "_min"), cfg.numValue(std::string(c.key) + "_max"),
                             cfg.numValue("cant_trade_step"), c.cfgDeg))
            vars.push_back({c.name, d,
                            !(std::string(c.name) == "cant_aft" && d < aftMin - 1e-9) &&
                                !(std::string(c.name) == "cant_fwd" && d < fwdMin - 1e-9),
                            std::fabs(d - c.cfgDeg) < 1e-6});

    std::vector<Point> out;
    for (const Variant& v : vars) {
        PodGeometry podGeom = veh.podGeom;
        double xf = xFront;
        if (v.trade == "station") { 
            if (twoRing) 
                xf = v.value; 
            else 
                podGeom.pod_x = v.value; 
        }
        if (v.trade == "cant_fwd") 
            podGeom.cant_fwd = v.value * cnst::DEG;
        if (v.trade == "cant_aft") 
            podGeom.cant_aft = v.value * cnst::DEG;
        if (v.trade == "cant_roll") 
            podGeom.cant_roll = v.value * cnst::DEG;

        const LayoutResult LayRes = analyseLayout(cfg, veh, buildArchitecture(arch.id, arch.text, podGeom, xBehind, xf));   

        Point pnt {v, faultLevel(LayRes), rankOf(LayRes), 0, 0, 0, 0, "", ""};
        for (const auto& c : classes) {
            const ClassResult result = evaluateClass(cfg, veh, LayRes, c);                                            // steps 2-3
            if (result.verdict != "no") 
                pnt.fits += (pnt.fits.empty() ? "" : " ") + c.id;
            if (result.verdict != "no" && c.available) 
                pnt.fitsBuyable += c.id;
            if (c.baseline) { 
                pnt.fmin = result.fmin; 
                pnt.fmax = result.window.F_max; 
                pnt.kg = result.budget.propellant_loaded; 
            }
        }
        pnt.ratio = pnt.fmax / pnt.fmin;
        if (pnt.fits.empty()) pnt.fits = "none";
        out.push_back(pnt);
    }
    return out;
}

}  // namespace

void runGeometrySweeps(const Config& cfg, const Vehicle& veh, const std::vector<Architecture>& archs,
                       const std::vector<LayoutResult>& results, const std::string& outDir) {

    if (cfg.strValue("run_geometry_sweeps", "yes") != "yes") { 
        std::cout << "  skipped (run_geometry_sweeps = no)\n"; 
        return; 
    }

    const std::vector<ThrusterClass> classes = loadThrusterClasses(cfg, veh);

    // ---- One thread per architecture

    std::vector<std::future<std::vector<Point>>> jobs;
    for (const auto& a : archs)
        jobs.push_back(std::async(std::launch::async, sweepArchitecture, std::cref(cfg), std::cref(veh), std::cref(a),
                                  std::cref(classes)));

    std::vector<std::vector<Point>> sweeps;
    for (auto& j : jobs) 
        sweeps.push_back(j.get());

    Csv stationSweep(outPath(outDir, "station_sweep.csv"));
    stationSweep.header({"arch_id", "layout", "family", "swept_ring", "x_m", "fault_tolerance", "F_min_N", "F_max_N", "window_ratio",
               "window_open", "loaded_propellant_kg", "classes_in_window", "config_value", "best"});
               
    Csv cantSweep(outPath(outDir, "cant_sweep.csv"));
    cantSweep.header({"arch_id", "layout", "family", "nozzle", "cant_deg", "admissible", "fault_tolerance", "F_min_N", "F_max_N",
               "window_ratio", "window_open", "loaded_propellant_kg", "classes_in_window", "config_value", "best"});

    Csv geomBest(outPath(outDir, "geometry_best.csv"));
    geomBest.header({"arch_id", "layout", "family", "parameter", "config_value", "best_value", "config_is_best",
               "window_ratio_config", "window_ratio_best", "loaded_propellant_kg_config", "loaded_propellant_kg_best",
               "classes_in_window_config", "classes_in_window_best"});

    const double gain = cfg.numValue("sweep_clear_gain");
    int sameCount = 0;
    int total = 0;
    for (size_t i = 0; i < archs.size(); ++i) {
        const auto& pts = sweeps[i];
        const int r0 = rankOf(results[i]);                              // fault-tolerance level with the config geometry
        std::string line;
        for (const char* trade : {"station", "cant_fwd", "cant_aft", "cant_roll"}) {
            const Point* cfgPt = nullptr;
            for (const auto& p : pts) 
                if (p.var.trade == trade && p.var.isConfig) 
                    cfgPt = &p;
            if (!cfgPt) 
                continue;

            // A value may replace the config value only if it is admissible, keeps the fault tolerance and the
            // buyable thrusters that fit now, and is clearly better (wider window, or less propellant)

            auto keepsFits = [&](const Point& p) {
                for (char c : cfgPt->fitsBuyable) 
                    if (p.fitsBuyable.find(c) == std::string::npos) 
                        return false;
                return true;
            };
            
            auto clearlyBetter = [&](const Point& p) {
                return p.var.admissible && p.faultTolRank <= r0 && keepsFits(p) &&
                       (p.ratio >= (1.0 + gain) * cfgPt->ratio ||
                        (p.ratio >= 0.99 * cfgPt->ratio && p.kg <= (1.0 - gain) * cfgPt->kg));
            };

            const Point* best = cfgPt;
            for (const auto& p : pts) {
                if (p.var.trade != trade || !clearlyBetter(p)) 
                    continue;
                if (best == cfgPt || p.ratio > 1.01 * best->ratio ||
                    (p.ratio >= 0.99 * best->ratio && p.kg < best->kg))
                    best = &p;
            }

            for (const auto& p : pts) {
                if (p.var.trade != trade) 
                    continue;

                const bool isBest = &p == best;
                if (std::string(trade) == "station") {
                    stationSweep << archs[i].id << results[i].name << results[i].family
                       << (results[i].family == "Double ring" ? "front ring (aft ring at " + cfg.strValue("ring_behind_cog_x", "") + " m)"
                           : results[i].family == "Dragon-like" ? std::string("sidewall ring") : std::string("Single ring"))
                       << p.var.value << p.faultTol << p.fmin << p.fmax << p.ratio << (p.ratio >= 1.0 ? "yes" : "no") << p.kg
                       << p.fits << (p.var.isConfig ? "yes" : "no") << (isBest ? "yes" : "no");
                    stationSweep.endl();
                } 
                else {
                    cantSweep << archs[i].id << results[i].name << results[i].family << std::string(trade).substr(5) << p.var.value
                       << (p.var.admissible ? "yes" : "no") << p.faultTol << p.fmin << p.fmax << p.ratio
                       << (p.ratio >= 1.0 ? "yes" : "no") << p.kg << p.fits << (p.var.isConfig ? "yes" : "no")
                       << (isBest ? "yes" : "no");
                    cantSweep.endl();
                }
            }

            const bool same = best == cfgPt;
            sameCount += same; ++total;
            geomBest << archs[i].id << results[i].name << results[i].family << trade << cfgPt->var.value << best->var.value
               << (same ? "yes" : "no") << cfgPt->ratio << best->ratio << cfgPt->kg << best->kg << cfgPt->fits << best->fits;
            geomBest.endl();
            line += std::string("  ") + (std::string(trade) == "station" ? "x" : std::string(trade).substr(5)) + " " +
                    num(cfgPt->var.value) + (same ? "" : "->" + num(best->var.value)) + " (" +
                    num(std::round(best->ratio * 100) / 100) + ")";
        }
        std::cout << "  " << archs[i].id << " " << results[i].name << ":" << line << "\n";
    }
    std::cout << "  (config value -> best value, window width F_max/F_min at the best value; > 1 = open)\n"
              << "  config value is the best in " << sameCount << " of " << total << " sweeps\n";
}


