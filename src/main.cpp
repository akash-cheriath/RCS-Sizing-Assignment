// PRELIMINARY SIZING OF RCS FOR A CAPSULE + SERVICE MODULE

#include <chrono>
#include <cmath>
#include <iostream>
#include "analyses.hpp"
#include "csv_export.hpp"

// The service module is the difference between the combined vehicle and the capsule, so its CoG follows from a
// mass balance:  m_sm cg_sm = m_comb cg_comb - m_cap cg_cap.  If that point lies outside the SM body, the two
// given CoGs cannot both be right (or they refer to different reference points).
static void massPropertyCheck(const Config& cfg, const Vehicle& veh, const std::string& outDir) {
    const double mSm = veh.combined_vehicle.mass - veh.capsule.mass;
    const Vec3 cgSm = (veh.combined_vehicle.cog * veh.combined_vehicle.mass - veh.capsule.cog * veh.capsule.mass) / mSm;
    const double len = cfg.numValue("service_module_height"), rad = 0.5 * cfg.numValue("service_module_diameter");
    const bool inside = cgSm.x <= 0.0 && cgSm.x >= -len && std::hypot(cgSm.y, cgSm.z) <= rad;
    Csv out(outPath(outDir, "mass_properties_check.csv"));
    out.header({"sm_mass_kg", "sm_cg_x_m", "sm_cg_y_m", "sm_cg_z_m", "sm_body_x_min_m", "sm_body_x_max_m", "inside_sm_body"});
    out << mSm << cgSm.x << cgSm.y << cgSm.z << -len << 0.0 << (inside ? "yes" : "no");
    out.endl();
    std::cout << "  implied service-module CoG: x = " << cgSm.x << " m, radial " << std::hypot(cgSm.y, cgSm.z)
              << " m -> " << (inside ? "inside the SM body: consistent" : "OUTSIDE the SM body: CHECK the CoG inputs") << "\n";
}

int main(int argc, char** argv) {
    const std::string cfgFile = argc > 1 ? argv[1] : "../config/config.cfg";
    const std::string outDir = argc > 2 ? argv[2] : "../results";
    try {
        const auto t0 = std::chrono::steady_clock::now();
        const Config cfg(cfgFile);
        const Vehicle veh = loadVehicle(cfg);
        std::cout << "=== Preliminary RCS sizing ===\nconfig: " << cfgFile << "\nresults: " << outDir
                  << "\nbaseline thruster: " << veh.thr.F << " N, Isp " << veh.thr.isp << " s, t_min " << veh.thr.t_min * 1e3
                  << " ms\n\n[check] Mass properties\n";
        massPropertyCheck(cfg, veh, outDir);

        std::cout << "\n[1] Layout analysis (authority per newton, 1 and 2 failures)\n";
        const std::vector<Architecture> archs = loadArchitectures(cfg, veh);
        const std::vector<LayoutResult> results = runLayoutAnalysis(cfg, veh, archs, outDir);

        const std::string base = cfg.strValue("baseline_layout", "");
        std::cout << "\n[2] Delta-v and propellant budget (" << base << ", baseline thruster)\n";
        bool found = false;
        for (const auto& R : results)
            if (R.name == base) { runDvBudget(cfg, veh, R, outDir); found = true; }
        if (!found) throw std::runtime_error("baseline_layout '" + base + "' is not in architectures.cfg");

        std::cout << "\n[2b] Proximity operations: keep-out sphere, approach corridor, departure cone (closed-form CW)\n";
        runProximityChecks(cfg, outDir);

        std::cout << "\n[3b] Geometry sweeps (pod station, cant angles; steps 1-3 re-run, all architectures in parallel)\n";
        runGeometrySweeps(cfg, veh, archs, results, outDir);

        std::cout << "\n[3-4] Thrust window and comparison (every architecture x every thruster class)\n";
        runComparison(cfg, veh, archs, results, outDir);

        std::cout << "\nDone in " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s.\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
