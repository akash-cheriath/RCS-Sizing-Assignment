// Entry points of the analyses done in this preliminary RCS sizing investigation.
//
// Analysis Flow : Same logic as the Orion orbit-control RCS trade:
// 1. runLayoutAnalysis     :   Checks what each architecture layout can do? 
//                              Largest pure force / torque axis (per newton of thrust) : all healthy, worst single, and worst double failure
//                              Stuck-on thruster check; Fault-tolerance level check
// 2a. runDvBudget          :   Checks how much delta-v / propellant is required for the mission, phase by phase, keeping ESA margins
// 2b. runProximityChecks   :   Checks the safety during approach and departure maneuvers (whether selected design dock precisely and departs safely)
// 3a. runThrustWindow      :   Measures the window of thrusts (smallest thrust meeting every "enough authority" driver and largest thrust
//                              that still allows fine control at docking)
// 3b. runGeometrySweeps    :   Run steps 1 - 3(a) again with pod ring moved axially and with each cant angle sweep : to check which position-
//                              cant angle combination gives the widest window
// 4. runComparison         :   Compare every architecture with thruster class to see if it fits the window? (Ranking and shortlisting of candidate
//                              architecture)

#pragma once
#include<string>
#include <vector>
#include "config_map.hpp"
#include "vehicle.hpp"

// The 12 "pure" demand directions used throughout : +-Fx, +-Fy, +-Fz, +-Tx, +-Ty, +-Tz
struct Demand {
    std::string name;
    Vec6 dem;
    bool torque;
};

inline std::vector<Demand> pureDemands() {
    std::vector<Demand> v;
    const char* nm[6] = {"Fx", "Fy", "Fz", "Tx", "Ty", "Tz"};
    for (int k = 0; k < 6; ++k) {
        for (int s = 1; s >= -1; s -= 2) {
            Vec6 dem{0, 0, 0, 0, 0, 0};
            dem[k] = s;
            v.push_back({std::string(s > 0 ? "+" : "-") + nm[k], dem, k >= 3});
        }
    }
    return v;
}

// Reentry velocity direction in body axes for trim angle alpha (lift towards +Z). The bank reversal is a
// rotation about this vector, not about the body X axis :
inline Vec3 reentryVelocityDir(double alpha) { 
    return Vec3(-std::cos(alpha), 0.0, std::sin(alpha)); 
}

// One candidate of the architecture.cfg built with F = 1 N per thruster
struct Architecture {
    std::string id;             // "01",...,"24"
    std::string family;         // "Single ring", "Double ring", "Dragon-like"
    std::string text;           // The config lines (for the tables)
    Layout layout;              // layout.name = e.g. "Q24"
};

// Step 1 Result, all "per newton of thrust" (authority scales linearly with the thrust)
struct LayoutResult {
    std::string name, family;
    int nThr = 0;
    // [mass configuration: 0 = combined vehicle, 1 = capsule][demand 0 ... 11 of pureDemands()]
    double nominal[2][12]     = {};
    double fail1[2][12]       = {};     // Worst case over all single failures
    double fail2[2][12]       = {};     // Worst case over all double failures
    double effBurn            = 0;      // -X efficiency of a main burn, combined vehicle (useful / fired thrust)
    double effAvg[2]          = {};     // Mean efficiency of the six pure forces (prox-ops pulsing)
    double torqueFree[2][6]   = {};     // Torque-only authority (forces free), nominal
    double torqueFree1[2][6]  = {};     // ... worst single failure
    double torqueFree2[2][6]  = {};     // ... worst double failure
    double bank[3]            = {};     // Capsule torque about the entry velocity vector: nominal, 1 fail, 2 fail
    double maxArm[3]          = {};     // Largest single-thruster torque per newton about X, Y, Z (combined vehicle)
    double maxRatePitchYaw    = 0;      // Largest single-thruster pitch/yaw rate change per N s, vector sum (combined vehicle)
    double axialBest[3]       = {};     // Better of pure +X / -X force: 0, 1, 2 failures (relaxed criterion)
    double stuckOnRatio       = 0;      // Smallest (counter-torque of the others / torque of a stuck-on thruster)
    bool oneFT                = false;  // Pure force and torque on all 12 axes after any 1 failure
    bool strict2FT            = false;  // ... after any 2 failures
    bool relaxed2FT           = false;  // Any 2 failed: torque in all 6 senses + bank + a pure axial force in one sense
};

// Step 2 Result
struct BudgetResult {
    double dv_total             = 0;    // [m/s] translational delta-v incl. margins
    double propellant_total     = 0;    // [kg] usable propellant
    double propellant_loaded    = 0;    // [kg] incl. residuals
    double mass_at_proxops      = 0;    // [kg] combined mass on arrival at WP0
    double mass_at_deorbit      = 0;    // [kg] combined mass at deorbit ignition
    double deorbit_burn_min     = 0;    // [min] nominal deorbit burn duration
};

// Step 3 Result: one thrust driver, and the window of one layout
struct ThrustDriver {
    std::string name, failure_case, basis;
    double thrust;                      // [N] thrust per thruster this driver needs (lower bound) or allows (upper bound)
    bool isMin;                         // true = lower bound
    double thrust_relax;                // the same driver under the RELAXED fault-tolerance criterion (= F for most drivers)
};

struct ThrustWindow {
    double F_min = 0, F_min_relaxed = 0, F_max = 0;     // [N]
    std::string driver_min, driver_min_relaxed, driver_max;
    std::vector<ThrustDriver> drivers;
};

// Step 4: Thruster classes of thruster_classes.cfg
struct ThrusterClass {
    std::string id;         // "A", "B", ...
    std::string note;       // comment of the config line
    ThrusterPerf perf;
    bool available = true;  // can be bought
    bool baseline = false;  // class_X = baseline (thruster of baseline_thruster.cfg)
};

// Step 4 Result of one architecture with one thruster class
struct ClassResult {
    const ThrusterClass* c = nullptr;
    BudgetResult budget;
    ThrustWindow window;
    std::string verdict = "no";         // "strict", "relaxed" or "no"
    std::string reason;                 // why not
    double fmin = 0;                    // F_min of the fault-tolerance level reached (strict, else relaxed)
};

inline std::string faultLevel(const LayoutResult& R) {
    return R.strict2FT ? "strict 2FT" : R.relaxed2FT ? "relaxed 2FT" : R.oneFT ? "1FT only" : "none";
}

// Step 1 (layout_analysis.cpp)
Architecture buildArchitecture(const std::string& id, const std::string& line, PodGeometry podGeom, double xBehind,
                               double xFront);
std::vector<Architecture> loadArchitectures(const Config& cfg, const Vehicle& veh);
LayoutResult analyseLayout(const Config& cfg, const Vehicle& veh, const Architecture& arch);
std::vector<LayoutResult> runLayoutAnalysis(const Config& cfg, const Vehicle& veh,
                                            const std::vector<Architecture>& archs, const std::string& outDir);

// Step 2a (dv_budget.cpp), step 2b (proximity_checks.cpp) and step 3a (thrust_window.cpp)
BudgetResult runDvBudget(const Config& cfg, const Vehicle& veh, const LayoutResult& layResult, const std::string& outDir);
void runProximityChecks(const Config& cfg, const std::string& outDir);
ThrustWindow runThrustWindow(const Config& cfg, const Vehicle& veh, const LayoutResult& layResult, const BudgetResult& layBudget);

// Step 3b (geometry_sweeps.cpp)
void runGeometrySweeps(const Config& cfg, const Vehicle& veh, const std::vector<Architecture>& archs,
                       const std::vector<LayoutResult>& layResults, const std::string& outDir);
// Step 4 (comparison.cpp)
std::vector<ThrusterClass> loadThrusterClasses(const Config& cfg, const Vehicle& veh);
ClassResult evaluateClass(const Config& cfg, const Vehicle& veh, const LayoutResult& layResult, const ThrusterClass& thrust_class);
void runComparison(const Config& cfg, const Vehicle& veh, const std::vector<Architecture>& archs,
                   const std::vector<LayoutResult>& layResults, const std::string& outDir);
