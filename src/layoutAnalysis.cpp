// Checks what each candidate architecture can do? - How large a pure force / torque can be produced along each body axis
//
// Layouts built with F = 1 N
// 
//   * 12 pure demands (+-Fx ... +-Tz), all six rows constrained    -> Translation burns
//   * 6 torque-only demands, force rows free                       -> Attitude control
//   * Torque about the REENTRY VELOCITY VECTOR, capsule CoG        -> Bank reversal
//   * Every single and double failure (failed thruster held off)   -> Worst case kept
//   * Both mass configurations: combined vehicle (orbit) and capsule (entry)
//
// FAULT-TOLERANCE LEVELS
//   1FT     : Pure force and torque on all 12 axes after any 1 failure (prox-ops)
//   strict  : The same after any 2 failures (crew: CCT-REQ-1130 deorbit / entry,
//             SSP 51721 collision with the ISS = catastrophic -> 2-fault tolerant)
//   relaxed : After any 2 failures, torque in all 6 senses (forces free), bank
//             torque, and a pure axial force in at least ONE sense: the vehicle
//             turns to use what is left; docking is aborted

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "analyses.hpp"
#include "csv_export.hpp"

namespace {

const std::array<bool, 6> allForceTorque = {true, true, true, true, true, true};
const std::array<bool, 6> onlyTorque = {false, false, false, true, true, true};

// Parsing config lines 

// 1. "a ; b ; c" -> {"a", "b", "c"} (trimmed, empty pieces dropped)
std::vector<std::string> split(const std::string& str, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(str);
    for (std::string piece; std::getline(ss, piece, sep);) {
        const auto a = piece.find_first_not_of(" \t");
        const auto b = piece.find_last_not_of(" \t");
        if (a != std::string::npos)
            out.push_back(piece.substr(a, b - a +1));
    }
    return out;
}

// 2. "12 RCS thrusters at clocks=45,135 aft=1" -> {clocks: "45,135", aft: "1"}
std::map<std::string, std::string> keyValues(const std::string& str) {
    std::map<std::string, std::string> keyVal;
    for (const auto& token : split(str, ' ')) {
        const auto eq = token.find('=');
        if (eq != std::string::npos) 
            keyVal[token.substr(0, eq)] = token.substr(eq + 1);
    }
    return keyVal;
}

double numOf(const std::map<std::string, std::string>& keyVal, const std::string& key, double defaultVal) {
    auto it = keyVal.find(key);
    return it == keyVal.end() ? defaultVal : std::stod(it->second);
}

bool has(const std::string& str, const char* word) {
    return str.find(word) != std::string::npos;
}

// Failure loops

// 1. Worst case authority over all single (k = 1) or double (k = 2) failures
double worstCase(const Matrix& B, const Vec6& d, const std::array<bool, 6>& constrained, int k) {
    const int num_thr = static_cast<int>(B[0].size());
    double worst = 1e300;
    std::vector<bool> act(num_thr, true);
    for (int i = 0; i < num_thr; ++i) {
        act[i] = false;
        if (k == 1) {
            worst = std::min(worst, maxOutput(B, d, constrained, act));
        } else {
            for (int j = i + 1; j < num_thr; ++j) {
                act[j] = false;
                worst = std::min(worst, maxOutput(B, d, constrained, act));
                act[j] = true;
            }
        }
        act[i] = true;
    }
    return worst;
}

// 2. RELAXED criterion : after failures the vehicle may first TURN, so it needs only the BEST of the demands 'dem' (+X or -X).
// Worst case over all k failures (k = 0, 1, 2) of that best value
double worstCaseBest(const Matrix& B, const std::vector<Vec6>& dem, const std::array<bool, 6>& constrained, int k) {
    const int num_thr = static_cast<int>(B[0].size());
    std::vector<bool> act(num_thr, true);

    auto best = [&]() {
        double b = 0;
        for (const auto& d : dem)
            b = std::max(b, maxOutput(B, d, constrained, act));
        return b;
    };
    if (k == 0)
        return best();
    
    double worst = 1e300;
    for (int i = 0; i < num_thr; ++i) {
        act[i] = false;
        if (k == 1) {
            worst = std::min(worst, best());
        } else {
            for (int j = i + 1; j < num_thr; ++j) {
                act[j] = false;
                worst = std::min(worst, best());
                act[j] = true;
            }
        }
        act[i] = true;
    }
    return worst;
}

}   // namespace 

// xBehind and xFront are stations of the ring behind and infront of CoG respectively
Architecture buildArchitecture(const std::string& id, const std::string& line, PodGeometry podGeom,
                                double xBehind, double xFront) {
    const auto parts = split(line, ';');
    Architecture arch;
    arch.id = id;
    arch.text = line;
    arch.layout.name = parts.at(0);
    podGeom.force = 1.0;                            // Authority being checked per Newton
    const double pod_xSingle = podGeom.pod_x;
    bool twoRing = false;
    bool nose = false;
    int pods = 0;
    for (size_t p = 1; p < parts.size(); ++p) {
        const std::string& part = parts[p];
        const auto keyVal = keyValues(part);
        std::vector<double> clocks;
        for (const auto& c : split(keyVal.count("clocks") ? keyVal.at("clocks") : "", ','))
            clocks.push_back(std::stod(c) * cnst::DEG);
        if (clocks.empty())
            throw std::runtime_error("arch_" + id + " : a part without clocks = ");
        
        const size_t before = arch.layout.thr.size();
        if (has(part, "docking")) {
            nose = true;
            layouts::addNoseRing(arch.layout, podGeom, numOf(keyVal, "x", podGeom.nose_x), numOf(keyVal, "r", podGeom.nose_r),
                                clocks, static_cast<int>(numOf(keyVal, "fwd", 1)), pods);
        } 
        
        else {
            double x_pod = pod_xSingle;
            if (has(part, "behind")) {
                x_pod = xBehind;
                twoRing = true;
            }
            if (has(part, "in front")) {
                x_pod = xFront;
                twoRing = true;
            }
            podGeom.pod_x = numOf(keyVal, "x", x_pod);
            const int roll = static_cast<int> (numOf(keyVal, "roll", 0));
            const bool alternate = roll == 1;           // One roll jet per pod, ROLLP-ROLLN alternating
            if (alternate && clocks.size() % 2)
                throw std::runtime_error("arch_" + id + " : roll = 1 needs an even number of pods");
            const Layout Lay = layouts::Ring(podGeom, clocks, static_cast<int> (numOf(keyVal, "aft", 0)),
                                            static_cast<int> (numOf(keyVal, "fwd", 0)), alternate ? 1 : roll / 2,
                                            arch.layout.name, alternate);
            for (Thruster ithr : Lay.thr) {
                ithr.pod += pods;
                arch.layout.thr.push_back(ithr);
            }
        }
        // The parts start with its thruster count : Check it against what was built
        const int stated = std::atoi(part.c_str());
        const int built = static_cast<int>(arch.layout.thr.size() - before);
        if (stated != built) 
            throw std::runtime_error("arch_" + id + " (" + arch.layout.name + "): part says " + std::to_string(stated) +
                                     " thrusters, the keys build " + std::to_string(built));
        pods += static_cast<int>(clocks.size());
    }
    arch.family = twoRing ? "Double ring" : nose ? "Dragon-like" : "Single ring";
    return arch;
}

std::vector<Architecture> loadArchitectures(const Config& cfg, const Vehicle& veh) {
    std::vector<Architecture> archs;
    for (int i = 1; i <= 99; ++i) {
        const std::string id = (i < 10 ? "0" : "") + std::to_string(i);
        const std::string line = cfg.strValue("arch_" + id, "");
        if (!line.empty())
            archs.push_back(buildArchitecture(id, line, veh.podGeom, cfg.numValue("ring_behind_cog_x"),
                                              cfg.numValue("ring_front_cog_x")));
    }
    if (archs.empty()) 
        throw std::runtime_error("no arch_NN lines in the config");
    return archs;
}

LayoutResult analyseLayout(const Config& cfg, const Vehicle& veh, const Architecture& arch) {
    const Vec3 vReentry = reentryVelocityDir(cfg.numValue("alpha_trim") * cnst::DEG);
    const double throttle = cfg.numValue("main_burn_throttle");
    const auto demands = pureDemands();
    const MassProps* mp[2] = {&veh.combined_vehicle, &veh.capsule};
    const Layout& Lay = arch.layout;
    const int num_thr = Lay.size();
    LayoutResult LayRes;
    LayRes.name = Lay.name; 
    LayRes.family = arch.family; 
    LayRes.nThr = num_thr;

    const std::vector<double> thrustVec = thrustVector(Lay);
    const std::vector<bool> allActive(num_thr, true);
    double stuck = 1e300;
    for (int c = 0; c < 2; ++c) {
        const Matrix B = buildB(Lay, mp[c]->cog);

        // Pure force / torque demands (all six rows constrained)
        for (int k = 0; k < 12; ++k) {
            LayRes.nominal[c][k] = maxOutput(B, demands[k].dem, allForceTorque, allActive);
            LayRes.fail1[c][k] = worstCase(B, demands[k].dem, allForceTorque, 1);
            LayRes.fail2[c][k] = worstCase(B, demands[k].dem, allForceTorque, 2);
        }

        // Efficiency of the six pure forces: useful output / total thrust fired (cant + off-pulsing)
        double sum = 0;
        for (int k = 0; k < 6; ++k) 
            sum += efficientOutput(B, thrustVec, demands[k].dem, allForceTorque, LayRes.nominal[c][k]);
        LayRes.effAvg[c] = sum / 6.0;

        // Torque-only demands (forces free): attitude control
        for (int k = 6; k < 12; ++k) {
            LayRes.torqueFree[c][k - 6] = maxOutput(B, demands[k].dem, onlyTorque, allActive);
            LayRes.torqueFree1[c][k - 6] = worstCase(B, demands[k].dem, onlyTorque, 1);
            LayRes.torqueFree2[c][k - 6] = worstCase(B, demands[k].dem, onlyTorque, 2);
        }
        // Stuck-on thruster: can the others produce the opposite of its torque?
        for (int j = 0; j < num_thr; ++j) {
            const Vec3 trq(B[3][j], B[4][j], B[5][j]);
            const double tn = norm(trq);
            if (tn < 1e-9) 
                continue;
            std::vector<bool> act(num_thr, true);
            act[j] = false;
            const Vec6 d{0, 0, 0, -trq.x / tn, -trq.y / tn, -trq.z / tn};
            stuck = std::min(stuck, maxOutput(B, d, onlyTorque, act) / tn);
        }
    }

    LayRes.stuckOnRatio = stuck;

    {   
        // Main-burn efficiency: pure -X force at main_burn_throttle x its maximum (combined vehicle)
        const Matrix B = buildB(Lay, veh.combined_vehicle.cog);
        LayRes.effBurn = efficientOutput(B, thrustVec, demands[1].dem, allForceTorque, throttle * LayRes.nominal[0][1]);
        // Largest single-thruster torque arm about each axis (fine-control limit of step 3)
        for (int ax = 0; ax < 3; ++ax){
            for (int j = 0; j < num_thr; ++j) 
                LayRes.maxArm[ax] = std::max(LayRes.maxArm[ax], std::fabs(B[3 + ax][j]));
        }
               
        const double Iyy = veh.combined_vehicle.inertia.matrix[1][1];
        const double Izz = veh.combined_vehicle.inertia.matrix[2][2];
        for (int j = 0; j < num_thr; ++j)
            LayRes.maxRatePitchYaw = std::max(LayRes.maxRatePitchYaw, std::hypot(B[4][j] / Iyy, B[5][j] / Izz));
    
        // Relaxed criterion: better of pure +X / -X after 0, 1, 2 failures
        const std::vector<Vec6> axial = {Vec6{1, 0, 0, 0, 0, 0}, Vec6{-1, 0, 0, 0, 0, 0}};
        for (int k = 0; k < 3; ++k) 
            LayRes.axialBest[k] = worstCaseBest(B, axial, allForceTorque, k);
    }

    {   
        // Capsule torque about the reentry velocity vector (bank), both senses
        const Matrix B = buildB(Lay, veh.capsule.cog);
        double best[3] = {1e300, 1e300, 1e300};
        for (int s = 1; s >= -1; s -= 2) {
            const Vec6 d{0, 0, 0, s * vReentry.x, s * vReentry.y, s * vReentry.z};
            best[0] = std::min(best[0], maxOutput(B, d, onlyTorque, allActive));
            best[1] = std::min(best[1], worstCase(B, d, onlyTorque, 1));
            best[2] = std::min(best[2], worstCase(B, d, onlyTorque, 2));
        }
        for (int i = 0; i < 3; ++i) 
            LayRes.bank[i] = best[i];
    }

    // Fault-tolerance levels
    LayRes.oneFT = LayRes.strict2FT = true;
    for (int k = 0; k < 12; ++k) {
        LayRes.oneFT = LayRes.oneFT && LayRes.fail1[0][k] > 1e-6;
        LayRes.strict2FT = LayRes.strict2FT && LayRes.fail2[0][k] > 1e-6;
    }
    bool torques = true;
    for (int c = 0; c < 2; ++c)
        for (int k = 0; k < 6; ++k) 
            torques = torques && LayRes.torqueFree2[c][k] > 1e-6;
    LayRes.relaxed2FT = torques && LayRes.bank[2] > 1e-6 && LayRes.axialBest[2] > 1e-6;
    return LayRes;
}


std::vector<LayoutResult> runLayoutAnalysis(const Config& cfg, const Vehicle& veh,
                                            const std::vector<Architecture>& archs, const std::string& outDir) {
    const auto demands = pureDemands();
    const char* cfgName[2] = {"combined", "capsule"};

    Csv thrust(outPath(outDir, "layout_thrusters.csv"));
    thrust.header({"arch_id", "layout", "name", "type", "pod", "string", "x", "y", "z", "dir_x", "dir_y", "dir_z"});
    Csv auth(outPath(outDir, "layout_authority.csv"));
    auth.header({"arch_id", "layout", "config", "demand", "type", "nominal_per_N", "fail1_per_N", "fail2_per_N"});
    Csv summ(outPath(outDir, "layout_summary.csv"));
    summ.header({"arch_id", "layout", "family", "n_thr", "n_aft", "n_fwd", "n_roll", "n_nose", "fault_tolerance",
               "all_axes_1FT", "strict_2FT", "relaxed_2FT", "minusX_per_N", "minusX_1fail_per_N", "minusX_2fail_per_N",
               "best_axial_2fail_per_N", "bank_per_N", "bank_2fail_per_N", "main_burn_eff", "stuck_on_ratio",
               "max_arm_x", "max_arm_y", "max_arm_z"});

    std::vector<LayoutResult> results;
    for (const Architecture& a : archs) {
        const Layout& Lay = a.layout;
        std::cout << "  " << a.id << " " << Lay.name << " (" << Lay.size() << " thrusters, " << a.family << ") ..." << std::flush;
        const LayoutResult LayRes = analyseLayout(cfg, veh, a);

        int nType[4] = {0, 0, 0, 0};                    // AFT, FWD, ROLL, NOSE
        for (const auto& t : Lay.thr) {
            thrust << a.id << Lay.name << t.name << t.type << t.pod << t.string_id << t.pos.x << t.pos.y << t.pos.z
               << t.dir.x << t.dir.y << t.dir.z;
            thrust.endl();
            nType[t.type == "AFT" ? 0 : t.type == "FWD" ? 1 : t.type == "NOSE" ? 3 : 2]++;
        }
        for (int c = 0; c < 2; ++c) {
            for (int k = 0; k < 12; ++k) {
                auth << a.id << Lay.name << cfgName[c] << demands[k].name << "pure" << LayRes.nominal[c][k] << LayRes.fail1[c][k]
                   << LayRes.fail2[c][k];
                auth.endl();
            }
            for (int k = 6; k < 12; ++k) {
                auth << a.id << Lay.name << cfgName[c] << demands[k].name << "torque_only" << LayRes.torqueFree[c][k - 6]
                   << LayRes.torqueFree1[c][k - 6] << LayRes.torqueFree2[c][k - 6];
                auth.endl();
            }
        }
        const std::string level = LayRes.strict2FT ? "strict 2FT" : LayRes.relaxed2FT ? "relaxed 2FT" : LayRes.oneFT ? "1FT only" : "none";
        summ << a.id << Lay.name << a.family << Lay.size() << nType[0] << nType[1] << nType[2] << nType[3] << level
           << (LayRes.oneFT ? "yes" : "no") << (LayRes.strict2FT ? "yes" : "no") << (LayRes.relaxed2FT ? "yes" : "no")
           << LayRes.nominal[0][1] << LayRes.fail1[0][1] << LayRes.fail2[0][1] << LayRes.axialBest[2] << LayRes.bank[0] << LayRes.bank[2]
           << LayRes.effBurn << LayRes.stuckOnRatio << LayRes.maxArm[0] << LayRes.maxArm[1] << LayRes.maxArm[2];
        summ.endl();
        std::cout << " " << level << "; -X " << LayRes.nominal[0][1] << " N/N (2 fail " << LayRes.fail2[0][1] << "), bank "
                  << LayRes.bank[0] << " Nm/N (2 fail " << LayRes.bank[2] << "), stuck-on ratio " << LayRes.stuckOnRatio << "\n";
        results.push_back(LayRes);
    }
    return results;
}

