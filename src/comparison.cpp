// Compares different thruster architecture

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "analyses.hpp"
#include "csv_export.hpp"

namespace {

std::string show(double F) {                    // thrust for the console: "inf" for a lost axis
    if (F > 1e6) return "inf";
    std::ostringstream o;
    o << static_cast<long>(std::lround(F));
    return o.str();
}

struct Row {
    int a;                                      // index of the architecture
    ClassResult result;
};

}  // namespace

// Thruster classes: "class_A = F=111 isp=280 isp_pulse=252 t_min=8 available=True", or "class_B = baseline ..."
std::vector<ThrusterClass> loadThrusterClasses(const Config& cfg, const Vehicle& veh) {
    std::vector<ThrusterClass> out;
    for (char c = 'A'; c <= 'Z'; ++c) {
        const std::string line = cfg.strValue(std::string("class_") + c, "");
        if (line.empty()) 
            continue;
        ThrusterClass thrClass;
        thrClass.id = std::string(1, c);
        std::istringstream is(line);
        for (std::string tok; is >> tok;) {
            const auto eq = tok.find('=');
            if (eq == std::string::npos) { 
                thrClass.baseline = thrClass.baseline || tok == "baseline"; 
                continue; 
            }

            const std::string k = tok.substr(0, eq), v = tok.substr(eq + 1);
            if (k == "F") 
                thrClass.perf.F = std::stod(v);
            else if (k == "isp") 
                thrClass.perf.isp = std::stod(v);
            else if (k == "isp_pulse") 
                thrClass.perf.isp_pulse = std::stod(v);
            else if (k == "t_min") 
                thrClass.perf.t_min = std::stod(v) * 1e-3;
            else if (k == "available") 
                thrClass.available = v == "True" || v == "true" || v == "yes";
        }
        if (thrClass.baseline) thrClass.perf = veh.thr;     // thrust, Isp and t_min of baseline_thruster.cfg
        out.push_back(thrClass);
    }
    if (out.empty()) 
        throw std::runtime_error("no class_X lines in the config");
    return out;
}

// One architecture with one thruster class: Budget (step 2), Window (step 3), Verdict.

ClassResult evaluateClass(const Config& cfg, const Vehicle& veh, const LayoutResult& LayRes, const ThrusterClass& ThrClass) {
    Vehicle v = veh;
    v.thr = ThrClass.perf;
    v.podGeom.force = ThrClass.perf.F;
    ClassResult res;
    res.c = &ThrClass;
    res.budget = runDvBudget(cfg, v, LayRes, "");
    res.window = runThrustWindow(cfg, v, LayRes, res.budget);
    res.fmin = LayRes.strict2FT ? res.window.F_min : res.window.F_min_relaxed;
    const double F = ThrClass.perf.F;
    const double stuckMin = cfg.numValue("stuck_on_ratio_min", 1.5);                // default 1.5 if the key is absent
    if (!LayRes.strict2FT && !LayRes.relaxed2FT) 
        res.reason = "not 2-fault tolerant";
    else if (LayRes.stuckOnRatio < stuckMin) {
        std::ostringstream out;
        out << "stuck-on ratio below " << stuckMin;
        res.reason = out.str();
    }
    else if (F > res.window.F_max) 
        res.reason = "above F_max " + show(res.window.F_max) + " N";
    else if (F < res.fmin) 
        res.reason = "below F_min " + show(res.fmin) + " N";
    else 
        res.verdict = LayRes.strict2FT ? "strict" : "relaxed";
    return res;
}

void runComparison(const Config& cfg, const Vehicle& veh, const std::vector<Architecture>& archs,
                   const std::vector<LayoutResult>& results, const std::string& outDir) {
    const std::vector<ThrusterClass> classes = loadThrusterClasses(cfg, veh);
    const ThrusterClass* base = &classes[0];
    for (const auto& c : classes) 
        if (c.baseline) 
            base = &c;
    std::cout << "  stuck-on ratio limit " << cfg.numValue("stuck_on_ratio_min", 1.5)
              << (cfg.strValue("stuck_on_ratio_min", "").empty() ? " (default: key stuck_on_ratio_min not in the config)" : "")
              << "\n";

    // ---- Every architecture x every class

    std::vector<Row> rows;
    for (int i = 0; i < static_cast<int>(archs.size()); ++i)
        for (const auto& c : classes) 
            rows.push_back({i, evaluateClass(cfg, veh, results[i], c)});

    auto rowOf = [&](int a, const ThrusterClass* c) -> const ClassResult& {
        for (const auto& r : rows) 
            if (r.a == a && r.result.c == c) 
                return r.result;
        throw std::runtime_error("comparison row not found");
    };

    // ---- Comparison.csv
    Csv cmp(outPath(outDir, "comparison.csv"));
    cmp.header({"arch_id", "layout", "family", "n_thr", "fault_tolerance", "class", "F_N", "isp_s", "t_min_ms",
                "available", "F_min_N", "driver_min", "F_min_relaxed_N", "driver_min_relaxed", "F_max_N", "driver_max",
                "verdict", "reason", "thrust_margin_over_F_min", "loaded_propellant_kg", "deorbit_burn_min"});
    for (const auto& row : rows) {
        const LayoutResult& LayRes = results[row.a];
        const ClassResult& res = row.result;
        cmp << archs[row.a].id << LayRes.name << LayRes.family << LayRes.nThr << faultLevel(LayRes) << res.c->id << res.c->perf.F << res.c->perf.isp
            << res.c->perf.t_min * 1e3 << (res.c->available ? "yes" : "no") << res.window.F_min << res.window.driver_min
            << res.window.F_min_relaxed << res.window.driver_min_relaxed << res.window.F_max << res.window.driver_max << res.verdict << res.reason
            << res.c->perf.F / res.fmin - 1.0 << res.budget.propellant_loaded << res.budget.deorbit_burn_min;
        cmp.endl();
    }

    // ---- Thrust_drivers.csv (baseline class)
    Csv drv(outPath(outDir, "thrust_drivers.csv"));
    drv.header({"arch_id", "layout", "driver", "bound", "failure_case", "F_N", "F_relaxed_N", "basis"});
    for (int i = 0; i < static_cast<int>(archs.size()); ++i)
        for (const auto& d : rowOf(i, base).window.drivers) {
            drv << archs[i].id << results[i].name << d.name << (d.isMin ? "min" : "max") << d.failure_case << d.thrust
                << d.thrust_relax << d.basis;
            drv.endl();
        }

    // ---- Console: the matrix architecture x class
    std::cout << "  cell = verdict and loaded propellant [kg] (S strict, R relaxed) or why not\n  layout   n   ";
    for (const auto& c : classes) std::cout << std::left << std::setw(20) << ("class " + c.id + " (" + show(c.perf.F) + " N)");
    std::cout << "\n";
    for (int i = 0; i < static_cast<int>(archs.size()); ++i) {
        std::cout << "  " << std::left << std::setw(7) << results[i].name << std::setw(4) << results[i].nThr;
        for (const auto& c : classes) {
            const ClassResult& r = rowOf(i, &c);
            std::string cell = r.verdict == "no" ? "- " + r.reason
                                                 : (r.verdict == "strict" ? "S " : "R ") + show(r.budget.propellant_loaded);
            std::cout << std::setw(20) << cell.substr(0, 19);
        }
        std::cout << "\n";
    }
    std::cout << std::right;

    // ---- Ranking of the qualified combinations (classes that can be bought)
    std::vector<const Row*> q;
    for (const auto& row : rows) 
        if (row.result.verdict != "no" && row.result.c->available) 
            q.push_back(&row);
    
    std::sort(q.begin(), q.end(), [&](const Row* a, const Row* b) {
        if ((a->result.verdict == "strict") != (b->result.verdict == "strict")) 
            return a->result.verdict == "strict";
        if (results[a->a].nThr != results[b->a].nThr) 
            return results[a->a].nThr < results[b->a].nThr;
        return a->result.budget.propellant_loaded < b->result.budget.propellant_loaded;
    });

    std::cout << "  qualified combinations (available classes): " << q.size() << "\n";
    for (size_t k = 0; k < q.size() && k < 8; ++k)
        std::cout << "    " << k + 1 << ". " << results[q[k]->a].name << " x class " << q[k]->result.c->id << " ("
                  << show(q[k]->result.c->perf.F) << " N): " << q[k]->result.verdict << ", " << results[q[k]->a].nThr
                  << " thrusters, window " << show(q[k]->result.fmin) << "-" << show(q[k]->result.window.F_max) << " N, "
                  << show(q[k]->result.budget.propellant_loaded) << " kg\n";

    
    // ---- Architecture ranking (baseline class: window and propellant independent of the class list)
    std::vector<int> order(archs.size());
    for (size_t i = 0; i < order.size(); ++i) 
        order[i] = static_cast<int>(i);

    auto ftRank = [&](int i) { 
        return results[i].strict2FT ? 0 : results[i].relaxed2FT ? 1 : 2; 
    };

    auto isOpen = [&](int i) { 
        const ClassResult& res = rowOf(i, base); 
        return res.fmin <= res.window.F_max; 
    };

    auto windowRatio = [&](int i) {         // F_max / F_min with the baseline thruster (>= 1: window open)
        const ClassResult& res = rowOf(i, base);
        return res.window.F_max / std::max(res.fmin, 1e-9);
    };

    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (ftRank(a) != ftRank(b)) 
            return ftRank(a) < ftRank(b);
        if (isOpen(a) != isOpen(b)) 
            return isOpen(a);
        if (!isOpen(a)) {                   // both closed: the one closest to opening first (1 % steps)
            const long ra = std::lround(100.0 * windowRatio(a));
            const long rb = std::lround(100.0 * windowRatio(b));
            if (ra != rb)
                return ra > rb;
        }
        if (results[a].nThr != results[b].nThr) 
            return results[a].nThr < results[b].nThr;
        return rowOf(a, base).budget.propellant_loaded < rowOf(b, base).budget.propellant_loaded;
    });

    auto fitting = [&](int i) {             // classes inside the window ("*" = cannot be bought)
        std::string s;
        for (const auto& c : classes)
            if (rowOf(i, &c).verdict != "no") 
                s += (s.empty() ? "" : " ") + c.id + (c.available ? "" : "*");
        return s.empty() ? std::string("none") : s;
    };

    Csv rk(outPath(outDir, "architecture_ranking.csv"));
    rk.header({"rank", "arch_id", "layout", "family", "n_thr", "fault_tolerance", "F_min_N", "F_max_N", "window_open",
               "driver_min", "driver_max", "loaded_propellant_kg_baseline_thruster", "stuck_on_ratio", "classes_in_window"});
    for (size_t k = 0; k < order.size(); ++k) {
        const int i = order[k];
        const ClassResult& r = rowOf(i, base);
        rk << static_cast<int>(k + 1) << archs[i].id << results[i].name << results[i].family << results[i].nThr
           << faultLevel(results[i]) << r.fmin << r.window.F_max << (isOpen(i) ? "yes" : "no")
           << (results[i].strict2FT ? r.window.driver_min : r.window.driver_min_relaxed) << r.window.driver_max
           << r.budget.propellant_loaded << results[i].stuckOnRatio << fitting(i);
        rk.endl();
    }

    // ---- Shortlist: best of each family + the architecture of the best qualified combination
    std::vector<std::pair<int, std::string>> shortlist;

    auto inList = [&](int i) { 
        for (const auto& s : shortlist) 
            if (s.first == i) 
                return true; 
        return false; 
    };

    for (const char* fam : {"Single ring", "Double ring", "Dragon-like"})
        for (int i : order)
            if (results[i].family == fam) { 
                shortlist.push_back({i, std::string("best ") + fam}); 
                break; 
            }

    if (!q.empty()) {
        const int i = q.front()->a;
        if (!inList(i)) 
            shortlist.push_back({i, "best qualified combination"});
        else 
            for (auto& s : shortlist) 
                if (s.first == i) 
                    s.second += " + best qualified combination";
    }

    Csv sl(outPath(outDir, "shortlist.csv"));
    sl.header({"arch_id", "layout", "family", "selected_as", "n_thr", "fault_tolerance", "layout_text",
               "F_min_N", "F_max_N", "window_open", "driver_min", "driver_max", "loaded_propellant_kg_baseline_thruster",
               "best_class", "best_class_F_N", "best_class_F_max_N", "best_class_loaded_propellant_kg", "classes_in_window",
               "stuck_on_ratio"});
    std::cout << "  SHORTLIST (window with the baseline thruster, t_min " << base->perf.t_min * 1e3 << " ms):\n";
    for (const auto& s : shortlist) {
        const int i = s.first;
        const ClassResult& r = rowOf(i, base);
        const ClassResult* best = nullptr;                  // Least propellant among the fitting available classes
        for (const auto& c : classes) {
            const ClassResult& x = rowOf(i, &c);
            if (x.verdict != "no" && c.available && (!best || x.budget.propellant_loaded < best->budget.propellant_loaded)) 
                best = &x;
        }
        const std::string text = archs[i].text.substr(archs[i].text.find(';') + 1);

        sl << archs[i].id << results[i].name << results[i].family << s.second << results[i].nThr << faultLevel(results[i])
           << text << r.fmin << r.window.F_max << (isOpen(i) ? "yes" : "no")
           << (results[i].strict2FT ? r.window.driver_min : r.window.driver_min_relaxed) << r.window.driver_max << r.budget.propellant_loaded
           << (best ? best->c->id : std::string("none")) << (best ? best->c->perf.F : 0.0) << (best ? best->window.F_max : 0.0)
           << (best ? best->budget.propellant_loaded : 0.0) << fitting(i) << results[i].stuckOnRatio;

        sl.endl();
        std::cout << "    " << std::left << std::setw(6) << results[i].name << std::right << " (" << results[i].nThr << ", "
                  << faultLevel(results[i]) << ", " << s.second << "): window " << show(r.fmin) << "-" << show(r.window.F_max)
                  << " N" << (isOpen(i) ? "" : " (closed)") << ", " << show(r.budget.propellant_loaded)
                  << " kg with the baseline thruster; fitting classes: " << fitting(i);
        if (best) std::cout << " -> class " << best->c->id << " (" << show(best->c->perf.F) << " N): window "
                            << show(best->fmin) << "-" << show(best->window.F_max) << " N, " << show(best->budget.propellant_loaded) << " kg";
        std::cout << "\n";
    }
}
