#pragma once
#include <string>
#include <map>
#include <fstream>
#include <stdexcept>


class Config{
public:
    // To avoid causing ERROR when multiple config files have same key names. The keys will be OVERRIDEN by the
    // file called further
    explicit Config(const std::string& file) { load(file, false); }
    void merge(const std::string& file) { load(file, true); }

    double numValue(const std::string& key) const {
        auto it = configMap.find(key);
        if (it != configMap.end()) {
            return std::stod(it->second);
        }
        throw std::runtime_error("Key not found in configuration: " + key);
    }

    double numValue(const std::string& key, double defaultValue) const {
        auto it = configMap.find(key);
        if (it != configMap.end()) {
            return std::stod(it->second);
        }
        return defaultValue;
    }

    std::string strValue(const std::string& key, const std::string& defaultValue) const {
        auto it = configMap.find(key);
        if (it != configMap.end()) {
            return it->second;
        }
        return defaultValue;
    }

private:

    std::map<std::string, std::string> configMap;
    std::map<std::string, std::string> from_ ;      // file each key was read from
    
    void load(const std::string& file, bool isOverride){
        std::ifstream in(file);
        if (!in) throw std::runtime_error("Cannot open config file: " + file);
        const auto slash = file.find_last_of("/\\");
        const std::string folder = slash == std::string::npos ? "" : file.substr(0, slash + 1);
        std::string line;
        while (std::getline(in, line)){
            auto hash = line.find('#');
            if (hash != std::string::npos) line = line.substr(0, hash);
            auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string key = trim(line.substr(0, eq));
            std::string val = trim(line.substr(eq + 1));
            if (key.empty()) continue;
            if (key == "include") { 
                load(folder + val, isOverride); 
                continue;
            }
            auto it = from_.find(key);
            if (!isOverride && it != from_.end())
                throw std::runtime_error("Config key '" + key + "' is defined twice in " + it->second + " and in " + file);
            configMap[key] = val;
            from_[key] = file;
        }

}
    static std::string trim(const std::string& str) {
        auto a = str.find_first_not_of(" \t\n\r");
        auto b = str.find_last_not_of(" \t\n\r");
        return (a == std::string::npos) ? "" : str.substr(a, b - a + 1);
    }
};

inline double clamp(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }
inline double sgn(double v) { return (v > 0) - (v < 0); }