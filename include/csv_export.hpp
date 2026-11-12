// Defines the class that stores the output into a csv file

#pragma once
#include <iostream>
#include <string>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <vector>

inline std::string outPath(const std::string& dir, const std::string& file) {
    return dir.empty() ? std::string() : dir + "/" + file;
}

class Csv {
private:
    bool first_ = true;
    std::ofstream out_;
    void sep() { 
        if (!first_) 
        out_ << ",";
        first_ = false;
    }

public:
    explicit Csv(const std::string& path) {
        if (path.empty()) return;
        out_.open(path);
        if (!out_)
            throw std::runtime_error("Cannot write " + path);
        out_ << std::setprecision(8);
    }

    Csv& header(const std::vector<std::string>& cols) {
        for (size_t i = 0; i < cols.size(); ++i)
            out_ << ( i ? "," : "" ) << cols[i];
            out_ << "\n";
            return *this;
    }

    // Generic row : accept numbers and strings mixed via the << operator
    template <typename T> Csv& operator<<(const T& v) {
        sep();
        out_ << v;
        return *this;
    }

    // Text fields are quoted when they contain a comma or a quote, otherwise a label might shift the columns
    Csv& operator<<(const std::string& s) {
        sep();
        if (s.find_first_of(",\"") == std::string::npos) {
            out_ << s;
            return *this;
        }
        out_ << '"';
        for (char c : s) {
            if (c == '"')
                out_ << '"';
            out_ << c;
        }
        out_ << '"';
        return *this;
    }

    Csv& operator<<(const char* s) { return *this << std::string(s); }
    Csv& endl() {
        out_ << "\n";
        first_ = true;
        return *this;
    }
};