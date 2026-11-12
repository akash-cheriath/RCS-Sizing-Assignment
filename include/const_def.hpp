/*
Defines the constants used in the RCS sizing assignment. These constants are used throughout the project to 
ensure consistency and accuracy in calculations.
*/

#pragma once

namespace cnst {
    constexpr double PI                     = 3.14159265358979323846;       // Pi
    constexpr double DEG                    = PI / 180.0;                   // Degrees to radians conversion factor
    constexpr double MU_EARTH               = 3.986004418e14;               // Standard gravitational parameter for Earth [m^3 s^-2]
    constexpr double G0                     = 9.80665;                      // Standard gravity [m/s^2]
    constexpr double RADIUS_EARTH           = 6378137;                      // Equatorial radius of Earth [m]
}

