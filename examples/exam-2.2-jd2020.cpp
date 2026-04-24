/**
 * Exercise 5: Define JD2020 and convert with different time systems
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include "TimeStruct.h"
#include "TimeConvert.h"
#include "Const.h"

class JD2020 {
public:
    JD2020(long double d = 0.0, TimeSystem ts = TimeSystem::GPS)
        : days(d), timeSystem(ts) {}

    std::string toString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(10) << days << " " << timeSystem.toString();
        return oss.str();
    }

    long double days; // days since 2020-01-01 00:00:00
    TimeSystem timeSystem;
};

inline std::ostream& operator<<(std::ostream& os, const JD2020& jd) {
    os << jd.toString();
    return os;
}

static inline long double jd2020EpochMJD() {
    static const long double epoch = convertYMD2JD(2020, 1, 1) - MJD_TO_JD;
    return epoch;
}

JD2020 CommonTime2JD2020(const CommonTime& ct, const TimeSystem& outTS) {
    CommonTime ctAdj = convertTimeSystem(ct, outTS);

    MJD mjd;
    CommonTime2MJD(ctAdj, mjd);

    long double days = mjd.mjd - jd2020EpochMJD();
    return JD2020(days, outTS);
}

CommonTime JD20202CommonTime(const JD2020& jd, const TimeSystem& outTS) {
    MJD mjd;
    mjd.mjd = jd2020EpochMJD() + jd.days;
    mjd.timeSystem = jd.timeSystem;

    CommonTime ct;
    MJD2CommonTime(mjd, ct);
    return convertTimeSystem(ct, outTS);
}

int main() {
    CivilTime gpsCivil(2025, 1, 4, 9, 0, 0.0, TimeSystem::GPS);
    CivilTime utcCivil(2025, 1, 4, 9, 0, 0.0, TimeSystem::UTC);

    CommonTime gpsCT = CivilTime2CommonTime(gpsCivil);
    CommonTime utcCT = CivilTime2CommonTime(utcCivil);

    JD2020 gpsJD = CommonTime2JD2020(gpsCT, TimeSystem::GPS);
    JD2020 utcJD = CommonTime2JD2020(utcCT, TimeSystem::UTC);

    CommonTime gpsBack = JD20202CommonTime(gpsJD, TimeSystem::GPS);
    CommonTime utcBack = JD20202CommonTime(utcJD, TimeSystem::UTC);

    std::cout << "GPS CivilTime: " << gpsCivil << std::endl;
    std::cout << "GPS JD2020  : " << gpsJD << std::endl;
    std::cout << "GPS back CT: " << gpsBack << std::endl;

    std::cout << "UTC CivilTime: " << utcCivil << std::endl;
    std::cout << "UTC JD2020  : " << utcJD << std::endl;
    std::cout << "UTC back CT: " << utcBack << std::endl;

    return 0;
}
