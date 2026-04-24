//
// Created by 1 on 2026/3/25.
//

#ifndef GNSSLAB_NAVREADHEADER_H
#define GNSSLAB_NAVREADHEADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct NavIonoCorrRecord {
    std::string type;
    std::vector<double> coeffs;
    std::string timeMark;
    int svId = -1;
    bool hasTimeMark = false;
    bool hasSvId = false;
};

struct NavTimeSysCorrRecord {
    std::string type;
    double a0 = 0.0;
    double a1 = 0.0;
    int refSow = 0;
    int refWeek = 0;
    std::string source;
    int utcId = -1;
    bool hasSource = false;
    bool hasUtcId = false;
};

struct NavGloCodePhaseBiasRecord {
    std::map<std::string, double> codeBias;
};

struct NavPhaseShiftRecord {
    std::string system;
    std::string obsCode;
    double correction = 0.0;
    int numSat = -1; // 0 or blank means all satellites of the system
    std::vector<std::string> satellites;
};

struct NavHeaderAll {
    double version = 0.0;
    std::string fileType;
    std::string satSystem;

    std::string pgm;
    std::string runBy;
    std::string date;

    std::vector<std::string> comments;

    std::vector<NavIonoCorrRecord> ionoCorrs;
    std::vector<NavTimeSysCorrRecord> timeSysCorrs;
    NavGloCodePhaseBiasRecord gloCodePhaseBias;
    std::vector<NavPhaseShiftRecord> phaseShifts;

    int leapSeconds = 0;
    int futureLeapSeconds = 0;
    int leapRefWeek = 0;
    int leapRefDay = 0;
    bool hasLeapSeconds = false;
};

class NavReadHeader {
public:
    NavReadHeader() : pStream(nullptr), headerRead(false) {}

    void setFileStream(std::fstream* stream);
    void parseHeader();
    void printHeader(std::ostream& os = std::cout) const;

    const NavHeaderAll& getHeader() const { return header; }

private:
    std::fstream* pStream;
    NavHeaderAll header;
    bool headerRead;

    static std::string trim(const std::string& s);
    static std::string safeSubstr(const std::string& s, size_t pos, size_t len);
    static int toInt(const std::string& s, int def = 0);
    static double toDouble(const std::string& s, double def = 0.0);
    static std::vector<std::string> splitTokens(const std::string& s);
    static bool looksLikePhaseShiftCont(const std::string& line);
};

#endif // GNSSLAB_NAVREADHEADER_H
