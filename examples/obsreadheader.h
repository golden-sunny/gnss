//
// Created by 1 on 2026/3/20.
//

#ifndef GNSSLAB_OBSREADHEADER_H
#define GNSSLAB_OBSREADHEADER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>

struct RnxTime {
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    double second = 0.0;
    std::string system;
    bool valid = false;
};

struct PhaseShiftRecord {
    std::string system;
    std::string obsCode;
    double correction = 0.0;
    int numSat = -1; // 0 or -1 means all
    std::vector<std::string> satellites;
};

struct ScaleFactorRecord {
    std::string system;
    int factor = 1;      // 1,10,100,1000
    int numObs = -1;     // 0 or -1 means all
    std::vector<std::string> obsTypes;
};

struct SysAppliedRecord {
    std::string system;
    std::string program;
    std::string source;
};

struct GlonassSlotFreqRecord {
    int numSat = 0;
    std::vector<std::pair<std::string, int>> satFreq;
};

struct PrnObsCountRecord {
    std::string sat;
    std::vector<int> counts;
};

struct RinexObsHeaderAll {
    double version = 0.0;
    std::string fileType;
    std::string satSystem;

    std::string pgm;
    std::string runBy;
    std::string date;

    std::vector<std::string> comments;

    std::string markerName;
    std::string markerNumber;
    std::string markerType;

    std::string observer;
    std::string agency;

    std::string recNumber;
    std::string recType;
    std::string recVersion;

    std::string antNumber;
    std::string antType;

    double approxX = 0.0;
    double approxY = 0.0;
    double approxZ = 0.0;

    double antH = 0.0;
    double antE = 0.0;
    double antN = 0.0;

    double antX = 0.0;
    double antY = 0.0;
    double antZ = 0.0;
    bool hasAntDeltaXYZ = false;

    struct AntPhaseCenter {
        std::string system;
        std::string obsCode;
        double d1 = 0.0;
        double d2 = 0.0;
        double d3 = 0.0;
    };
    std::vector<AntPhaseCenter> antPhaseCenters;

    std::vector<std::vector<double>> antBSightXYZ;
    std::vector<double> antZeroDirAzi;
    std::vector<std::vector<double>> antZeroDirXYZ;

    std::vector<std::vector<double>> centerOfMassXYZ;

    std::map<std::string, std::vector<std::string>> sysObsTypes;

    std::string signalStrengthUnit;

    double interval = 0.0;
    bool hasInterval = false;

    RnxTime timeFirstObs;
    RnxTime timeLastObs;

    int rcvClockOffsAppl = -1;

    std::vector<SysAppliedRecord> sysDcbsApplied;
    std::vector<SysAppliedRecord> sysPcvsApplied;

    std::vector<ScaleFactorRecord> sysScaleFactors;
    std::vector<PhaseShiftRecord> sysPhaseShifts;

    GlonassSlotFreqRecord gloSlotFreq;

    int leapSeconds = 0;
    int futureLeapSeconds = 0;
    int weekLeapSeconds = 0;
    int dayLeapSeconds = 0;
    bool hasLeapSeconds = false;

    int numSatellites = 0;
    bool hasNumSatellites = false;

    std::vector<PrnObsCountRecord> prnObsCounts;
};

class ObsReadHeader {
public:
    ObsReadHeader() : pStream(nullptr), headerRead(false) {}

    void setFileStream(std::fstream* stream);
    void parseHeader();
    void printHeader(std::ostream& os = std::cout) const;

    const RinexObsHeaderAll& getHeader() const { return header; }

private:
    std::fstream* pStream;
    RinexObsHeaderAll header;
    bool headerRead;

    static std::string trim(const std::string& s);
    static std::string safeSubstr(const std::string& s, size_t pos, size_t len);
    static int toInt(const std::string& s, int def = 0);
    static double toDouble(const std::string& s, double def = 0.0);

    static std::vector<std::string> splitTokens(const std::string& s);
    static std::vector<std::string> parseObsTypesFromLine(const std::string& line);
    static std::vector<int> parseCountsFromLine(const std::string& line, size_t start);
    static std::vector<std::string> parseSatListByChunk(const std::string& line, size_t start);
    static RnxTime parseRnxTime(const std::string& line);
};

#endif // GNSSLAB_OBSREADHEADER_H
