//
// Created by 1 on 2026/3/25.
//

#ifndef GNSSLAB_NAVREADER_H
#define GNSSLAB_NAVREADER_H

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <stdexcept>

#include "GnssStruct.h"
#include "NavEphBDS.hpp"
#include "NavEphGPS.hpp"
#include "navreadheader.h"

struct NavGpsRecord {
    SatID sat;
    NavEphGPS eph;
};

struct NavBdsRecord {
    SatID sat;
    NavEphBDS eph;
};

enum class NavRecordSystem {
    Unknown = 0,
    GPS,
    BDS
};

class EndOfNavFile : public std::runtime_error {
public:
    explicit EndOfNavFile(const std::string& msg) : std::runtime_error(msg) {}
};

class NavReader {
public:
    NavReader() : pStream(nullptr), pHeader(nullptr) {}

    void setFileStream(std::fstream* stream);
    void setHeader(const NavHeaderAll* header);

    NavGpsRecord readNextGpsRecord();
    NavBdsRecord readNextBdsRecord();
    const std::map<SatID, std::map<CommonTime, NavEphGPS>>& getGpsEphData() const { return gpsEphData; }
    std::set<SatID> collectUniqueSats();

private:
    std::fstream* pStream;
    const NavHeaderAll* pHeader;
    std::map<SatID, std::map<CommonTime, NavEphGPS>> gpsEphData;

    static std::string trim(const std::string& s);
    static std::string safeSubstr(const std::string& s, size_t pos, size_t len);
    static int toInt(const std::string& s, int def = 0);
    static double toDouble(const std::string& s, double def = 0.0);
    static std::string normalizeD(const std::string& s);
};

#endif // GNSSLAB_NAVREADER_H
