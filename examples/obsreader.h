//
// Created by 1 on 2026/3/23.
//

#ifndef GNSSLAB_OBSREADER_H
#define GNSSLAB_OBSREADER_H

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

#include "obsreadheader.h"

struct ObsEpochData {
    RnxTime time;
    int epochFlag = 0;
    int numSat = 0;
    double rcvrClockOffset = 0.0;
    bool hasClockOffset = false;

    // sat -> obsType -> value
    std::map<std::string, std::map<std::string, double>> satTypeValue;
};

class EndOfObsFile : public std::runtime_error {
public:
    explicit EndOfObsFile(const std::string& msg) : std::runtime_error(msg) {}
};

class ObsReader {
public:
    ObsReader() : pStream(nullptr), pHeader(nullptr) {}

    void setFileStream(std::fstream* stream) { pStream = stream; }
    void setHeader(const RinexObsHeaderAll* header) { pHeader = header; }

    // Read next epoch observation block; throw EndOfObsFile on EOF
    ObsEpochData readNextEpoch();

private:
    std::fstream* pStream;
    const RinexObsHeaderAll* pHeader;

    static std::string trim(const std::string& s);
    static std::string safeSubstr(const std::string& s, size_t pos, size_t len);
    static int toInt(const std::string& s, int def = 0);
    static double toDouble(const std::string& s, double def = 0.0);

    static void padTo80(std::string& line);
};

#endif // GNSSLAB_OBSREADER_H
