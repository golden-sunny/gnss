//
// Created by 1 on 2026/3/23.
//

#include "obsreader.h"

#include <sstream>
#include <limits>
#include <cctype>

std::string ObsReader::trim(const std::string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

std::string ObsReader::safeSubstr(const std::string& s, size_t pos, size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

int ObsReader::toInt(const std::string& s, int def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stoi(t); } catch (...) { return def; }
}

double ObsReader::toDouble(const std::string& s, double def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stod(t); } catch (...) { return def; }
}

void ObsReader::padTo80(std::string& line) {
    if (line.size() < 80) line += std::string(80 - line.size(), ' ');
}

ObsEpochData ObsReader::readNextEpoch() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("Obs reader: stream not ready.");
    }
    if (!pHeader) {
        throw std::runtime_error("Obs reader: header not set.");
    }

    std::string line;
    while (true) {
        if (!std::getline(*pStream, line)) {
            throw EndOfObsFile("Obs reader: end of file.");
        }
        if (line.empty()) continue;
        if (line[0] != '>') continue;
        padTo80(line);

        ObsEpochData epoch;
        {
            std::istringstream iss(line.substr(1));
            iss >> epoch.time.year
                >> epoch.time.month
                >> epoch.time.day
                >> epoch.time.hour
                >> epoch.time.minute
                >> epoch.time.second
                >> epoch.epochFlag
                >> epoch.numSat;
            epoch.time.valid = (epoch.time.year > 0);
            if (iss >> epoch.rcvrClockOffset) {
                epoch.hasClockOffset = true;
            }
        }

        // Only epoch flag 0 or 1: observation records follow
        if (epoch.epochFlag != 0 && epoch.epochFlag != 1) {
            // Special event records: skip the following lines
            for (int i = 0; i < epoch.numSat; ++i) {
                if (!std::getline(*pStream, line)) {
                    throw EndOfObsFile("Obs reader: end of file in event block.");
                }
            }
            continue;
        }

        auto looksLikeSatLine = [](const std::string& s) -> bool {
            if (s.size() < 3) return false;
            if (!(std::isalpha(static_cast<unsigned char>(s[0])))) return false;
            if (!std::isdigit(static_cast<unsigned char>(s[1]))) return false;
            if (!std::isdigit(static_cast<unsigned char>(s[2]))) return false;
            return true;
        };

        for (int i = 0; i < epoch.numSat; ++i) {
            if (!std::getline(*pStream, line)) {
                throw EndOfObsFile("Obs reader: end of file in obs block.");
            }
            std::string sat = trim(safeSubstr(line, 0, 3));
            if (sat.empty()) continue;
            std::string sys = sat.substr(0, 1);

            std::vector<std::string> types;
            auto it = pHeader->sysObsTypes.find(sys);
            if (it != pHeader->sysObsTypes.end()) types = it->second;

            int totalTypes = static_cast<int>(types.size());
            int idx = 0;

            // parse first line (after satellite id)
            size_t start = 3;
            while (idx < totalTypes) {
                size_t pos = start;
                while (idx < totalTypes && pos < line.size()) {
                    std::string v = safeSubstr(line, pos, 14);
                    double val = trim(v).empty() ? std::numeric_limits<double>::quiet_NaN()
                                                 : toDouble(v, std::numeric_limits<double>::quiet_NaN());
                    epoch.satTypeValue[sat][types[idx]] = val;
                    idx++;
                    pos += 16;
                }
                if (idx >= totalTypes) break;
                // peek next line: if it looks like a new satellite record, stop and fill remaining as NA
                std::streampos sp = pStream->tellg();
                std::string nextLine;
                if (!std::getline(*pStream, nextLine)) {
                    throw EndOfObsFile("Obs reader: end of file in obs continuation.");
                }
                if (looksLikeSatLine(nextLine)) {
                    pStream->seekg(sp);
                    for (; idx < totalTypes; ++idx) {
                        epoch.satTypeValue[sat][types[idx]] = std::numeric_limits<double>::quiet_NaN();
                    }
                    break;
                }
                line = nextLine;
                start = 0; // continuation lines have no satellite id
            }
        }

        return epoch;
    }
}
