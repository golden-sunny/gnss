//
// Created by 1 on 2026/3/20.
//

#include "obsreadheader.h"

#include <sstream>
#include <iomanip>
#include <stdexcept>

void ObsReadHeader::setFileStream(std::fstream* stream) {
    pStream = stream;
}

std::string ObsReadHeader::trim(const std::string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

std::string ObsReadHeader::safeSubstr(const std::string& s, size_t pos, size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

int ObsReadHeader::toInt(const std::string& s, int def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stoi(t); } catch (...) { return def; }
}

double ObsReadHeader::toDouble(const std::string& s, double def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stod(t); } catch (...) { return def; }
}

std::vector<std::string> ObsReadHeader::splitTokens(const std::string& s) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string t;
    while (iss >> t) tokens.push_back(t);
    return tokens;
}

std::vector<std::string> ObsReadHeader::parseObsTypesFromLine(const std::string& line) {
    std::vector<std::string> types;
    for (size_t i = 7; i + 3 <= 60; i += 4) {
        std::string t = trim(safeSubstr(line, i, 3));
        if (!t.empty()) types.push_back(t);
    }
    return types;
}

std::vector<int> ObsReadHeader::parseCountsFromLine(const std::string& line, size_t start) {
    std::vector<int> counts;
    for (size_t i = start; i + 6 <= 60; i += 6) {
        std::string v = trim(safeSubstr(line, i, 6));
        if (!v.empty()) counts.push_back(toInt(v, 0));
    }
    return counts;
}

std::vector<std::string> ObsReadHeader::parseSatListByChunk(const std::string& line, size_t start) {
    std::vector<std::string> sats;
    for (size_t i = start; i + 3 <= 60; i += 4) {
        std::string t = trim(safeSubstr(line, i, 3));
        if (!t.empty()) sats.push_back(t);
    }
    return sats;
}

RnxTime ObsReadHeader::parseRnxTime(const std::string& line) {
    RnxTime t;
    t.year = toInt(safeSubstr(line, 0, 6), 0);
    t.month = toInt(safeSubstr(line, 6, 6), 0);
    t.day = toInt(safeSubstr(line, 12, 6), 0);
    t.hour = toInt(safeSubstr(line, 18, 6), 0);
    t.minute = toInt(safeSubstr(line, 24, 6), 0);
    t.second = toDouble(safeSubstr(line, 30, 13), 0.0);
    t.system = trim(safeSubstr(line, 48, 3));
    t.valid = (t.year > 0);
    return t;
}

void ObsReadHeader::parseHeader() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("RINEX obs header: stream not ready.");
    }

    std::string line;
    while (std::getline(*pStream, line)) {
        if (line.size() < 80) line += std::string(80 - line.size(), ' ');
        std::string label = trim(safeSubstr(line, 60, 20));

        if (label == "RINEX VERSION / TYPE") {
            header.version = toDouble(safeSubstr(line, 0, 9), 0.0);
            header.fileType =  trim(safeSubstr(line, 20, 20));
            header.satSystem = trim(safeSubstr(line, 40, 1));
        } else if (label == "PGM / RUN BY / DATE") {
            header.pgm = trim(safeSubstr(line, 0, 20));
            header.runBy = trim(safeSubstr(line, 20, 20));
            header.date = trim(safeSubstr(line, 40, 20));
        } else if (label == "COMMENT") {
            header.comments.push_back(trim(safeSubstr(line, 0, 60)));
        } else if (label == "MARKER NAME") {
            header.markerName = trim(safeSubstr(line, 0, 60));
        } else if (label == "MARKER NUMBER") {
            header.markerNumber = trim(safeSubstr(line, 0, 20));
        } else if (label == "MARKER TYPE") {
            header.markerType = trim(safeSubstr(line, 0, 20));
        } else if (label == "OBSERVER / AGENCY") {
            header.observer = trim(safeSubstr(line, 0, 20));
            header.agency = trim(safeSubstr(line, 20, 40));
        } else if (label == "REC # / TYPE / VERS") {
            header.recNumber = trim(safeSubstr(line, 0, 20));
            header.recType = trim(safeSubstr(line, 20, 20));
            header.recVersion = trim(safeSubstr(line, 40, 20));
        } else if (label == "ANT # / TYPE") {
            header.antNumber = trim(safeSubstr(line, 0, 20));
            header.antType = trim(safeSubstr(line, 20, 20));
        } else if (label == "APPROX POSITION XYZ") {
            header.approxX = toDouble(safeSubstr(line, 0, 14), 0.0);
            header.approxY = toDouble(safeSubstr(line, 14, 14), 0.0);
            header.approxZ = toDouble(safeSubstr(line, 28, 14), 0.0);
        } else if (label == "ANTENNA: DELTA H/E/N") {
            header.antH = toDouble(safeSubstr(line, 0, 14), 0.0);
            header.antE = toDouble(safeSubstr(line, 14, 14), 0.0);
            header.antN = toDouble(safeSubstr(line, 28, 14), 0.0);
        } else if (label == "ANTENNA: DELTA X/Y/Z") {
            header.antX = toDouble(safeSubstr(line, 0, 14), 0.0);
            header.antY = toDouble(safeSubstr(line, 14, 14), 0.0);
            header.antZ = toDouble(safeSubstr(line, 28, 14), 0.0);
            header.hasAntDeltaXYZ = true;
        } else if (label == "ANTENNA: PHASECENTER") {
            RinexObsHeaderAll::AntPhaseCenter apc;
            apc.system = trim(safeSubstr(line, 0, 1));
            apc.obsCode = trim(safeSubstr(line, 2, 3));
            apc.d1 = toDouble(safeSubstr(line, 5, 9), 0.0);
            apc.d2 = toDouble(safeSubstr(line, 14, 9), 0.0);
            apc.d3 = toDouble(safeSubstr(line, 23, 9), 0.0);
            header.antPhaseCenters.push_back(apc);
        } else if (label == "ANTENNA: B.SIGHT XYZ") {
            std::vector<double> v(3, 0.0);
            v[0] = toDouble(safeSubstr(line, 0, 14), 0.0);
            v[1] = toDouble(safeSubstr(line, 14, 14), 0.0);
            v[2] = toDouble(safeSubstr(line, 28, 14), 0.0);
            header.antBSightXYZ.push_back(v);
        } else if (label == "ANTENNA: ZERODIR AZI") {
            header.antZeroDirAzi.push_back(toDouble(safeSubstr(line, 0, 14), 0.0));
        } else if (label == "ANTENNA: ZERODIR XYZ") {
            std::vector<double> v(3, 0.0);
            v[0] = toDouble(safeSubstr(line, 0, 14), 0.0);
            v[1] = toDouble(safeSubstr(line, 14, 14), 0.0);
            v[2] = toDouble(safeSubstr(line, 28, 14), 0.0);
            header.antZeroDirXYZ.push_back(v);
        } else if (label == "CENTER OF MASS: XYZ") {
            std::vector<double> v(3, 0.0);
            v[0] = toDouble(safeSubstr(line, 0, 14), 0.0);
            v[1] = toDouble(safeSubstr(line, 14, 14), 0.0);
            v[2] = toDouble(safeSubstr(line, 28, 14), 0.0);
            header.centerOfMassXYZ.push_back(v);
        } else if (label == "SYS / # / OBS TYPES") {
            std::string sys = trim(safeSubstr(line, 0, 1));
            int num = toInt(safeSubstr(line, 3, 3), 0);
            std::vector<std::string> types = parseObsTypesFromLine(line);
            while ((int)types.size() < num) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != "SYS / # / OBS TYPES") {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<std::string> more = parseObsTypesFromLine(cont);
                types.insert(types.end(), more.begin(), more.end());
            }
            if (!sys.empty()) header.sysObsTypes[sys] = types;
        } else if (label == "SIGNAL STRENGTH UNIT") {
            header.signalStrengthUnit = trim(safeSubstr(line, 0, 20));
        } else if (label == "INTERVAL") {
            header.interval = toDouble(safeSubstr(line, 0, 10), 0.0);
            header.hasInterval = true;
        } else if (label == "TIME OF FIRST OBS") {
            header.timeFirstObs = parseRnxTime(line);
        } else if (label == "TIME OF LAST OBS") {
            header.timeLastObs = parseRnxTime(line);
        } else if (label == "RCV CLOCK OFFS APPL") {
            header.rcvClockOffsAppl = toInt(safeSubstr(line, 0, 6), -1);
        } else if (label == "SYS / DCBS APPLIED") {
            SysAppliedRecord rec;
            rec.system = trim(safeSubstr(line, 0, 1));
            rec.program = trim(safeSubstr(line, 2, 17));
            rec.source = trim(safeSubstr(line, 20, 40));
            header.sysDcbsApplied.push_back(rec);
        } else if (label == "SYS / PCVS APPLIED") {
            SysAppliedRecord rec;
            rec.system = trim(safeSubstr(line, 0, 1));
            rec.program = trim(safeSubstr(line, 2, 17));
            rec.source = trim(safeSubstr(line, 20, 40));
            header.sysPcvsApplied.push_back(rec);
        } else if (label == "SYS / SCALE FACTOR") {
            ScaleFactorRecord rec;
            rec.system = trim(safeSubstr(line, 0, 1));
            rec.factor = toInt(safeSubstr(line, 2, 4), 1);
            rec.numObs = toInt(safeSubstr(line, 8, 2), -1);
            rec.obsTypes = parseObsTypesFromLine(line);
            while (rec.numObs > 0 && (int)rec.obsTypes.size() < rec.numObs) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != "SYS / SCALE FACTOR") {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<std::string> more = parseObsTypesFromLine(cont);
                rec.obsTypes.insert(rec.obsTypes.end(), more.begin(), more.end());
            }
            header.sysScaleFactors.push_back(rec);
        } else if (label == "SYS / PHASE SHIFTS") {
            PhaseShiftRecord rec;
            rec.system = trim(safeSubstr(line, 0, 1));
            rec.obsCode = trim(safeSubstr(line, 2, 3));
            rec.correction = toDouble(safeSubstr(line, 6, 8), 0.0);
            // 2X, I2.2: blank or 0 means all satellites of the system
            std::string numField = safeSubstr(line, 16, 2);
            rec.numSat = trim(numField).empty() ? 0 : toInt(numField, 0);
            std::vector<std::string> sats = parseSatListByChunk(line, 20);
            rec.satellites.insert(rec.satellites.end(), sats.begin(), sats.end());
            while (rec.numSat > 0 && (int)rec.satellites.size() < rec.numSat) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != "SYS / PHASE SHIFTS") {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<std::string> more = parseSatListByChunk(cont, 18);
                rec.satellites.insert(rec.satellites.end(), more.begin(), more.end());
            }
            header.sysPhaseShifts.push_back(rec);
        } else if (label == "GLONASS SLOT / FRQ #") {
            if (header.gloSlotFreq.numSat == 0) {
                header.gloSlotFreq.numSat = toInt(safeSubstr(line, 0, 3), 0);
            }
            std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
            size_t start = 0;
            if (!tokens.empty() && tokens[0].size() <= 3 && std::isdigit(tokens[0][0])) {
                start = 1;
            }
            for (size_t i = start; i + 1 < tokens.size(); i += 2) {
                std::string sat = tokens[i];
                int frq = toInt(tokens[i + 1], 0);
                if (!sat.empty()) header.gloSlotFreq.satFreq.emplace_back(sat, frq);
            }
            while (header.gloSlotFreq.numSat > 0 &&
                   (int)header.gloSlotFreq.satFreq.size() < header.gloSlotFreq.numSat) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != "GLONASS SLOT / FRQ #") {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<std::string> moreTokens = splitTokens(safeSubstr(cont, 0, 60));
                for (size_t i = 0; i + 1 < moreTokens.size(); i += 2) {
                    std::string sat = moreTokens[i];
                    int frq = toInt(moreTokens[i + 1], 0);
                    if (!sat.empty()) header.gloSlotFreq.satFreq.emplace_back(sat, frq);
                }
            }
        } else if (label == "LEAP SECONDS") {
            header.leapSeconds = toInt(safeSubstr(line, 0, 6), 0);
            header.futureLeapSeconds = toInt(safeSubstr(line, 6, 6), 0);
            header.weekLeapSeconds = toInt(safeSubstr(line, 12, 6), 0);
            header.dayLeapSeconds = toInt(safeSubstr(line, 18, 6), 0);
            header.hasLeapSeconds = true;
        } else if (label == "# OF SATELLITES") {
            header.numSatellites = toInt(safeSubstr(line, 0, 6), 0);
            header.hasNumSatellites = true;
        } else if (label == "PRN / # OF OBS") {
            PrnObsCountRecord rec;
            rec.sat = trim(safeSubstr(line, 0, 3));
            rec.counts = parseCountsFromLine(line, 3);
            std::string sys = rec.sat.empty() ? "" : rec.sat.substr(0, 1);
            int expected = 0;
            if (!sys.empty() && header.sysObsTypes.count(sys)) {
                expected = (int)header.sysObsTypes[sys].size();
            }
            while (expected > 0 && (int)rec.counts.size() < expected) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != "PRN / # OF OBS") {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<int> more = parseCountsFromLine(cont, 6);
                rec.counts.insert(rec.counts.end(), more.begin(), more.end());
            }
            header.prnObsCounts.push_back(rec);
        } else if (label == "END OF HEADER") {
            headerRead = true;
            break;
        }
    }
}

void ObsReadHeader::printHeader(std::ostream& os) const {
    os << "RINEX VERSION / TYPE\n";
    os << "  Version: " << header.version << "\n";
    os << "  File type: " << header.fileType << "\n";
    os << "  Satellite system: " << header.satSystem << "\n";

    os << "PGM / RUN BY / DATE\n";
    os << "  Program: " << header.pgm << "\n";
    os << "  Run by: " << header.runBy << "\n";
    os << "  Date: " << header.date << "\n";

    os << "COMMENT (" << header.comments.size() << ")\n";
    for (const auto& c : header.comments) os << "  " << c << "\n";

    os << "MARKER\n";
    os << "  Name: " << header.markerName << "\n";
    os << "  Number: " << header.markerNumber << "\n";
    os << "  Type: " << header.markerType << "\n";

    os << "OBSERVER / AGENCY\n";
    os << "  Observer: " << header.observer << "\n";
    os << "  Agency: " << header.agency << "\n";

    os << "REC # / TYPE / VERS\n";
    os << "  Number: " << header.recNumber << "\n";
    os << "  Type: " << header.recType << "\n";
    os << "  Version: " << header.recVersion << "\n";

    os << "ANT # / TYPE\n";
    os << "  Number: " << header.antNumber << "\n";
    os << "  Type: " << header.antType << "\n";

    os << "APPROX POSITION XYZ\n";
    os << "  X: " << header.approxX << "  Y: " << header.approxY << "  Z: " << header.approxZ << "\n";

    os << "ANTENNA: DELTA H/E/N\n";
    os << "  H: " << header.antH << "  E: " << header.antE << "  N: " << header.antN << "\n";

    if (header.hasAntDeltaXYZ) {
        os << "ANTENNA: DELTA X/Y/Z\n";
        os << "  X: " << header.antX << "  Y: " << header.antY << "  Z: " << header.antZ << "\n";
    }

    os << "ANTENNA: PHASECENTER (" << header.antPhaseCenters.size() << ")\n";
    for (const auto& apc : header.antPhaseCenters) {
        os << "  " << apc.system << " " << apc.obsCode
           << "  " << apc.d1 << " " << apc.d2 << " " << apc.d3 << "\n";
    }

    os << "ANTENNA: B.SIGHT XYZ (" << header.antBSightXYZ.size() << ")\n";
    for (const auto& v : header.antBSightXYZ) {
        os << "  " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    os << "ANTENNA: ZERODIR AZI (" << header.antZeroDirAzi.size() << ")\n";
    for (double v : header.antZeroDirAzi) os << "  " << v << "\n";

    os << "ANTENNA: ZERODIR XYZ (" << header.antZeroDirXYZ.size() << ")\n";
    for (const auto& v : header.antZeroDirXYZ) {
        os << "  " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    os << "CENTER OF MASS: XYZ (" << header.centerOfMassXYZ.size() << ")\n";
    for (const auto& v : header.centerOfMassXYZ) {
        os << "  " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    os << "SYS / # / OBS TYPES\n";
    for (const auto& kv : header.sysObsTypes) {
        os << "  " << kv.first << " (" << kv.second.size() << "): ";
        for (const auto& t : kv.second) os << t << " ";
        os << "\n";
    }

    if (!header.signalStrengthUnit.empty()) {
        os << "SIGNAL STRENGTH UNIT: " << header.signalStrengthUnit << "\n";
    }

    if (header.hasInterval) {
        os << "INTERVAL: " << header.interval << " s\n";
    }

    if (header.timeFirstObs.valid) {
        os << "TIME OF FIRST OBS: " << header.timeFirstObs.year << " "
           << header.timeFirstObs.month << " " << header.timeFirstObs.day << " "
           << header.timeFirstObs.hour << " " << header.timeFirstObs.minute << " "
           << std::fixed << std::setprecision(7) << header.timeFirstObs.second
           << " " << header.timeFirstObs.system << "\n";
    }
    if (header.timeLastObs.valid) {
        os << "TIME OF LAST OBS: " << header.timeLastObs.year << " "
           << header.timeLastObs.month << " " << header.timeLastObs.day << " "
           << header.timeLastObs.hour << " " << header.timeLastObs.minute << " "
           << std::fixed << std::setprecision(7) << header.timeLastObs.second
           << " " << header.timeLastObs.system << "\n";
    }

    if (header.rcvClockOffsAppl >= 0) {
        os << "RCV CLOCK OFFS APPL: " << header.rcvClockOffsAppl << "\n";
    }

    os << "SYS / DCBS APPLIED (" << header.sysDcbsApplied.size() << ")\n";
    for (const auto& r : header.sysDcbsApplied) {
        os << "  " << r.system << " " << r.program << " " << r.source << "\n";
    }

    os << "SYS / PCVS APPLIED (" << header.sysPcvsApplied.size() << ")\n";
    for (const auto& r : header.sysPcvsApplied) {
        os << "  " << r.system << " " << r.program << " " << r.source << "\n";
    }

    os << "SYS / SCALE FACTOR (" << header.sysScaleFactors.size() << ")\n";
    for (const auto& r : header.sysScaleFactors) {
        os << "  " << r.system << " factor=" << r.factor << " num=" << r.numObs << " types=";
        for (const auto& t : r.obsTypes) os << t << " ";
        os << "\n";
    }

    os << "SYS / PHASE SHIFTS (" << header.sysPhaseShifts.size() << ")\n";
    for (const auto& r : header.sysPhaseShifts) {
        os << "  " << r.system << " " << r.obsCode << " corr=" << r.correction;
        if (r.numSat <= 0) {
            os << " numSat=ALL";
        } else {
            os << " numSat=" << r.numSat << " sats=";
            for (const auto& s : r.satellites) os << s << " ";
        }
        os << "\n";
    }

    if (header.gloSlotFreq.numSat > 0) {
        os << "GLONASS SLOT / FRQ # (" << header.gloSlotFreq.numSat << ")\n";
        for (const auto& p : header.gloSlotFreq.satFreq) {
            os << "  " << p.first << " " << p.second << "\n";
        }
    }

    if (header.hasLeapSeconds) {
        os << "LEAP SECONDS: " << header.leapSeconds
           << "  future=" << header.futureLeapSeconds
           << "  week=" << header.weekLeapSeconds
           << "  day=" << header.dayLeapSeconds << "\n";
    }

    if (header.hasNumSatellites) {
        os << "# OF SATELLITES: " << header.numSatellites << "\n";
    }

    os << "PRN / # OF OBS (" << header.prnObsCounts.size() << ")\n";
    for (const auto& r : header.prnObsCounts) {
        os << "  " << r.sat << " : ";
        for (int v : r.counts) os << v << " ";
        os << "\n";
    }
}
