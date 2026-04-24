//
// Created by 1 on 2026/3/25.
//

#include "navreadheader.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

void NavReadHeader::setFileStream(std::fstream* stream) {
    pStream = stream;
}

std::string NavReadHeader::trim(const std::string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

std::string NavReadHeader::safeSubstr(const std::string& s, size_t pos, size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

int NavReadHeader::toInt(const std::string& s, int def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stoi(t); } catch (...) { return def; }
}

double NavReadHeader::toDouble(const std::string& s, double def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try {
        for (char& c : t) {
            if (c == 'D') c = 'E';
        }
        return std::stod(t);
    } catch (...) {
        return def;
    }
}

std::vector<std::string> NavReadHeader::splitTokens(const std::string& s) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) tokens.push_back(tok);
    return tokens;
}

bool NavReadHeader::looksLikePhaseShiftCont(const std::string& line) {
    std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
    if (tokens.empty()) return false;
    if (tokens.size() == 1) return true;
    return tokens[0].size() == 3 && std::isalpha(static_cast<unsigned char>(tokens[0][0]));
}

void NavReadHeader::parseHeader() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("RINEX nav header: stream not ready.");
    }

    std::string line;
    while (std::getline(*pStream, line)) {
        if (line.size() < 80) line += std::string(80 - line.size(), ' ');
        std::string label = trim(safeSubstr(line, 60, 20));

        if (label == "RINEX VERSION / TYPE") {
            header.version = toDouble(safeSubstr(line, 0, 9), 0.0);
            header.fileType = trim(safeSubstr(line, 20, 1));
            header.satSystem = trim(safeSubstr(line, 40, 1));
        } else if (label == "PGM / RUN BY / DATE") {
            header.pgm = trim(safeSubstr(line, 0, 20));
            header.runBy = trim(safeSubstr(line, 20, 20));
            header.date = trim(safeSubstr(line, 40, 20));
        } else if (label == "COMMENT") {
            header.comments.push_back(trim(safeSubstr(line, 0, 60)));
        } else if (label == "IONOSPHERIC CORR") {
            NavIonoCorrRecord rec;
            std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
            if (!tokens.empty()) {
                rec.type = tokens[0];
                size_t idx = 1;
                while (idx < tokens.size() && rec.coeffs.size() < 4) {
                    if (std::isdigit(static_cast<unsigned char>(tokens[idx][0])) ||
                        tokens[idx][0] == '-' || tokens[idx][0] == '+') {
                        rec.coeffs.push_back(toDouble(tokens[idx], 0.0));
                    }
                    ++idx;
                }
                for (; idx < tokens.size(); ++idx) {
                    if (tokens[idx].size() == 1 &&
                        std::isalpha(static_cast<unsigned char>(tokens[idx][0]))) {
                        rec.timeMark = tokens[idx];
                        rec.hasTimeMark = true;
                    } else if (std::isdigit(static_cast<unsigned char>(tokens[idx][0])) ||
                               tokens[idx][0] == '-' || tokens[idx][0] == '+') {
                        rec.svId = toInt(tokens[idx], -1);
                        rec.hasSvId = true;
                    }
                }
            }
            header.ionoCorrs.push_back(rec);
        } else if (label == "TIME SYSTEM CORR") {
            NavTimeSysCorrRecord rec;
            std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
            if (!tokens.empty()) {
                rec.type = tokens[0];
                if (tokens.size() > 1) rec.a0 = toDouble(tokens[1], 0.0);
                if (tokens.size() > 2) rec.a1 = toDouble(tokens[2], 0.0);
                if (tokens.size() > 3) rec.refSow = toInt(tokens[3], 0);
                if (tokens.size() > 4) rec.refWeek = toInt(tokens[4], 0);
                if (tokens.size() > 5) {
                    rec.source = tokens[5];
                    rec.hasSource = true;
                }
                if (tokens.size() > 6) {
                    rec.utcId = toInt(tokens[6], -1);
                    rec.hasUtcId = true;
                }
            }
            header.timeSysCorrs.push_back(rec);
        } else if (label == "GLONASS COD/PHS/BIS") {
            std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
            for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
                if (tokens[i].empty()) continue;
                if (std::isalpha(static_cast<unsigned char>(tokens[i][0]))) {
                    header.gloCodePhaseBias.codeBias[tokens[i]] = toDouble(tokens[i + 1], 0.0);
                }
            }
        } else if (label == "SYS / PHASE SHIFT" || label == "SYS / PHASE SHIFTS") {
            NavPhaseShiftRecord rec;
            std::vector<std::string> tokens = splitTokens(safeSubstr(line, 0, 60));
            if (!tokens.empty() && tokens[0].size() == 1) {
                rec.system = tokens[0];
                if (tokens.size() > 1) rec.obsCode = tokens[1];
                if (tokens.size() > 2) rec.correction = toDouble(tokens[2], 0.0);
                if (tokens.size() > 3) rec.numSat = toInt(tokens[3], 0);
                for (size_t i = 4; i < tokens.size(); ++i) {
                    rec.satellites.push_back(tokens[i]);
                }
            }
            while (rec.numSat > 0 && static_cast<int>(rec.satellites.size()) < rec.numSat) {
                std::streampos pos = pStream->tellg();
                std::string cont;
                if (!std::getline(*pStream, cont)) break;
                if (cont.size() < 80) cont += std::string(80 - cont.size(), ' ');
                std::string contLabel = trim(safeSubstr(cont, 60, 20));
                if (contLabel != label) {
                    pStream->seekg(pos);
                    break;
                }
                std::vector<std::string> more = splitTokens(safeSubstr(cont, 0, 60));
                for (const auto& sat : more) {
                    if (sat.size() == 3 && std::isalpha(static_cast<unsigned char>(sat[0]))) {
                        rec.satellites.push_back(sat);
                    }
                }
            }
            header.phaseShifts.push_back(rec);
        } else if (label == "LEAP SECONDS") {
            header.leapSeconds = toInt(safeSubstr(line, 0, 6), 0);
            header.futureLeapSeconds = toInt(safeSubstr(line, 6, 6), 0);
            header.leapRefWeek = toInt(safeSubstr(line, 12, 6), 0);
            header.leapRefDay = toInt(safeSubstr(line, 18, 6), 0);
            header.hasLeapSeconds = true;
        } else if (label == "END OF HEADER") {
            headerRead = true;
            break;
        }
    }
}

void NavReadHeader::printHeader(std::ostream& os) const {
    os << "NAV RINEX HEADER\n";
    os << "  Version: " << header.version << "\n";
    os << "  File type: " << header.fileType << "\n";
    os << "  Satellite system: " << header.satSystem << "\n";

    os << "PGM / RUN BY / DATE\n";
    os << "  Program: " << header.pgm << "\n";
    os << "  Run by: " << header.runBy << "\n";
    os << "  Date: " << header.date << "\n";

    os << "COMMENT (" << header.comments.size() << ")\n";
    for (const auto& c : header.comments) {
        os << "  " << c << "\n";
    }

    os << "IONOSPHERIC CORR (" << header.ionoCorrs.size() << ")\n";
    for (const auto& rec : header.ionoCorrs) {
        os << "  " << rec.type << " :";
        for (double c : rec.coeffs) {
            os << " " << c;
        }
        if (rec.hasTimeMark) os << "  timeMark=" << rec.timeMark;
        if (rec.hasSvId) os << "  svId=" << rec.svId;
        os << "\n";
    }

    os << "TIME SYSTEM CORR (" << header.timeSysCorrs.size() << ")\n";
    for (const auto& rec : header.timeSysCorrs) {
        os << "  " << rec.type
           << " a0=" << rec.a0
           << " a1=" << rec.a1
           << " refSow=" << rec.refSow
           << " refWeek=" << rec.refWeek;
        if (rec.hasSource) os << " source=" << rec.source;
        if (rec.hasUtcId) os << " utcId=" << rec.utcId;
        os << "\n";
    }

    os << "GLONASS COD/PHS/BIS (" << header.gloCodePhaseBias.codeBias.size() << ")\n";
    for (const auto& kv : header.gloCodePhaseBias.codeBias) {
        os << "  " << kv.first << " : " << kv.second << "\n";
    }

    os << "SYS / PHASE SHIFT (" << header.phaseShifts.size() << ")\n";
    for (const auto& rec : header.phaseShifts) {
        os << "  " << rec.system << " " << rec.obsCode
           << " corr=" << rec.correction;
        if (rec.numSat <= 0) {
            os << " numSat=ALL";
        } else {
            os << " numSat=" << rec.numSat << " sats=";
            for (const auto& sat : rec.satellites) {
                os << sat << " ";
            }
        }
        os << "\n";
    }

    if (header.hasLeapSeconds) {
        os << "LEAP SECONDS\n";
        os << "  leapSeconds=" << header.leapSeconds
           << " futureLeapSeconds=" << header.futureLeapSeconds
           << " refWeek=" << header.leapRefWeek
           << " refDay=" << header.leapRefDay << "\n";
    }
}
