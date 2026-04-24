//
// Created by 1 on 2026/3/25.
//

#include "navreader.h"
#include "BDSWeekSecond.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

void NavReader::setFileStream(std::fstream* stream) {
    pStream = stream;
}

void NavReader::setHeader(const NavHeaderAll* header) {
    pHeader = header;
    (void)pHeader;
}

std::string NavReader::trim(const std::string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

std::string NavReader::safeSubstr(const std::string& s, size_t pos, size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

int NavReader::toInt(const std::string& s, int def) {
    std::string t = trim(s);
    if (t.empty()) return def;
    try { return std::stoi(t); } catch (...) { return def; }
}

double NavReader::toDouble(const std::string& s, double def) {
    std::string t = normalizeD(trim(s));
    if (t.empty()) return def;
    try { return std::stod(t); } catch (...) { return def; }
}

std::string NavReader::normalizeD(const std::string& s) {
    std::string out = s;
    for (char& c : out) {
        if (c == 'D') c = 'E';
    }
    return out;
}

NavGpsRecord NavReader::readNextGpsRecord() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("RINEX nav reader: stream not ready.");
    }

    std::string line;
    while (std::getline(*pStream, line)) {
        if (line.empty()) continue;
        if (line[0] != 'G' && line[0] != 'g') continue;

        if (line.size() < 80) line += std::string(80 - line.size(), ' ');

        NavGpsRecord rec;
        rec.sat = SatID(trim(safeSubstr(line, 0, 3)));

        int yr = toInt(safeSubstr(line, 4, 4), 0);
        int mo = toInt(safeSubstr(line, 9, 2), 0);
        int day = toInt(safeSubstr(line, 12, 2), 0);
        int hr = toInt(safeSubstr(line, 15, 2), 0);
        int min = toInt(safeSubstr(line, 18, 2), 0);
        double sec = toDouble(safeSubstr(line, 21, 2), 0.0);

        short ds = 0;
        if (sec >= 60.) {
            ds = static_cast<short>(sec);
            sec = 0.0;
        }

        CivilTime cvt(yr, mo, day, hr, min, sec);
        rec.eph.CivilToc = cvt;
        rec.eph.ctToc = CivilTime2CommonTime(cvt);
        if (ds != 0) rec.eph.ctToc += ds;
        rec.eph.ctToc.setTimeSystem(TimeSystem::GPS);

        GPSWeekSecond gws;
        CommonTime2WeekSecond(rec.eph.ctToc, gws);
        rec.eph.Toc = gws.sow;

        rec.eph.af0 = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.af1 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.af2 = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 1
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-1.");
        rec.eph.IODE = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Crs = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.Delta_n = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.M0 = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 2
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-2.");
        rec.eph.Cuc = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.ecc = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.Cus = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.sqrt_A = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 3
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-3.");
        rec.eph.Toe = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Cic = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.OMEGA_0 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.Cis = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 4
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-4.");
        rec.eph.i0 = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Crc = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.omega = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.OMEGA_DOT = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 5
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-5.");
        rec.eph.IDOT = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.L2Codes = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.GPSWeek = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.L2Pflag = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 6
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-6.");
        rec.eph.URA = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.SV_health = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.TGD = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.IODC = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 7
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-7.");
        rec.eph.HOWtime = static_cast<long>(toDouble(safeSubstr(line, 4, 19), 0.0));
        rec.eph.fitInterval = toDouble(safeSubstr(line, 23, 19), 0.0);

        while (rec.eph.HOWtime < 0) {
            rec.eph.HOWtime += static_cast<long>(FULLWEEK);
            rec.eph.GPSWeek--;
        }

        if (rec.eph.HOWtime - rec.eph.Toe > HALFWEEK) {
            rec.eph.GPSWeek--;
        } else if (rec.eph.HOWtime - rec.eph.Toe < -HALFWEEK) {
            rec.eph.GPSWeek++;
        }

        long adjHOWtime = rec.eph.HOWtime;
        long lToc = static_cast<long>(rec.eph.Toc);
        if ((rec.eph.HOWtime % SEC_PER_DAY) == 0 &&
            (lToc % SEC_PER_DAY) == 0 &&
            rec.eph.HOWtime == lToc) {
            adjHOWtime = rec.eph.HOWtime - 30;
            if (adjHOWtime < 0) {
                adjHOWtime += FULLWEEK;
            }
        }

        double dt = rec.eph.Toc - adjHOWtime;
        int week = static_cast<int>(rec.eph.GPSWeek);
        if (dt < -HALFWEEK) week++;
        else if (dt > HALFWEEK) week--;

        // Toe is a broadcast ephemeris reference epoch, so it should be
        // converted from GPS week + seconds-of-week rather than set manually.
        GPSWeekSecond toeGws(week, rec.eph.Toe, TimeSystem::GPS);
        WeekSecond2CommonTime(toeGws, rec.eph.ctToe);
        rec.eph.ctToe.setTimeSystem(TimeSystem::GPS);

        GPSWeekSecond gws2 = GPSWeekSecond(week, rec.eph.Toc, TimeSystem::GPS);
        WeekSecond2CommonTime(gws2, rec.eph.ctToc);
        rec.eph.ctToc.setTimeSystem(TimeSystem::GPS);

        gpsEphData[rec.sat][rec.eph.ctToe] = rec.eph;
        return rec;
    }

    throw EndOfNavFile("RINEX nav reader: end of file.");
}

NavBdsRecord NavReader::readNextBdsRecord() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("RINEX nav reader: stream not ready.");
    }

    std::string line;
    while (std::getline(*pStream, line)) {
        if (line.empty()) continue;
        if (line[0] != 'C' && line[0] != 'c') continue;

        if (line.size() < 80) line += std::string(80 - line.size(), ' ');

        NavBdsRecord rec;
        rec.sat = SatID(trim(safeSubstr(line, 0, 3)));

        int yr = toInt(safeSubstr(line, 4, 4), 0);
        int mo = toInt(safeSubstr(line, 9, 2), 0);
        int day = toInt(safeSubstr(line, 12, 2), 0);
        int hr = toInt(safeSubstr(line, 15, 2), 0);
        int min = toInt(safeSubstr(line, 18, 2), 0);
        double sec = toDouble(safeSubstr(line, 21, 2), 0.0);

        short ds = 0;
        if (sec >= 60.) {
            ds = static_cast<short>(sec);
            sec = 0.0;
        }

        CivilTime cvt(yr, mo, day, hr, min, sec, TimeSystem::BDT);
        rec.eph.CivilToc = cvt;
        rec.eph.ctToc = CivilTime2CommonTime(cvt);
        if (ds != 0) rec.eph.ctToc += ds;
        rec.eph.ctToc.setTimeSystem(TimeSystem::BDT);

        BDSWeekSecond bws;
        CommonTime2WeekSecond(rec.eph.ctToc, bws);
        rec.eph.Toc = bws.sow;

        rec.eph.af0 = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.af1 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.af2 = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 1
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-1.");
        rec.eph.IODE = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Crs = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.Delta_n = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.M0 = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 2
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-2.");
        rec.eph.Cuc = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.ecc = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.Cus = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.sqrt_A = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 3
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-3.");
        rec.eph.Toe = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Cic = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.OMEGA_0 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.Cis = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 4
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-4.");
        rec.eph.i0 = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.Crc = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.omega = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.OMEGA_DOT = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 5
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-5.");
        rec.eph.IDOT = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.L2Codes = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.BDTWeek = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.L2Pflag = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 6
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-6.");
        rec.eph.URA = toDouble(safeSubstr(line, 4, 19), 0.0);
        rec.eph.SV_health = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.TGD1 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.TGD2 = toDouble(safeSubstr(line, 61, 19), 0.0);

        // BROADCAST ORBIT - 7
        if (!std::getline(*pStream, line)) throw EndOfNavFile("RINEX nav reader: EOF in orbit-7.");
        rec.eph.HOWtime = static_cast<long>(toDouble(safeSubstr(line, 4, 19), 0.0));
        rec.eph.fitInterval = toDouble(safeSubstr(line, 23, 19), 0.0);
        rec.eph.spare1 = toDouble(safeSubstr(line, 42, 19), 0.0);
        rec.eph.spare2 = toDouble(safeSubstr(line, 61, 19), 0.0);

        while (rec.eph.HOWtime < 0) {
            rec.eph.HOWtime += static_cast<long>(FULLWEEK);
            rec.eph.BDTWeek--;
        }

        if (rec.eph.HOWtime - rec.eph.Toe > HALFWEEK) {
            rec.eph.BDTWeek--;
        } else if (rec.eph.HOWtime - rec.eph.Toe < -HALFWEEK) {
            rec.eph.BDTWeek++;
        }

        long adjHOWtime = rec.eph.HOWtime;
        long lToc = static_cast<long>(rec.eph.Toc);
        if ((rec.eph.HOWtime % SEC_PER_DAY) == 0 &&
            (lToc % SEC_PER_DAY) == 0 &&
            rec.eph.HOWtime == lToc) {
            adjHOWtime = rec.eph.HOWtime - 30;
            if (adjHOWtime < 0) {
                adjHOWtime += FULLWEEK;
            }
        }

        double dt = rec.eph.Toc - adjHOWtime;
        int week = static_cast<int>(rec.eph.BDTWeek);
        if (dt < -HALFWEEK) week++;
        else if (dt > HALFWEEK) week--;

        BDSWeekSecond toeBws(week, rec.eph.Toe, TimeSystem::BDT);
        WeekSecond2CommonTime(toeBws, rec.eph.ctToe);
        rec.eph.ctToe.setTimeSystem(TimeSystem::BDT);

        BDSWeekSecond tocBws(week, rec.eph.Toc, TimeSystem::BDT);
        WeekSecond2CommonTime(tocBws, rec.eph.ctToc);
        rec.eph.ctToc.setTimeSystem(TimeSystem::BDT);

        return rec;
    }

    throw EndOfNavFile("RINEX nav reader: end of file.");
}

std::set<SatID> NavReader::collectUniqueSats() {
    if (!pStream || !(*pStream)) {
        throw std::runtime_error("RINEX nav reader: stream not ready.");
    }

    std::set<SatID> uniqueSats;
    pStream->clear();
    pStream->seekg(0, std::ios::beg);

    std::string line;
    bool inHeader = true;
    while (std::getline(*pStream, line)) {
        if (inHeader) {
            std::string label;
            if (line.size() >= 80) {
                label = trim(safeSubstr(line, 60, 20));
            } else {
                line += std::string(80 - line.size(), ' ');
                label = trim(safeSubstr(line, 60, 20));
            }
            if (label == "END OF HEADER") {
                inHeader = false;
            }
            continue;
        }

        if (line.size() < 3) continue;
        if (!std::isalpha(static_cast<unsigned char>(line[0]))) continue;
        if (!std::isdigit(static_cast<unsigned char>(line[1])) ||
            !std::isdigit(static_cast<unsigned char>(line[2]))) {
            continue;
        }

        uniqueSats.insert(SatID(trim(safeSubstr(line, 0, 3))));
    }

    pStream->clear();
    pStream->seekg(0, std::ios::beg);
    return uniqueSats;
}
