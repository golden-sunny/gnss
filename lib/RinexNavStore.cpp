/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: Shoujian Zhang锛宻hjzhang@sgg.whu.edu.cn锛?2024-10-10
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hern谩ndez-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */

#include <cmath>
#include <limits>

#include "RinexNavStore.hpp"
#include "StringUtils.h"

using namespace std;
#define debug 0

const string RinexNavStore::stringVersion = "RINEX VERSION / TYPE";
const string RinexNavStore::stringRunBy = "PGM / RUN BY / DATE";
const string RinexNavStore::stringComment = "COMMENT";
const string RinexNavStore::stringIonoCorr = "IONOSPHERIC CORR";
const string RinexNavStore::stringTimeSysCorr = "TIME SYSTEM CORR";
const string RinexNavStore::stringLeapSeconds = "LEAP SECONDS";
//R2.10GLO
const string RinexNavStore::stringCorrSysTime = "CORR TO SYSTEM TIME";
//R2.11GPS
const string RinexNavStore::stringDeltaUTC = "DELTA-UTC: A0,A1,T,W";
//R2.11GEO
const string RinexNavStore::stringDUTC = "D-UTC A0,A1,T,W,S,U";
//R2.11
const string RinexNavStore::stringIonAlpha = "ION ALPHA";
//R2.11
const string RinexNavStore::stringIonBeta = "ION BETA";
const string RinexNavStore::stringEoH = "END OF HEADER";


void RinexNavStore::loadGPSEph(NavEphGPS &gpsEph, string &line, fstream &navFileStream) {

    SatID sat(line.substr(0,3));

    ///add each sat into the satTable
    vector<SatID>::iterator result = find(satTable.begin(), satTable.end(), sat);
    if (result == satTable.end()) {
        satTable.push_back(sat);
    }

    int yr = safeStoi(line.substr(4, 4));
    int mo = safeStoi(line.substr(9, 2));
    int day = safeStoi(line.substr(12, 2));
    int hr = safeStoi(line.substr(15, 2));
    int min = safeStoi(line.substr(18, 2));
    double sec = safeStod(line.substr(21, 2));

    /// Fix RINEX epochs of the form 'yy mm dd hr 59 60.0'
    short ds = 0;
    if (sec >= 60.) {
        ds = sec;
        sec = 0;
    }

    CivilTime cvt(yr, mo, day, hr, min, sec);
    gpsEph.CivilToc = cvt;
//      gpsEph.ctToe = cvt.convertToCommonTime();
    gpsEph.ctToe = CivilTime2CommonTime(cvt);;

    if (ds != 0) gpsEph.ctToe += ds;
    gpsEph.ctToe.setTimeSystem(TimeSystem::GPS);

    GPSWeekSecond gws;
    CommonTime2WeekSecond(gpsEph.ctToe, gws);     // sow is system-independent

    gpsEph.Toc = gws.sow;
    gpsEph.af0 = safeStod(line.substr(23, 19));
    gpsEph.af1 = safeStod(line.substr(42, 19));
    gpsEph.af2 = safeStod(line.substr(61, 19));

    ///orbit-1
    int n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.IODE = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Crs = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Delta_n = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.M0 = safeStod(line.substr(n, 19));
    ///orbit-2
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.Cuc = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.ecc = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Cus = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.sqrt_A = safeStod(line.substr(n, 19));
    ///orbit-3
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.Toe = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Cic = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.OMEGA_0 = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Cis = safeStod(line.substr(n, 19));
    ///orbit-4
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.i0 = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.Crc = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.omega = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.OMEGA_DOT = safeStod(line.substr(n, 19));
    ///orbit-5
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.IDOT = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.L2Codes = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.GPSWeek = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.L2Pflag = safeStod(line.substr(n, 19));
    ///orbit-6
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.URA = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.SV_health = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.TGD = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.IODC = safeStod(line.substr(n, 19));
    ///orbit-7
    n = 4;
    getline(navFileStream, line);
    replace(line.begin(), line.end(), 'D', 'e');
    gpsEph.HOWtime = safeStod(line.substr(n, 19));
    n += 19;
    gpsEph.fitInterval = safeStod(line.substr(n, 19));
    n += 19;

    /// some process
    /// Some RINEX files have HOW < 0.
    while (gpsEph.HOWtime < 0) {
        gpsEph.HOWtime += (long) FULLWEEK;
        gpsEph.GPSWeek--;
    }

    /// In RINEX *files*, weeknum is the week of TOE.
    /// Internally (Rx3NavData), weeknum is week of HOW
    if (gpsEph.HOWtime - gpsEph.Toe > HALFWEEK)
        gpsEph.GPSWeek--;
    else if (gpsEph.HOWtime - gpsEph.Toe < -HALFWEEK)
        gpsEph.GPSWeek++;

    /// Get week for clock, to build Toc
    long adjHOWtime = gpsEph.HOWtime;
    short adjWeeknum = gpsEph.GPSWeek;
    long lToc = (long) gpsEph.Toc;
    if ((gpsEph.HOWtime % SEC_PER_DAY) == 0 &&
        ((lToc) % SEC_PER_DAY) == 0 &&
        gpsEph.HOWtime == lToc) {
        adjHOWtime = gpsEph.HOWtime - 30;
        if (adjHOWtime < 0) {
            adjHOWtime += FULLWEEK;
            adjWeeknum--;
        }
    }

    double dt = gpsEph.Toc - adjHOWtime;
    int week = gpsEph.GPSWeek;
    if (dt < -HALFWEEK) week++; else if (dt > HALFWEEK) week--;
    GPSWeekSecond gws2 = GPSWeekSecond(week, gpsEph.Toc, TimeSystem::GPS);
//      gpsEph.ctToc = GPSWeekSecond(week, gpsEph.Toc, TimeSystem::GPS);
    WeekSecond2CommonTime(gws2, gpsEph.ctToc);

    gpsEph.ctToc.setTimeSystem(TimeSystem::GPS);

    gpsEphData[sat][gpsEph.ctToe] = gpsEph;
}

void RinexNavStore::loadFile(string &file) {
    rx3NavFile = file;
    if (rx3NavFile.size() == 0) {
        cout << "the nav file path is empty!" << endl;
        exit(-1);
    }

    if (debug)
        cout << "RinexNavStore: fileName:" << rx3NavFile << endl;

    fstream navFileStream(rx3NavFile.c_str(), ios::in);
    if (!navFileStream) {
        cerr << "can't open file:" << rx3NavFile << endl;
        exit(-1);
    }

    int lineNumber(0);

    ///first, we should read nav head
    while (1) {
        string line;
        getline(navFileStream, line);

        if (debug)
            cout << "RinexNavStore:" << line << endl;

        stripTrailing(line);

        if (line.length() == 0) continue;
        else if (line.length() < 60) {
            cout << line << endl;
            cout << "line.length is fault" << line.length() << endl;
            FFStreamError e("Invalid line length, \n"
                            "may be the file is generated by windows, \n"
                            "please use dos2unix to convert the file!");
            throw (e);
        }

        lineNumber++;

        string thisLabel(line, 60, 20);

        /// following is huge if else else ... endif for each record type
        if (thisLabel == stringVersion) {
            /// "RINEX VERSION / TYPE"
            version = safeStod(line.substr(0, 20));
            fileType = strip(line.substr(20, 20));
            if(version<3.0)
            {
                FileMissingException e("don't support navigation file with version less than 3.0");
                throw(e);
            }
            if (version >= 3) {                        // ver 3
                if (fileType[0] != 'N' && fileType[0] != 'n') {
                    FFStreamError e("File type is not NAVIGATION: " + fileType);
                    throw(e);
                }
                fileSys = strip(line.substr(40, 20));   // not in ver 2
            }
            fileType = "NAVIGATION";
        } else if (thisLabel == stringRunBy) {
            /// "PGM / RUN BY / DATE"
            fileProgram = strip(line.substr(0, 20));
            fileAgency = strip(line.substr(20, 20));
            // R2 may not have 'UTC' at end
            date = strip(line.substr(40, 20));
        } else if (thisLabel == stringComment) {
            /// "COMMENT"
            commentList.push_back(strip(line.substr(0, 60)));
        } else if (thisLabel == stringIonoCorr) {
            /// "IONOSPHERIC CORR"
            string ionoCorrType = strip(line.substr(0, 4));
            vector<double> ionoCorrCoeff;
            for (int i = 0; i < 4; i++) {
                double ionoCorr = safeStod(line.substr(5 + 12 * i, 12));
                ionoCorrCoeff.push_back(ionoCorr);
            }
            ionoCorrData[ionoCorrType].clear();
            ionoCorrData[ionoCorrType] = ionoCorrCoeff;
        } else if (thisLabel == stringTimeSysCorr) {
            /// "TIME SYSTEM CORR"
            string timeSysCorrType = strip(line.substr(0, 4));

            TimeSysCorr timeSysCorrValue;
            timeSysCorrValue.A0 = safeStod(line.substr(5, 17));
            timeSysCorrValue.A1 = safeStod(line.substr(22, 16));
            timeSysCorrValue.refSOW = safeStoi(line.substr(38, 7));
            timeSysCorrValue.refWeek = safeStoi(line.substr(45, 5));
            timeSysCorrValue.geoProvider = string(" ");
            timeSysCorrValue.geoUTCid = 0;

            timeSysCorrData[timeSysCorrType] = timeSysCorrValue;
        } else if (thisLabel == stringLeapSeconds) {
            /// "LEAP SECONDS"
            leapSeconds = safeStoi(line.substr(0, 6));
            leapDelta = safeStoi(line.substr(6, 6));
            leapWeek = safeStoi(line.substr(12, 6));
            leapDay = safeStoi(line.substr(18, 6));
        } else if (thisLabel == stringEoH) {
            /// "END OF HEADER"
            break;
        }
    }

    ///now, start read nav data
    while (navFileStream.peek() != EOF) {
        string line;
        getline(navFileStream, line);

        if (debug)
            cout << "RinexNavStore:" << line << endl;

        replace(line.begin(), line.end(), 'D', 'e');

        if (line[0] == 'G') {
            NavEphGPS gpsEph;
            loadGPSEph(gpsEph, line, navFileStream);
        }
    }
}

Xvt RinexNavStore::getXvt(const SatID &sat, const CommonTime &epoch) {
    Xvt xvt;
    CommonTime realEpoch;
    TimeSystem ts;
    if (debug)
        cout << sat << endl;
    if (sat.system == "G") {
        ts = TimeSystem::GPS;
        realEpoch = convertTimeSystem(epoch, ts);

        if (debug)
            cout << CommonTime2CivilTime(epoch) << endl;

        NavEphGPS gpsEph = findGPSEph(sat, realEpoch);

        if (debug) {
            cout << "RinexNavStore::GPS eph:" << endl;
            gpsEph.printData();
        }

        xvt = gpsEph.svXvt(realEpoch);

        if (debug) {
            cout << "RinexNavStore::xvt:" << endl;
            cout << xvt << endl;
        }
    } else {
        InvalidRequest e("RinexNavStore: don't support the input satellite system!");
        throw (e);
    }

    return xvt;
}

NavEphGPS RinexNavStore::findGPSEph(const SatID &sat, const CommonTime &epoch) {
    NavEphGPS bestEph;
    const auto satIt = gpsEphData.find(sat);
    if (satIt == gpsEphData.end() || satIt->second.empty()) {
        return bestEph;
    }

    bool found = false;
    double bestAbsDiff = std::numeric_limits<double>::infinity();

    for (const auto &ephEntry : satIt->second) {
        const double diff = std::fabs(epoch - ephEntry.first);
        if (!found || diff < bestAbsDiff) {
            bestAbsDiff = diff;
            bestEph = ephEntry.second;
            found = true;
        }
    }

    return bestEph;
}







