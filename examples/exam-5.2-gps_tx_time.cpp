/**
 * Exercise 5.1 and 5.2 example.
 *
 * This program compares GPS satellite transmit time computed from:
 * - raw C1 pseudorange
 * - IF-combined pseudorange
 *
 * It also outputs the GPS TGD correction term so the variation can be plotted
 * later from the generated CSV file.
 */

#include <cmath>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <limits>

#include "Const.h"
#include "GnssFunc.h"
#include "NavEphGPS.hpp"
#include "RinexNavStore.hpp"
#include "RinexObsReader.h"
#include "TimeConvert.h"

using namespace std;

namespace {

struct TxStats {
    size_t n{0};
    double sum{0.0};
    double sumSq{0.0};
    double minVal{std::numeric_limits<double>::infinity()};
    double maxVal{-std::numeric_limits<double>::infinity()};

    void add(double v) {
        if (!std::isfinite(v)) return;
        ++n;
        sum += v;
        sumSq += v * v;
        if (v < minVal) minVal = v;
        if (v > maxVal) maxVal = v;
    }

    double mean() const {
        return n ? sum / static_cast<double>(n) : 0.0;
    }

    double rms() const {
        return n ? std::sqrt(sumSq / static_cast<double>(n)) : 0.0;
    }
};

static double rawGpsClockBiasSec(const NavEphGPS &eph, const CommonTime &t) {
    const double dt = t - eph.ctToc;
    return eph.af0 + dt * (eph.af1 + dt * eph.af2);
}

static CommonTime computeTransmitTime(const NavEphGPS &eph,
                                      const CommonTime &recvGpsTime,
                                      double pseudorangeMeters,
                                      bool applyTgd) {
    CommonTime transmit = recvGpsTime;
    transmit -= pseudorangeMeters / C_MPS;

    CommonTime tt = transmit;
    for (int i = 0; i < 2; ++i) {
        double clkBias = rawGpsClockBiasSec(eph, tt);
        if (applyTgd) {
            clkBias -= eph.TGD;
        }
        double relCorr = eph.svRelativity(tt);
        tt = transmit;
        tt -= (clkBias + relCorr);
    }

    return tt;
}

static double computeIfPseudorange(const TypeValueMap &tv, const string &sys) {
    auto itC1 = tv.find("C1");
    auto itC2 = tv.find("C2");
    if (itC1 == tv.end() || itC2 == tv.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double f1 = getFreq(sys, "C1");
    double f2 = getFreq(sys, "C2");
    return (f1 * f1 * itC1->second - f2 * f2 * itC2->second) / (f1 * f1 - f2 * f2);
}

static string formatCivilTime(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    ostringstream oss;
    oss << civ.year << "-"
        << setw(2) << setfill('0') << civ.month << "-"
        << setw(2) << setfill('0') << civ.day << " "
        << setw(2) << setfill('0') << civ.hour << ":"
        << setw(2) << setfill('0') << civ.minute << ":"
        << fixed << setprecision(6) << setw(0) << civ.second;
    return oss.str();
}

static CommonTime toGpsTime(const CommonTime &ct) {
    CommonTime gpsTime = convertTimeSystem(ct, TimeSystem::GPS);
    gpsTime.setTimeSystem(TimeSystem::GPS);
    return gpsTime;
}

static void printRowSummary(const SatID &sat,
                            const CommonTime &epochGps,
                            const NavEphGPS &eph,
                            double c1,
                            double ifPr,
                            const CommonTime &txRaw,
                            const CommonTime &txTgd,
                            const CommonTime &txIf,
                            ostream &csv,
                            TxStats &tgdStats,
                            TxStats &ifStats) {
    const double tgdSec = eph.TGD;
    const double tgdMeters = C_MPS * tgdSec;
    const double dtTgdNs = (txTgd - txRaw) * 1.0e9;
    const double dtIfNs = (txIf - txTgd) * 1.0e9;
    const double txRawOffset = txRaw - epochGps;
    const double txTgdOffset = txTgd - epochGps;
    const double txIfOffset = txIf - epochGps;

    tgdStats.add(dtTgdNs);
    ifStats.add(dtIfNs);

    csv << formatCivilTime(epochGps) << ","
        << sat << ","
        << fixed << setprecision(3)
        << c1 << ","
        << ifPr << ","
        << scientific << setprecision(12)
        << tgdSec << ","
        << tgdMeters << ","
        << txRawOffset << ","
        << txTgdOffset << ","
        << txIfOffset << ","
        << dtTgdNs << ","
        << dtIfNs << "\n";
}

} // namespace

int main() {
    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\Zero-baseline\\";
    string obsFile = dataDir + "oem719-202203031500-1.obs";
    string navFile = dataDir + "BRDC00IGS_R_20220620000_01D_MN.rnx";
    string csvFile = dataDir + "exam-5.2-gps_tx_time.csv";

    fstream roverObsStream(obsFile);
    if (!roverObsStream) {
        cerr << "rover file open error!" << strerror(errno) << endl;
        return -1;
    }

    RinexNavStore navStore;
    navStore.loadFile(navFile);

    map<string, set<string>> selectedTypes;
    selectedTypes["G"].insert("C1C");
    selectedTypes["G"].insert("C2W");

    RinexObsReader readObsRover;
    readObsRover.setFileStream(&roverObsStream);
    readObsRover.setSelectedTypes(selectedTypes);

    ofstream csv(csvFile);
    if (!csv) {
        cerr << "csv file open error!" << strerror(errno) << endl;
        return -1;
    }

    csv << "epoch,sat,c1_m,if_m,tgd_s,tgd_m,tx_raw_offset_s,tx_tgd_offset_s,tx_if_offset_s,dt_tgd_ns,dt_if_ns\n";

    TxStats tgdStats;
    TxStats ifStats;
    size_t printedRows = 0;

    while (true) {
        ObsData roverData;
        try {
            roverData = readObsRover.parseRinexObs();
        } catch (EndOfFile &) {
            break;
        }

        convertObsType(roverData);

        CommonTime epochGps = toGpsTime(roverData.epoch);

        for (const auto &satEntry : roverData.satTypeValueData) {
            const SatID &sat = satEntry.first;
            if (sat.system != "G") {
                continue;
            }

            const TypeValueMap &tv = satEntry.second;
            auto itC1 = tv.find("C1");
            auto itC2 = tv.find("C2");
            if (itC1 == tv.end() || itC2 == tv.end()) {
                continue;
            }

            NavEphGPS eph = navStore.findGPSEph(sat, epochGps);
            if (!eph.isValid(epochGps)) {
                continue;
            }

            const double c1 = itC1->second;
            const double ifPr = computeIfPseudorange(tv, sat.system);
            if (!std::isfinite(ifPr)) {
                continue;
            }

            const CommonTime txRaw = computeTransmitTime(eph, epochGps, c1, false);
            const CommonTime txTgd = computeTransmitTime(eph, epochGps, c1, true);
            const CommonTime txIf = computeTransmitTime(eph, epochGps, ifPr, true);

            printRowSummary(sat, epochGps, eph, c1, ifPr, txRaw, txTgd, txIf,
                            csv, tgdStats, ifStats);

            if (printedRows < 12) {
                cout << "epoch=" << formatCivilTime(epochGps)
                     << " sat=" << sat
                     << " C1=" << fixed << setprecision(3) << c1
                     << " IF=" << ifPr
                     << " TGD(s)=" << scientific << setprecision(12) << eph.TGD
                     << " dTx_TGD(ns)=" << (txTgd - txRaw) * 1.0e9
                     << " dTx_IF(ns)=" << (txIf - txTgd) * 1.0e9
                     << endl;
                ++printedRows;
            }
        }
    }

    cout << "\n--- Summary ---" << endl;
    cout << "TGD transmit-time shift(ns): mean=" << tgdStats.mean()
         << " rms=" << tgdStats.rms()
         << " min=" << tgdStats.minVal
         << " max=" << tgdStats.maxVal << endl;
    cout << "IF transmit-time shift(ns): mean=" << ifStats.mean()
         << " rms=" << ifStats.rms()
         << " min=" << ifStats.minVal
         << " max=" << ifStats.maxVal << endl;
    cout << "CSV written to: " << csvFile << endl;

    roverObsStream.close();
    csv.close();

    return 0;
}
