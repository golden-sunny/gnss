/**
 * Exercise 5.3 test.
 *
 * Reproduce the textbook GPS G10 case at 2025-01-01 00:00:00 GPST:
 * - use the WUH2 rover observation file to fetch the C1C pseudorange
 * - use the broadcast navigation file to iterate satellite transmit time twice
 * - compare the computed values with the textbook reference values
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

#include "Const.h"
#include "GnssFunc.h"
#include "NavEphGPS.hpp"
#include "RinexNavStore.hpp"
#include "RinexObsReader.h"
#include "SPPIFCode.h"
#include "TimeConvert.h"

using namespace std;

namespace {

struct TxIterResult {
    CommonTime tx;
    Xvt xvt;
};

struct BookCalcResult {
    CommonTime recvGps;
    double c1{0.0};
    Eigen::Vector3d recXYZ{Eigen::Vector3d::Zero()};
    TxIterResult iter1;
    TxIterResult iter2;
    Xvt earthRotXvt;
};

struct BookTheory {
    CommonTime recvGps;
    double c1{0.0};
    CommonTime tx1;
    CommonTime tx2;
    Eigen::Vector3d x1;
    Eigen::Vector3d v1;
    double clk1{0.0};
    double rel1{0.0};
    Eigen::Vector3d x2;
    Eigen::Vector3d v2;
    double clk2{0.0};
    double rel2{0.0};
    Eigen::Vector3d recXYZ;
    Eigen::Vector3d earthRotX;
    Eigen::Vector3d earthRotV;
};

static CommonTime toGpsTime(const CommonTime &ct) {
    CommonTime gpsTime = convertTimeSystem(ct, TimeSystem::GPS);
    gpsTime.setTimeSystem(TimeSystem::GPS);
    return gpsTime;
}

static string formatBookCommonTime(const CommonTime &ct) {
    long mjd = 0;
    double sod = 0.0;
    TimeSystem ts;
    ct.get(mjd, sod, ts);

    const long jdDay = mjd + static_cast<long>(MJD_TO_JD);
    ostringstream oss;
    oss << "(" << jdDay << ", "
        << fixed << setprecision(18) << sod
        << ", " << ts.toString() << ")";
    return oss.str();
}

static string formatVec3(const Eigen::Vector3d &v, int precision = 3) {
    ostringstream oss;
    oss << fixed << setprecision(precision)
        << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return oss.str();
}

static string formatScalar(double v, int precision = 12) {
    ostringstream oss;
    oss << scientific << setprecision(precision) << v;
    return oss.str();
}

static string formatTimeDiffNs(const CommonTime &calc, const CommonTime &theory) {
    ostringstream oss;
    oss << scientific << setprecision(6)
        << (calc - theory) * 1.0e9;
    return oss.str();
}

static string formatVecDiff(const Eigen::Vector3d &calc, const Eigen::Vector3d &theory, int precision = 3) {
    ostringstream oss;
    oss << fixed << setprecision(precision)
        << "(" << (calc.x() - theory.x())
        << ", " << (calc.y() - theory.y())
        << ", " << (calc.z() - theory.z()) << ")";
    return oss.str();
}

static void printCompareHeader() {
    cout << left
         << setw(18) << "Item"
         << setw(44) << "Computed"
         << setw(44) << "Theory"
         << "Diff" << endl;
    cout << string(120, '-') << endl;
}

static void printTimeCompare(const string &label,
                             const CommonTime &calc,
                             const CommonTime &theory) {
    cout << left
         << setw(18) << label
         << setw(44) << formatBookCommonTime(calc)
         << setw(44) << formatBookCommonTime(theory)
         << formatTimeDiffNs(calc, theory) << " ns" << endl;
}

static void printScalarCompare(const string &label,
                               double calc,
                               double theory,
                               int precision = 12) {
    cout << left
         << setw(18) << label
         << setw(44) << formatScalar(calc, precision)
         << setw(44) << formatScalar(theory, precision)
         << formatScalar(calc - theory, precision) << endl;
}

static void printVecCompare(const string &label,
                            const Eigen::Vector3d &calc,
                            const Eigen::Vector3d &theory,
                            int precision = 3) {
    cout << left
         << setw(18) << label
         << setw(44) << formatVec3(calc, precision)
         << setw(44) << formatVec3(theory, precision)
         << formatVecDiff(calc, theory, precision) << endl;
}

static BookCalcResult computeBookCase(RinexNavStore &navStore,
                                      const SatID &sat,
                                      const CommonTime &recvGps,
                                      double c1Meters,
                                      const Eigen::Vector3d &recXYZ) {
    BookCalcResult result;
    result.recvGps = recvGps;
    result.c1 = c1Meters;
    result.recXYZ = recXYZ;

    result.iter1.tx = recvGps;
    result.iter1.tx -= c1Meters / C_MPS;
    result.iter1.xvt = navStore.getXvt(sat, result.iter1.tx);

    result.iter2.tx = result.iter1.tx;
    result.iter2.tx -= (result.iter1.xvt.clkbias + result.iter1.xvt.relcorr);
    result.iter2.xvt = navStore.getXvt(sat, result.iter2.tx);

    SPPIFCode spp;
    std::map<SatID, Xvt> satXvtTransTime;
    satXvtTransTime[sat] = result.iter2.xvt;
    Eigen::Vector3d recXYZCopy = recXYZ;
    std::map<SatID, Xvt> satXvtRecTime = spp.earthRotation(recXYZCopy, satXvtTransTime);
    result.earthRotXvt = satXvtRecTime.at(sat);

    return result;
}

static BookTheory getBookTheory() {
    BookTheory theory;
    theory.recvGps = CivilTime2CommonTime(CivilTime(2025, 1, 1, 0, 0, 0, TimeSystem::GPS));
    theory.c1 = 20183692.238000;
    // The screenshot uses a JD label that is one day ahead for the transmit epoch.
    // Here we keep the absolute epoch used by CommonTime so the comparison is meaningful.
    theory.tx1 = CommonTime(60675, 86399.9326744501185883, TimeSystem::GPS);
    theory.tx2 = CommonTime(60675, 86399.932932118477765471, TimeSystem::GPS);
    theory.x1 = Eigen::Vector3d(-7192606.606, 21829534.875, 13240734.408);
    theory.v1 = Eigen::Vector3d(-1207.655, 1140.117, -2613.429);
    theory.clk1 = -0.0002576895;
    theory.rel1 = 0.0000000230;
    theory.x2 = Eigen::Vector3d(-7192606.917, 21829535.169, 13240733.735);
    theory.v2 = Eigen::Vector3d(-1207.655, 1140.117, -2613.429);
    theory.clk2 = -0.0002576895;
    theory.rel2 = 0.0000000230;
    theory.recXYZ = Eigen::Vector3d(-2267749.0000, 5009154.0000, 3221290.0000);
    theory.earthRotX = Eigen::Vector3d(-7192499.721, 21829570.488, 13240733.735);
    theory.earthRotV = Eigen::Vector3d(-1207.649, 1140.123, -2613.429);
    return theory;
}

} // namespace

int main() {
    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\";
    const string obsFile = dataDir + "WUH200CHN_R_20250010000_01D_30S_MO.rnx";
    string navFile = dataDir + "BRDC00IGS_R_20250010000_01D_MN.rnx";

    fstream roverObsStream(obsFile);
    if (!roverObsStream) {
        cerr << "rover file open error: " << strerror(errno) << endl;
        return -1;
    }

    RinexNavStore navStore;
    navStore.loadFile(navFile);

    map<string, set<string>> selectedTypes;
    selectedTypes["G"].insert("C1C");

    RinexObsReader readObsRover;
    readObsRover.setFileStream(&roverObsStream);
    readObsRover.setSelectedTypes(selectedTypes);

    const BookTheory theory = getBookTheory();
    const SatID targetSat("G10");
    const double epochEpsSec = 1.0e-6;
    bool found = false;
    BookCalcResult calc;

    while (true) {
        ObsData roverData;
        try {
            roverData = readObsRover.parseRinexObs();
        } catch (EndOfFile &) {
            break;
        }

        convertObsType(roverData);
        CommonTime epochGps = toGpsTime(roverData.epoch);

        if (std::fabs(epochGps - theory.recvGps) > epochEpsSec) {
            continue;
        }

        auto satIt = roverData.satTypeValueData.find(targetSat);
        if (satIt == roverData.satTypeValueData.end()) {
            continue;
        }

        auto c1It = satIt->second.find("C1");
        if (c1It == satIt->second.end()) {
            continue;
        }

        calc = computeBookCase(navStore, targetSat, epochGps, c1It->second, roverData.antennaPosition);
        found = true;
        break;
    }

    if (!found) {
        cerr << "failed to find G10 at 2025-01-01 00:00:00 GPST" << endl;
        return -1;
    }

    cout << "Exercise 5.3 test: G10 transmit-time reproduction" << endl;
    cout << "Observation file: " << obsFile << endl;
    cout << "Navigation file : " << navFile << endl;
    cout << "Receive epoch   : " << formatBookCommonTime(calc.recvGps) << endl;
    cout << "C1C pseudorange : " << fixed << setprecision(3) << calc.c1 << " m" << endl;
    cout << "Receiver XYZ    : " << formatVec3(calc.recXYZ, 3) << endl;
    cout << endl;

    cout << "[Iteration 1]" << endl;
    printCompareHeader();
    printTimeCompare("Transmit time", calc.iter1.tx, theory.tx1);
    printVecCompare("Satellite XYZ", calc.iter1.xvt.x, theory.x1, 3);
    printVecCompare("Satellite VEL", calc.iter1.xvt.v, theory.v1, 3);
    printScalarCompare("Clock bias(s)", calc.iter1.xvt.clkbias, theory.clk1, 10);
    printScalarCompare("Relativity(s)", calc.iter1.xvt.relcorr, theory.rel1, 9);
    cout << endl;

    cout << "[Iteration 2]" << endl;
    printCompareHeader();
    printTimeCompare("Transmit time", calc.iter2.tx, theory.tx2);
    printVecCompare("Satellite XYZ", calc.iter2.xvt.x, theory.x2, 3);
    printVecCompare("Satellite VEL", calc.iter2.xvt.v, theory.v2, 3);
    printScalarCompare("Clock bias(s)", calc.iter2.xvt.clkbias, theory.clk2, 10);
    printScalarCompare("Relativity(s)", calc.iter2.xvt.relcorr, theory.rel2, 9);
    cout << endl;

    cout << "[Earth rotation]" << endl;
    printCompareHeader();
    printVecCompare("Receiver XYZ", calc.recXYZ, theory.recXYZ, 3);
    printVecCompare("Rotated XYZ", calc.earthRotXvt.x, theory.earthRotX, 3);
    printVecCompare("Rotated VEL", calc.earthRotXvt.v, theory.earthRotV, 3);

    roverObsStream.close();
    return 0;
}
