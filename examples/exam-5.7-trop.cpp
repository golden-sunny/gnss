#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Const.h"
#ifndef M_PI
#define M_PI PI
#endif
#include "CoordConvert.h"
#include "CoordStruct.h"
#include "GnssStruct.h"
#include "RinexNavStore.hpp"
#include "TimeConvert.h"

#include "navreadheader.h"
#include "navreader.h"
#include "obsreadheader.h"
#include "obsreader.h"

using namespace std;

namespace {

constexpr double kElevationCutoffDeg = 10.0;
constexpr double kRelativeHumidity = 0.5;

struct ExampleInput {
    XYZ receiverXyz;
    double azimuthDeg;
    double elevationDeg;
    double relativeHumidity;
};

struct Geodetic {
    double latRad{0.0};
    double lonRad{0.0};
    double hMeters{0.0};
};

struct AzEl {
    double azRad{0.0};
    double elRad{0.0};
};

struct TropResult {
    double elevationRad{0.0};
    double latitudeRad{0.0};
    double longitudeRad{0.0};
    double heightMeters{0.0};
    double pressureHpa{0.0};
    double temperatureKelvin{0.0};
    double saturationVaporHpa{0.0};
    double vaporPressureHpa{0.0};
    double saastamoinenFactor{0.0};
    double zhdMeters{0.0};
    double zwdMeters{0.0};
    double mappingFunction{0.0};
    double ztdMeters{0.0};
    double slantDelayMeters{0.0};
};

struct ReferenceValue {
    string name;
    double computed;
    double reference;
    string unit;
    int precision;
};

struct TropRow {
    CommonTime epochGps;
    SatID sat;
    string primaryCodeType;
    double elevationDeg{0.0};
    double azimuthDeg{0.0};
    double pressureHpa{0.0};
    double temperatureKelvin{0.0};
    double vaporPressureHpa{0.0};
    double zhdMeters{0.0};
    double zwdMeters{0.0};
    double mappingFunction{0.0};
    double ztdMeters{0.0};
    double slantDelayMeters{0.0};
};

struct EpochSummary {
    CommonTime epochGps;
    size_t satCount{0};
    double meanElevationDeg{0.0};
    double meanTropMeters{0.0};
    double minTropMeters{0.0};
    double maxTropMeters{0.0};
};

struct MetricStats {
    size_t n{0};
    double mean{0.0};
    double minVal{0.0};
    double maxVal{0.0};
    CommonTime minEpoch;
    CommonTime maxEpoch;
};

static bool finitePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

static double clampUnit(double x) {
    return max(-1.0, min(1.0, x));
}

static double wrapTwoPi(double x) {
    const double twoPi = 2.0 * PI;
    x = fmod(x, twoPi);
    if (x < 0.0) x += twoPi;
    return x;
}

static CommonTime toGpsTime(const RnxTime &rt) {
    CivilTime civ(rt.year, rt.month, rt.day, rt.hour, rt.minute, rt.second, TimeSystem::GPS);
    CommonTime ct = CivilTime2CommonTime(civ);
    ct.setTimeSystem(TimeSystem::GPS);
    return ct;
}

static string formatCivilTime(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    ostringstream oss;
    oss << civ.year << "-"
        << setw(2) << setfill('0') << civ.month << "-"
        << setw(2) << civ.day << " "
        << setw(2) << civ.hour << ":"
        << setw(2) << civ.minute << ":"
        << fixed << setprecision(3) << setw(6) << civ.second;
    return oss.str();
}

static Geodetic ecefToGeodetic(const Eigen::Vector3d &xyz) {
    WGS84 wgs84;
    const BLH blh = xyz2blh(XYZ(xyz), wgs84);

    Geodetic out;
    out.latRad = blh.B();
    out.lonRad = blh.L();
    out.hMeters = blh.H();
    return out;
}

static AzEl computeAzEl(const Eigen::Vector3d &recXYZ,
                        const Geodetic &recBLH,
                        const Eigen::Vector3d &satXYZ) {
    const Eigen::Vector3d los = (satXYZ - recXYZ).normalized();
    const double lat = recBLH.latRad;
    const double lon = recBLH.lonRad;

    const Eigen::Vector3d east(-sin(lon), cos(lon), 0.0);
    const Eigen::Vector3d north(-sin(lat) * cos(lon),
                                -sin(lat) * sin(lon),
                                cos(lat));
    const Eigen::Vector3d up(cos(lat) * cos(lon),
                             cos(lat) * sin(lon),
                             sin(lat));

    AzEl out;
    out.azRad = wrapTwoPi(atan2(los.dot(east), los.dot(north)));
    out.elRad = asin(clampUnit(los.dot(up)));
    return out;
}

static Eigen::Vector3d rotateSatelliteToReceive(const Eigen::Vector3d &recXYZ,
                                                const Xvt &xvtTx) {
    const double travelSec = (xvtTx.x - recXYZ).norm() / C_MPS;
    const double wt = OMEGA_EARTH * travelSec;
    Eigen::Vector3d rot;
    rot.x() = cos(wt) * xvtTx.x.x() + sin(wt) * xvtTx.x.y();
    rot.y() = -sin(wt) * xvtTx.x.x() + cos(wt) * xvtTx.x.y();
    rot.z() = xvtTx.x.z();
    return rot;
}

static void loadGpsNavIntoStore(NavReader &navReader, RinexNavStore &navStore, size_t &count) {
    count = 0;
    while (true) {
        try {
            NavGpsRecord rec = navReader.readNextGpsRecord();
            double fitHours = rec.eph.fitInterval;
            if (!std::isfinite(fitHours) || fitHours <= 0.0) {
                fitHours = 4.0;
            }
            const double halfFitSec = 0.5 * fitHours * 3600.0;
            rec.eph.beginValid = rec.eph.ctToe - halfFitSec;
            rec.eph.endValid = rec.eph.ctToe + halfFitSec;
            navStore.gpsEphData[rec.sat][rec.eph.ctToe] = rec.eph;
            ++count;
        } catch (const EndOfNavFile &) {
            break;
        }
    }
}

static Xvt computeTransmitXvt(RinexNavStore &navStore,
                              const SatID &sat,
                              const CommonTime &recvGps,
                              double pseudorangeMeters) {
    CommonTime tx = recvGps;
    tx -= pseudorangeMeters / C_MPS;
    Xvt xvt = navStore.getXvt(sat, tx);
    tx -= (xvt.clkbias + xvt.relcorr);
    return navStore.getXvt(sat, tx);
}

static string choosePrimaryCodeType(const map<string, double> &tv) {
    const vector<string> preferred = {"C1C", "C1W", "C1X"};
    for (const auto &type : preferred) {
        auto it = tv.find(type);
        if (it != tv.end() && finitePositive(it->second)) {
            return type;
        }
    }
    return "";
}

static double standardPressureHpa(double heightMeters) {
    return 1013.25 * pow(1.0 - 0.0000226 * heightMeters, 5.225);
}

static double standardTemperatureKelvin(double heightMeters) {
    const double temperatureCelsius = 15.0 - 0.0065 * heightMeters;
    return temperatureCelsius + 273.15;
}

static double saturationVaporPressureHpa(double temperatureKelvin) {
    return 6.108 * exp((17.15 * temperatureKelvin - 4684.0)
                       / (temperatureKelvin - 38.45));
}

static double saastamoinenFactor(double latitudeRad, double heightMeters) {
    const double heightKm = heightMeters / 1000.0;
    return 1.0 - 0.266e-2 * cos(2.0 * latitudeRad) - 0.28e-3 * heightKm;
}

static double zenithHydrostaticDelay(double pressureHpa,
                                     double latitudeRad,
                                     double heightMeters) {
    return 2.277e-3 * pressureHpa / saastamoinenFactor(latitudeRad, heightMeters);
}

static double zenithWetDelay(double temperatureKelvin, double vaporPressureHpa) {
    return 0.002277 * (1255.0 / temperatureKelvin + 0.05) * vaporPressureHpa;
}

static double simpleMappingFunction(double elevationRad) {
    double sinElev = sin(elevationRad);
    if (sinElev < 0.05) {
        sinElev = 0.05;
    }
    return 1.0 / sinElev;
}

static TropResult computeSaastamoinen(const Geodetic &recBLH,
                                      double elevationDeg,
                                      double relativeHumidity) {
    TropResult result{};
    result.elevationRad = elevationDeg * DEG_TO_RAD;
    result.latitudeRad = recBLH.latRad;
    result.longitudeRad = recBLH.lonRad;
    result.heightMeters = recBLH.hMeters;
    result.pressureHpa = standardPressureHpa(result.heightMeters);
    result.temperatureKelvin = standardTemperatureKelvin(result.heightMeters);
    result.saturationVaporHpa = saturationVaporPressureHpa(result.temperatureKelvin);
    result.vaporPressureHpa = relativeHumidity * result.saturationVaporHpa;
    result.saastamoinenFactor = saastamoinenFactor(result.latitudeRad, result.heightMeters);
    result.zhdMeters = zenithHydrostaticDelay(result.pressureHpa,
                                              result.latitudeRad,
                                              result.heightMeters);
    result.zwdMeters = zenithWetDelay(result.temperatureKelvin, result.vaporPressureHpa);
    result.mappingFunction = simpleMappingFunction(result.elevationRad);
    result.ztdMeters = result.zhdMeters + result.zwdMeters;
    result.slantDelayMeters = result.ztdMeters * result.mappingFunction;
    return result;
}

static TropResult reproduceExample53(const ExampleInput &input) {
    WGS84 wgs84;
    const BLH blh = xyz2blh(input.receiverXyz, wgs84);

    Geodetic recBLH;
    recBLH.latRad = blh.B();
    recBLH.lonRad = blh.L();
    recBLH.hMeters = blh.H();
    return computeSaastamoinen(recBLH, input.elevationDeg, input.relativeHumidity);
}

static void printInputs(const ExampleInput &input) {
    cout << "Exercise 5.7: tropospheric delay with the Saastamoinen model" << endl;
    cout << "Reproducing textbook example 5-3 and extending it to one-day processing" << endl << endl;

    cout << "Input data" << endl;
    cout << "  Receiver XYZ (m): "
         << fixed << setprecision(4)
         << "(" << input.receiverXyz.X() << ", "
         << input.receiverXyz.Y() << ", "
         << input.receiverXyz.Z() << ")" << endl;
    cout << "  Satellite azimuth (deg):   " << setprecision(6) << input.azimuthDeg << endl;
    cout << "  Satellite elevation (deg): " << setprecision(6) << input.elevationDeg << endl;
    cout << "  Relative humidity:         " << setprecision(3) << input.relativeHumidity << endl
         << endl;
}

static void printIntermediateResults(const TropResult &result) {
    cout << "Computed intermediate values" << endl;
    cout << "  Elevation angle (rad):       " << fixed << setprecision(6)
         << result.elevationRad << endl;
    cout << "  Latitude (rad):              " << result.latitudeRad << endl;
    cout << "  Longitude (deg):             " << setprecision(6)
         << result.longitudeRad * RAD_TO_DEG << endl;
    cout << "  Ellipsoidal height (m):      " << setprecision(6)
         << result.heightMeters << endl;
    cout << "  Pressure (hPa):              " << setprecision(6)
         << result.pressureHpa << endl;
    cout << "  Temperature (K):             " << setprecision(6)
         << result.temperatureKelvin << endl;
    cout << "  Saturation vapor pressure:   " << setprecision(6)
         << result.saturationVaporHpa << " hPa" << endl;
    cout << "  Actual vapor pressure (hPa): " << setprecision(6)
         << result.vaporPressureHpa << endl;
    cout << "  f(B,H):                      " << setprecision(9)
         << result.saastamoinenFactor << endl;
    cout << "  ZHD (m):                     " << setprecision(6)
         << result.zhdMeters << endl;
    cout << "  ZWD (m):                     " << setprecision(6)
         << result.zwdMeters << endl;
    cout << "  Mapping function:            " << setprecision(9)
         << result.mappingFunction << endl;
    cout << "  Zenith total delay (m):      " << setprecision(6)
         << result.ztdMeters << endl;
    cout << "  Slant tropospheric delay (m): " << setprecision(6)
         << result.slantDelayMeters << endl
         << endl;
}

static void printComparison(const TropResult &result) {
    const vector<ReferenceValue> items = {
        {"Elevation angle", result.elevationRad, 1.448750, "rad", 6},
        {"Latitude", result.latitudeRad, 0.532878, "rad", 6},
        {"Ellipsoidal height", result.heightMeters, 25.126182, "m", 6},
        {"Pressure", result.pressureHpa, 1010.247266, "hPa", 6},
        {"Temperature", result.temperatureKelvin, 287.986680, "K", 6},
        {"Vapor pressure", result.vaporPressureHpa, 8.484425, "hPa", 6},
        {"ZHD", result.zhdMeters, 2.303314, "m", 6},
        {"ZWD", result.zwdMeters, 0.085155, "m", 6},
        {"Relative humidity", 0.5, 0.5, "", 3},
        {"Slant trop delay", result.slantDelayMeters, 2.406368, "m", 6}
    };

    cout << "Comparison with Table 5-3" << endl;
    cout << left << setw(22) << "Parameter"
         << right << setw(18) << "Computed"
         << setw(18) << "Reference"
         << setw(18) << "Difference"
         << setw(10) << "Unit" << endl;
    cout << string(86, '-') << endl;

    for (const auto &item : items) {
        cout << left << setw(22) << item.name
             << right << setw(18) << fixed << setprecision(item.precision) << item.computed
             << setw(18) << fixed << setprecision(item.precision) << item.reference
             << setw(18) << scientific << setprecision(6) << (item.computed - item.reference)
             << setw(10) << item.unit << endl;
    }
    cout << endl;
}

static vector<EpochSummary> summarizeByEpoch(const vector<TropRow> &rows) {
    struct Accumulator {
        size_t satCount{0};
        double elevationSum{0.0};
        double tropSum{0.0};
        double tropMin{numeric_limits<double>::infinity()};
        double tropMax{-numeric_limits<double>::infinity()};
    };

    map<CommonTime, Accumulator> acc;
    for (const auto &row : rows) {
        auto &item = acc[row.epochGps];
        ++item.satCount;
        item.elevationSum += row.elevationDeg;
        item.tropSum += row.slantDelayMeters;
        item.tropMin = min(item.tropMin, row.slantDelayMeters);
        item.tropMax = max(item.tropMax, row.slantDelayMeters);
    }

    vector<EpochSummary> out;
    out.reserve(acc.size());
    for (const auto &entry : acc) {
        EpochSummary row;
        row.epochGps = entry.first;
        row.satCount = entry.second.satCount;
        row.meanElevationDeg = entry.second.elevationSum / static_cast<double>(entry.second.satCount);
        row.meanTropMeters = entry.second.tropSum / static_cast<double>(entry.second.satCount);
        row.minTropMeters = entry.second.tropMin;
        row.maxTropMeters = entry.second.tropMax;
        out.push_back(row);
    }
    return out;
}

static MetricStats computeMetricStats(const vector<TropRow> &rows) {
    MetricStats stats;
    if (rows.empty()) return stats;

    stats.n = rows.size();
    stats.minVal = rows.front().slantDelayMeters;
    stats.maxVal = rows.front().slantDelayMeters;
    stats.minEpoch = rows.front().epochGps;
    stats.maxEpoch = rows.front().epochGps;

    double sum = 0.0;
    for (const auto &row : rows) {
        sum += row.slantDelayMeters;
        if (row.slantDelayMeters < stats.minVal) {
            stats.minVal = row.slantDelayMeters;
            stats.minEpoch = row.epochGps;
        }
        if (row.slantDelayMeters > stats.maxVal) {
            stats.maxVal = row.slantDelayMeters;
            stats.maxEpoch = row.epochGps;
        }
    }
    stats.mean = sum / static_cast<double>(stats.n);
    return stats;
}

static void writeDetailCsv(const string &csvFile, const vector<TropRow> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open csv file: " + csvFile);

    csv << "epoch,sat,primary_code_type,elevation_deg,azimuth_deg,pressure_hpa,"
        << "temperature_k,vapor_pressure_hpa,zhd_m,zwd_m,mapping_function,ztd_m,slant_delay_m\n";
    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.sat << ","
            << row.primaryCodeType << ","
            << fixed << setprecision(6)
            << row.elevationDeg << ","
            << row.azimuthDeg << ","
            << row.pressureHpa << ","
            << row.temperatureKelvin << ","
            << row.vaporPressureHpa << ","
            << row.zhdMeters << ","
            << row.zwdMeters << ","
            << row.mappingFunction << ","
            << row.ztdMeters << ","
            << row.slantDelayMeters << "\n";
    }
}

static void writeSummaryCsv(const string &csvFile, const vector<EpochSummary> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open summary csv file: " + csvFile);

    csv << "epoch,sat_count,mean_elevation_deg,mean_trop_m,min_trop_m,max_trop_m\n";
    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.satCount << ","
            << fixed << setprecision(6)
            << row.meanElevationDeg << ","
            << row.meanTropMeters << ","
            << row.minTropMeters << ","
            << row.maxTropMeters << "\n";
    }
}

} // namespace

int main() {
    const ExampleInput input{
        XYZ(-2267749.0, 5009154.0, 3221290.0),
        265.659363,
        83.007251,
        0.5
    };

    printInputs(input);
    const TropResult book = reproduceExample53(input);
    printIntermediateResults(book);
    printComparison(book);

    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\";
    const string obsFile = dataDir + "WUH200CHN_R_20250010000_01D_30S_MO.rnx";
    const string navFile = dataDir + "BRDC00IGS_R_20250010000_01D_MN.rnx";
    const string detailCsvFile = dataDir + "exam-5.7-trop.csv";
    const string summaryCsvFile = dataDir + "exam-5.7-trop_summary.csv";

    fstream obsStream(obsFile);
    if (!obsStream) {
        cerr << "obs file open error: " << strerror(errno) << endl;
        return -1;
    }

    fstream navStream(navFile);
    if (!navStream) {
        cerr << "nav file open error: " << strerror(errno) << endl;
        return -1;
    }

    ObsReadHeader obsHeaderReader;
    obsHeaderReader.setFileStream(&obsStream);
    obsHeaderReader.parseHeader();
    const RinexObsHeaderAll &obsHeader = obsHeaderReader.getHeader();

    NavReadHeader navHeaderReader;
    navHeaderReader.setFileStream(&navStream);
    navHeaderReader.parseHeader();
    const NavHeaderAll &navHeader = navHeaderReader.getHeader();
    (void)navHeader;

    ObsReader obsReader;
    obsReader.setFileStream(&obsStream);
    obsReader.setHeader(&obsHeader);

    NavReader navReader;
    navReader.setFileStream(&navStream);
    navReader.setHeader(&navHeader);

    RinexNavStore navStore;
    size_t gpsRecordCount = 0;
    loadGpsNavIntoStore(navReader, navStore, gpsRecordCount);

    const Eigen::Vector3d recXYZ(obsHeader.approxX, obsHeader.approxY, obsHeader.approxZ);
    const Geodetic recBLH = ecefToGeodetic(recXYZ);

    cout << "--- Daily troposphere processing ---" << endl;
    cout << "Observation file: " << obsFile << endl;
    cout << "Navigation file : " << navFile << endl;
    cout << "Station         : " << obsHeader.markerName << endl;
    cout << "Receiver BLH    : lat=" << fixed << setprecision(8) << recBLH.latRad * RAD_TO_DEG
         << " deg, lon=" << recBLH.lonRad * RAD_TO_DEG
         << " deg, h=" << setprecision(3) << recBLH.hMeters << " m" << endl;
    cout << "GPS nav records : " << gpsRecordCount << endl;
    cout << "Elevation cutoff: " << kElevationCutoffDeg << " deg" << endl;
    cout << "Relative humidity used in Saastamoinen: " << kRelativeHumidity << endl << endl;

    vector<TropRow> rows;
    size_t totalEpochs = 0;
    size_t totalGpsSeen = 0;
    size_t rejectedNoCode = 0;
    size_t rejectedNav = 0;
    size_t rejectedLowElev = 0;

    while (true) {
        ObsEpochData epochData;
        try {
            epochData = obsReader.readNextEpoch();
        } catch (const EndOfObsFile &) {
            break;
        }

        ++totalEpochs;
        const CommonTime epochGps = toGpsTime(epochData.time);

        for (const auto &satEntry : epochData.satTypeValue) {
            const SatID sat(satEntry.first);
            if (sat.system != "G") continue;
            ++totalGpsSeen;

            const auto &tv = satEntry.second;
            const string primaryCodeType = choosePrimaryCodeType(tv);
            if (primaryCodeType.empty()) {
                ++rejectedNoCode;
                continue;
            }

            Xvt xvtTx;
            try {
                xvtTx = computeTransmitXvt(navStore, sat, epochGps, tv.at(primaryCodeType));
            } catch (...) {
                ++rejectedNav;
                continue;
            }

            const Eigen::Vector3d satXYZRec = rotateSatelliteToReceive(recXYZ, xvtTx);
            const AzEl azel = computeAzEl(recXYZ, recBLH, satXYZRec);
            const double elevDeg = azel.elRad * RAD_TO_DEG;
            if (elevDeg < kElevationCutoffDeg) {
                ++rejectedLowElev;
                continue;
            }

            const TropResult trop = computeSaastamoinen(recBLH, elevDeg, kRelativeHumidity);

            TropRow row;
            row.epochGps = epochGps;
            row.sat = sat;
            row.primaryCodeType = primaryCodeType;
            row.elevationDeg = elevDeg;
            row.azimuthDeg = azel.azRad * RAD_TO_DEG;
            row.pressureHpa = trop.pressureHpa;
            row.temperatureKelvin = trop.temperatureKelvin;
            row.vaporPressureHpa = trop.vaporPressureHpa;
            row.zhdMeters = trop.zhdMeters;
            row.zwdMeters = trop.zwdMeters;
            row.mappingFunction = trop.mappingFunction;
            row.ztdMeters = trop.ztdMeters;
            row.slantDelayMeters = trop.slantDelayMeters;
            rows.push_back(row);
        }
    }

    sort(rows.begin(), rows.end(), [](const TropRow &a, const TropRow &b) {
        const double dt = a.epochGps - b.epochGps;
        if (abs(dt) > 1.0e-9) return dt < 0.0;
        return a.sat < b.sat;
    });

    const vector<EpochSummary> summaryRows = summarizeByEpoch(rows);
    const MetricStats stats = computeMetricStats(rows);

    writeDetailCsv(detailCsvFile, rows);
    writeSummaryCsv(summaryCsvFile, summaryRows);

    cout << "Epochs read              : " << totalEpochs << endl;
    cout << "GPS observations seen    : " << totalGpsSeen << endl;
    cout << "Accepted GPS samples     : " << rows.size() << endl;
    cout << "Rejected no primary code : " << rejectedNoCode << endl;
    cout << "Rejected no nav solution : " << rejectedNav << endl;
    cout << "Rejected low elevation   : " << rejectedLowElev << endl;
    cout << fixed << setprecision(3)
         << "Slant trop mean/min/max  : " << stats.mean << " / " << stats.minVal
         << " / " << stats.maxVal << " m" << endl;
    if (stats.n > 0) {
        cout << "  min epoch: " << formatCivilTime(stats.minEpoch)
             << ", max epoch: " << formatCivilTime(stats.maxEpoch) << endl;
    }
    cout << "Detail CSV written to    : " << detailCsvFile << endl;
    cout << "Summary CSV written to   : " << summaryCsvFile << endl;

    obsStream.close();
    navStream.close();
    return 0;
}
