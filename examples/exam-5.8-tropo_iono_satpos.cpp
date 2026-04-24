/**
 * Exercise 5.8 tropo_iono_satpos.
 *
 * A compact GNSS observation-processing workflow that reuses the existing
 * navigation store, header readers, and textbook correction logic to:
 * - compute broadcast satellite state at transmit time
 * - apply Earth-rotation correction to the satellite position
 * - extract satellite clock / relativity / TGD terms
 * - estimate ionospheric delay with the GPS Klobuchar model
 * - estimate tropospheric delay with the Saastamoinen model
 * - optionally level dual-frequency GPS ionosphere for a model-vs-data check
 * - summarize and plot the daily ionosphere/troposphere variation
 */

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Const.h"
#ifndef M_PI
#define M_PI PI
#endif
#include "CoordConvert.h"
#include "CoordStruct.h"
#include "GnssFunc.h"
#include "GnssStruct.h"
#include "RinexNavStore.hpp"
#include "TimeConvert.h"

#include "navreadheader.h"
#include "navreader.h"
#include "obsreadheader.h"
#include "obsreader.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#endif

using namespace std;

namespace {

constexpr double kElevationCutoffDeg = 10.0;
constexpr double kRelativeHumidity = 0.5;
constexpr size_t kMinDualArcPoints = 10;
constexpr double kMaxDualModelDiffMeters = 15.0;

struct Geodetic {
    double latRad{0.0};
    double lonRad{0.0};
    double hMeters{0.0};
};

struct AzEl {
    double azRad{0.0};
    double elRad{0.0};
};

struct KlobucharCoeffs {
    array<double, 4> alpha{{0.0, 0.0, 0.0, 0.0}};
    array<double, 4> beta{{0.0, 0.0, 0.0, 0.0}};
};

struct KlobucharResult {
    double localTimeSec{0.0};
    double amplitudeSec{0.0};
    double periodSec{0.0};
    double mapping{0.0};
    double delaySec{0.0};
    double delayMeters{0.0};
};

struct TropResult {
    double pressureHpa{0.0};
    double temperatureKelvin{0.0};
    double saturationVaporHpa{0.0};
    double vaporPressureHpa{0.0};
    double zhdMeters{0.0};
    double zwdMeters{0.0};
    double ztdMeters{0.0};
    double mappingFunction{0.0};
    double slantDelayMeters{0.0};
};

struct TransmitSolution {
    CommonTime txGps;
    Xvt xvtTx;
    double recvMinusTxSec{0.0};
};

struct ObservationRow {
    CommonTime epochGps;
    SatID sat;
    string primaryCodeType;
    double primaryCodeRawMeters{0.0};
    double primaryCodeTgdCorrectedMeters{0.0};
    double tgdCorrectionMeters{0.0};
    double c1cMeters{numeric_limits<double>::quiet_NaN()};
    double c1wMeters{numeric_limits<double>::quiet_NaN()};
    double c2wMeters{numeric_limits<double>::quiet_NaN()};
    double l1wCycles{numeric_limits<double>::quiet_NaN()};
    double l2wCycles{numeric_limits<double>::quiet_NaN()};
    Eigen::Vector3d satPosTx{Eigen::Vector3d::Zero()};
    Eigen::Vector3d satVelTx{Eigen::Vector3d::Zero()};
    Eigen::Vector3d satPosRx{Eigen::Vector3d::Zero()};
    Eigen::Vector3d satVelRx{Eigen::Vector3d::Zero()};
    double transmitOffsetSec{0.0};
    double satClockBiasSec{0.0};
    double relativitySec{0.0};
    double tgdSec{0.0};
    double elevationDeg{0.0};
    double azimuthDeg{0.0};
    double modelIonoMeters{0.0};
    double modelTropMeters{0.0};
    bool hasDualRaw{false};
    double codeIonoMeters{numeric_limits<double>::quiet_NaN()};
    double phaseIonoMeters{numeric_limits<double>::quiet_NaN()};
    bool hasDualLeveled{false};
    double dualLeveledIonoMeters{numeric_limits<double>::quiet_NaN()};
    double dualMinusModelMeters{numeric_limits<double>::quiet_NaN()};
};

struct EpochSummary {
    CommonTime epochGps;
    size_t satCount{0};
    size_t dualCount{0};
    double meanElevationDeg{0.0};
    double meanIonoModelMeters{0.0};
    double meanDualIonoMeters{numeric_limits<double>::quiet_NaN()};
    double meanTropMeters{0.0};
    double meanClockBiasNs{0.0};
    double meanRelativityNs{0.0};
    double meanTgdNs{0.0};
};

struct MetricStats {
    size_t n{0};
    double mean{0.0};
    double minVal{0.0};
    double maxVal{0.0};
    CommonTime minEpoch;
    CommonTime maxEpoch;
};

struct SvgSeries {
    string label;
    string color;
    string fillColor;
    double lineWidth{2.0};
    double pointRadius{0.0};
    double opacity{1.0};
    vector<pair<double, double>> points;
};

struct ChartPanel {
    string title;
    string yLabel;
    vector<SvgSeries> seriesList;
    double yMin{0.0};
    double yMax{1.0};
};

static void padRange(double &mn, double &mx, double ratio);

static bool finiteValue(double value) {
    return std::isfinite(value);
}

static bool finitePositive(double value) {
    return finiteValue(value) && value > 0.0;
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

static double wrapPi(double x) {
    x = wrapTwoPi(x);
    if (x > PI) x -= 2.0 * PI;
    return x;
}

static double wrapSecondsOfDay(double sec) {
    sec = fmod(sec, static_cast<double>(SEC_PER_DAY));
    if (sec < 0.0) sec += SEC_PER_DAY;
    return sec;
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

static double secondsOfDay(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    return civ.hour * 3600.0 + civ.minute * 60.0 + civ.second;
}

static Geodetic ecefToGeodetic(const Eigen::Vector3d &xyz) {
    WGS84 wgs84;
    BLH blh = xyz2blh(XYZ(xyz), wgs84);

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

static Eigen::Vector3d rotateEarthVector(const Eigen::Vector3d &vec, double wt) {
    Eigen::Vector3d out = vec;
    out.x() = cos(wt) * vec.x() + sin(wt) * vec.y();
    out.y() = -sin(wt) * vec.x() + cos(wt) * vec.y();
    out.z() = vec.z();
    return out;
}

static Xvt rotateSatelliteToReceive(const Eigen::Vector3d &recXYZ,
                                    const Xvt &xvtTx) {
    const double travelSec = (xvtTx.x - recXYZ).norm() / C_MPS;
    const double wt = OMEGA_EARTH * travelSec;

    Xvt xvtRx = xvtTx;
    xvtRx.x = rotateEarthVector(xvtTx.x, wt);
    xvtRx.v = rotateEarthVector(xvtTx.v, wt);
    return xvtRx;
}

static KlobucharCoeffs findGpsKlobucharCoeffs(const NavHeaderAll &navHeader) {
    KlobucharCoeffs coeffs;
    bool hasAlpha = false;
    bool hasBeta = false;

    for (const auto &rec : navHeader.ionoCorrs) {
        if (rec.type == "GPSA" && rec.coeffs.size() >= 4) {
            copy_n(rec.coeffs.begin(), 4, coeffs.alpha.begin());
            hasAlpha = true;
        } else if (rec.type == "GPSB" && rec.coeffs.size() >= 4) {
            copy_n(rec.coeffs.begin(), 4, coeffs.beta.begin());
            hasBeta = true;
        }
    }

    if (!hasAlpha || !hasBeta) {
        throw runtime_error("GPSA/GPSB ionospheric coefficients were not found in the nav header.");
    }
    return coeffs;
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

static void addIfFiniteCode(ObsData &obsData,
                            const SatID &sat,
                            const map<string, double> &tv,
                            const string &obsType) {
    auto it = tv.find(obsType);
    if (it != tv.end() && finitePositive(it->second)) {
        obsData.satTypeValueData[sat][obsType] = it->second;
    }
}

static TransmitSolution computeTransmitSolution(RinexNavStore &navStore,
                                                const SatID &sat,
                                                const CommonTime &recvGps,
                                                double pseudorangeMeters) {
    TransmitSolution sol;
    sol.txGps = recvGps;
    sol.txGps -= pseudorangeMeters / C_MPS;

    for (int i = 0; i < 2; ++i) {
        sol.xvtTx = navStore.getXvt(sat, sol.txGps);
        sol.txGps = recvGps;
        sol.txGps -= pseudorangeMeters / C_MPS;
        sol.txGps -= (sol.xvtTx.clkbias + sol.xvtTx.relcorr);
    }

    sol.recvMinusTxSec = recvGps - sol.txGps;
    return sol;
}

static KlobucharResult computeKlobuchar(const Geodetic &recBLH,
                                        double azRad,
                                        double elevRad,
                                        double gpsSecondsOfDay,
                                        const KlobucharCoeffs &coeffs) {
    const double geomagPoleLat = 78.3 * DEG_TO_RAD;
    const double geomagPoleLon = 283.0 * DEG_TO_RAD;

    KlobucharResult out;
    const double elevSemiCircle = elevRad / PI;
    const double psi = (0.0137 / (elevSemiCircle + 0.11) - 0.022) * PI;

    const double ippLatRad = asin(clampUnit(sin(recBLH.latRad) * cos(psi)
                         + cos(recBLH.latRad) * sin(psi) * cos(azRad)));
    const double ippLonRad = wrapPi(recBLH.lonRad + psi * sin(azRad) / cos(ippLatRad));
    const double ippGeomagLatRad = asin(clampUnit(sin(ippLatRad) * sin(geomagPoleLat)
                                  + cos(ippLatRad) * cos(geomagPoleLat)
                                  * cos(ippLonRad - geomagPoleLon)));

    out.localTimeSec = wrapSecondsOfDay(43200.0 * ippLonRad / PI + gpsSecondsOfDay);

    const double phi = ippGeomagLatRad / PI;
    double powPhi = 1.0;
    for (int i = 0; i < 4; ++i) {
        out.amplitudeSec += coeffs.alpha[i] * powPhi;
        out.periodSec += coeffs.beta[i] * powPhi;
        powPhi *= phi;
    }
    if (out.amplitudeSec < 0.0) out.amplitudeSec = 0.0;
    if (out.periodSec < 72000.0) out.periodSec = 72000.0;

    const double phase = 2.0 * PI * (out.localTimeSec - 50400.0) / out.periodSec;
    out.mapping = 1.0 + 16.0 * pow(0.53 - elevSemiCircle, 3.0);

    if (abs(phase) < PI / 2.0) {
        out.delaySec = (5.0e-9 + out.amplitudeSec * cos(phase)) * out.mapping;
    } else {
        out.delaySec = 5.0e-9 * out.mapping;
    }

    out.delayMeters = out.delaySec * C_MPS;
    return out;
}

static TropResult computeSaastamoinen(const Geodetic &recBLH, double elevRad) {
    TropResult out;
    const double heightKm = recBLH.hMeters / 1000.0;

    out.pressureHpa = 1013.25 * pow(1.0 - 0.0000226 * recBLH.hMeters, 5.225);
    out.temperatureKelvin = 15.0 - 0.0065 * recBLH.hMeters + 273.15;
    out.saturationVaporHpa = 6.108 * exp((17.15 * out.temperatureKelvin - 4684.0)
                                         / (out.temperatureKelvin - 38.45));
    out.vaporPressureHpa = kRelativeHumidity * out.saturationVaporHpa;

    const double denom = 1.0 - 0.266e-2 * cos(2.0 * recBLH.latRad) - 0.28e-3 * heightKm;
    out.zhdMeters = 2.277e-3 * out.pressureHpa / denom;
    out.zwdMeters = 0.002277 * (1255.0 / out.temperatureKelvin + 0.05) * out.vaporPressureHpa;
    out.ztdMeters = out.zhdMeters + out.zwdMeters;

    double sinElev = sin(elevRad);
    if (sinElev < 0.05) {
        sinElev = 0.05;
    }
    out.mappingFunction = 1.0 / sinElev;
    out.slantDelayMeters = out.ztdMeters * out.mappingFunction;
    return out;
}

static double dualFreqL1IonoMeters(double p1Meters, double p2Meters) {
    const double gamma = (L1_FREQ_GPS * L1_FREQ_GPS) / (L2_FREQ_GPS * L2_FREQ_GPS);
    return (p2Meters - p1Meters) / (gamma - 1.0);
}

static double dualFreqL1PhaseIonoMeters(double l1Cycles, double l2Cycles) {
    const double gamma = (L1_FREQ_GPS * L1_FREQ_GPS) / (L2_FREQ_GPS * L2_FREQ_GPS);
    const double l1Meters = l1Cycles * L1_WAVELENGTH_GPS;
    const double l2Meters = l2Cycles * L2_WAVELENGTH_GPS;
    return (l1Meters - l2Meters) / (gamma - 1.0);
}

static double median(vector<double> values) {
    if (values.empty()) return 0.0;
    const size_t mid = values.size() / 2;
    nth_element(values.begin(), values.begin() + mid, values.end());
    double med = values[mid];
    if (values.size() % 2 == 0) {
        nth_element(values.begin(), values.begin() + mid - 1, values.end());
        med = 0.5 * (med + values[mid - 1]);
    }
    return med;
}

static void applyDualIonoLeveling(map<SatID, vector<ObservationRow>> &rowsBySat) {
    auto flushArc = [](vector<size_t> &arc, vector<ObservationRow> &rows) {
        if (arc.size() < kMinDualArcPoints) {
            arc.clear();
            return;
        }

        vector<double> leveling;
        leveling.reserve(arc.size());
        for (size_t idx : arc) {
            const auto &row = rows[idx];
            leveling.push_back(row.codeIonoMeters - row.phaseIonoMeters);
        }
        const double bias = median(leveling);

        for (size_t idx : arc) {
            auto &row = rows[idx];
            const double dual = row.phaseIonoMeters + bias;
            const double diff = dual - row.modelIonoMeters;
            if (!finiteValue(diff) || abs(diff) > kMaxDualModelDiffMeters) {
                continue;
            }
            row.hasDualLeveled = true;
            row.dualLeveledIonoMeters = dual;
            row.dualMinusModelMeters = diff;
        }
        arc.clear();
    };

    for (auto &satRows : rowsBySat) {
        auto &rows = satRows.second;
        vector<size_t> arc;
        CommonTime prevEpoch;
        double prevPhaseIono = 0.0;
        bool hasPrev = false;

        for (size_t idx = 0; idx < rows.size(); ++idx) {
            auto &row = rows[idx];
            const bool valid = row.hasDualRaw;
            if (!valid) {
                flushArc(arc, rows);
                hasPrev = false;
                continue;
            }

            bool newArc = false;
            if (hasPrev) {
                const double dt = row.epochGps - prevEpoch;
                if (dt > 120.0 || abs(row.phaseIonoMeters - prevPhaseIono) > 5.0) {
                    newArc = true;
                }
            }
            if (newArc) {
                flushArc(arc, rows);
            }

            arc.push_back(idx);
            prevEpoch = row.epochGps;
            prevPhaseIono = row.phaseIonoMeters;
            hasPrev = true;
        }

        flushArc(arc, rows);
    }
}

static vector<ObservationRow> flattenRows(const map<SatID, vector<ObservationRow>> &rowsBySat) {
    vector<ObservationRow> rows;
    for (const auto &entry : rowsBySat) {
        rows.insert(rows.end(), entry.second.begin(), entry.second.end());
    }
    sort(rows.begin(), rows.end(), [](const ObservationRow &a, const ObservationRow &b) {
        const double dt = a.epochGps - b.epochGps;
        if (abs(dt) > 1.0e-9) return dt < 0.0;
        return a.sat < b.sat;
    });
    return rows;
}

static vector<EpochSummary> summarizeByEpoch(const vector<ObservationRow> &rows) {
    struct Accumulator {
        size_t satCount{0};
        size_t dualCount{0};
        double elevationSum{0.0};
        double ionoModelSum{0.0};
        double dualIonoSum{0.0};
        double tropSum{0.0};
        double clkBiasNsSum{0.0};
        double relativityNsSum{0.0};
        double tgdNsSum{0.0};
    };

    map<CommonTime, Accumulator> acc;
    for (const auto &row : rows) {
        auto &item = acc[row.epochGps];
        ++item.satCount;
        item.elevationSum += row.elevationDeg;
        item.ionoModelSum += row.modelIonoMeters;
        item.tropSum += row.modelTropMeters;
        item.clkBiasNsSum += row.satClockBiasSec * 1.0e9;
        item.relativityNsSum += row.relativitySec * 1.0e9;
        item.tgdNsSum += row.tgdSec * 1.0e9;
        if (row.hasDualLeveled) {
            ++item.dualCount;
            item.dualIonoSum += row.dualLeveledIonoMeters;
        }
    }

    vector<EpochSummary> out;
    out.reserve(acc.size());
    for (const auto &entry : acc) {
        EpochSummary row;
        row.epochGps = entry.first;
        row.satCount = entry.second.satCount;
        row.dualCount = entry.second.dualCount;
        row.meanElevationDeg = entry.second.elevationSum / static_cast<double>(entry.second.satCount);
        row.meanIonoModelMeters = entry.second.ionoModelSum / static_cast<double>(entry.second.satCount);
        row.meanTropMeters = entry.second.tropSum / static_cast<double>(entry.second.satCount);
        row.meanClockBiasNs = entry.second.clkBiasNsSum / static_cast<double>(entry.second.satCount);
        row.meanRelativityNs = entry.second.relativityNsSum / static_cast<double>(entry.second.satCount);
        row.meanTgdNs = entry.second.tgdNsSum / static_cast<double>(entry.second.satCount);
        if (entry.second.dualCount > 0) {
            row.meanDualIonoMeters = entry.second.dualIonoSum / static_cast<double>(entry.second.dualCount);
        }
        out.push_back(row);
    }
    return out;
}

static vector<ObservationRow> filterRowsBySat(const vector<ObservationRow> &rows, const SatID &sat) {
    vector<ObservationRow> out;
    for (const auto &row : rows) {
        if (row.sat == sat) {
            out.push_back(row);
        }
    }
    return out;
}

static pair<double, double> computeSeriesRange(const vector<SvgSeries> &seriesList) {
    double mn = numeric_limits<double>::infinity();
    double mx = -numeric_limits<double>::infinity();
    for (const auto &series : seriesList) {
        for (const auto &pt : series.points) {
            if (!finiteValue(pt.second)) continue;
            mn = min(mn, pt.second);
            mx = max(mx, pt.second);
        }
    }
    padRange(mn, mx, 0.12);
    return {mn, mx};
}

static MetricStats computeMetricStats(const vector<ObservationRow> &rows,
                                     const function<bool(const ObservationRow &)> &predicate,
                                     const function<double(const ObservationRow &)> &getter) {
    MetricStats stats;
    bool first = true;
    for (const auto &row : rows) {
        if (!predicate(row)) continue;
        const double value = getter(row);
        if (!finiteValue(value)) continue;
        if (first) {
            stats.minVal = stats.maxVal = value;
            stats.minEpoch = stats.maxEpoch = row.epochGps;
            first = false;
        }
        ++stats.n;
        stats.mean += value;
        if (value < stats.minVal) {
            stats.minVal = value;
            stats.minEpoch = row.epochGps;
        }
        if (value > stats.maxVal) {
            stats.maxVal = value;
            stats.maxEpoch = row.epochGps;
        }
    }
    if (stats.n > 0) {
        stats.mean /= static_cast<double>(stats.n);
    }
    return stats;
}

static void padRange(double &mn, double &mx, double ratio) {
    if (!finiteValue(mn) || !finiteValue(mx)) {
        mn = -1.0;
        mx = 1.0;
        return;
    }
    if (mx <= mn) {
        const double delta = max(1.0, abs(mx) * ratio);
        mn -= delta;
        mx += delta;
        return;
    }
    const double span = mx - mn;
    mn -= span * ratio;
    mx += span * ratio;
}

static string svgPolyline(const SvgSeries &series,
                          double xMin,
                          double xMax,
                          double yMin,
                          double yMax,
                          double x0,
                          double y0,
                          double width,
                          double height) {
    if (series.points.empty()) return "";

    auto mapX = [&](double x) {
        if (xMax <= xMin) return x0;
        return x0 + (x - xMin) / (xMax - xMin) * width;
    };
    auto mapY = [&](double y) {
        if (yMax <= yMin) return y0 + height * 0.5;
        return y0 + height - (y - yMin) / (yMax - yMin) * height;
    };

    ostringstream oss;
    if (series.lineWidth > 0.0) {
        oss << "<polyline fill=\"none\" stroke=\"" << series.color
            << "\" stroke-width=\"" << series.lineWidth
            << "\" opacity=\"" << series.opacity
            << "\" points=\"";
        for (const auto &pt : series.points) {
            oss << mapX(pt.first) << "," << mapY(pt.second) << " ";
        }
        oss << "\"/>\n";
    }
    if (series.pointRadius > 0.0) {
        for (const auto &pt : series.points) {
            oss << "<circle cx=\"" << mapX(pt.first)
                << "\" cy=\"" << mapY(pt.second)
                << "\" r=\"" << series.pointRadius
                << "\" fill=\"" << (series.fillColor.empty() ? series.color : series.fillColor)
                << "\" opacity=\"" << series.opacity << "\"/>\n";
        }
    }
    return oss.str();
}

static void writeDetailCsv(const string &csvFile, const vector<ObservationRow> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open csv file: " + csvFile);

    csv << "epoch,sat,primary_code_type,primary_code_raw_m,primary_code_tgd_corrected_m,"
        << "tgd_correction_m,c1c_m,c1w_m,c2w_m,l1w_cycles,l2w_cycles,tx_offset_s,"
        << "sat_x_tx_m,sat_y_tx_m,sat_z_tx_m,sat_vx_tx_m,sat_vy_tx_m,sat_vz_tx_m,"
        << "sat_x_rx_m,sat_y_rx_m,sat_z_rx_m,sat_vx_rx_m,sat_vy_rx_m,sat_vz_rx_m,"
        << "sv_clock_bias_s,relativity_s,tgd_s,elevation_deg,azimuth_deg,"
        << "model_iono_m,model_trop_m,code_iono_m,phase_iono_m,dual_iono_m,dual_minus_model_m\n";

    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.sat << ","
            << row.primaryCodeType << ","
            << fixed << setprecision(3)
            << row.primaryCodeRawMeters << ","
            << row.primaryCodeTgdCorrectedMeters << ","
            << row.tgdCorrectionMeters << ","
            << row.c1cMeters << ","
            << row.c1wMeters << ","
            << row.c2wMeters << ","
            << setprecision(6)
            << row.l1wCycles << ","
            << row.l2wCycles << ","
            << scientific << setprecision(10)
            << row.transmitOffsetSec << ","
            << fixed << setprecision(3)
            << row.satPosTx.x() << ","
            << row.satPosTx.y() << ","
            << row.satPosTx.z() << ","
            << row.satVelTx.x() << ","
            << row.satVelTx.y() << ","
            << row.satVelTx.z() << ","
            << row.satPosRx.x() << ","
            << row.satPosRx.y() << ","
            << row.satPosRx.z() << ","
            << row.satVelRx.x() << ","
            << row.satVelRx.y() << ","
            << row.satVelRx.z() << ","
            << scientific << setprecision(10)
            << row.satClockBiasSec << ","
            << row.relativitySec << ","
            << row.tgdSec << ","
            << fixed << setprecision(6)
            << row.elevationDeg << ","
            << row.azimuthDeg << ","
            << row.modelIonoMeters << ","
            << row.modelTropMeters << ","
            << row.codeIonoMeters << ","
            << row.phaseIonoMeters << ","
            << row.dualLeveledIonoMeters << ","
            << row.dualMinusModelMeters << "\n";
    }
}

static void writeSummaryCsv(const string &csvFile, const vector<EpochSummary> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open summary csv file: " + csvFile);

    csv << "epoch,sat_count,dual_count,mean_elevation_deg,mean_iono_model_m,"
        << "mean_dual_iono_m,mean_trop_m,mean_clock_bias_ns,mean_relativity_ns,mean_tgd_ns\n";
    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.satCount << ","
            << row.dualCount << ","
            << fixed << setprecision(6)
            << row.meanElevationDeg << ","
            << row.meanIonoModelMeters << ","
            << row.meanDualIonoMeters << ","
            << row.meanTropMeters << ","
            << row.meanClockBiasNs << ","
            << row.meanRelativityNs << ","
            << row.meanTgdNs << "\n";
    }
}

static void writeOverviewSvg(const string &svgFile,
                             const vector<ObservationRow> &rows,
                             const vector<EpochSummary> &summary,
                             const MetricStats &ionoStats,
                             const MetricStats &tropStats,
                             const MetricStats &dualStats) {
    ofstream svg(svgFile);
    if (!svg) throw runtime_error("failed to open svg file: " + svgFile);

    vector<pair<double, double>> ionoMean;
    vector<pair<double, double>> dualMean;
    vector<pair<double, double>> tropMean;

    double ionoMin = numeric_limits<double>::infinity();
    double ionoMax = -numeric_limits<double>::infinity();
    double tropMin = numeric_limits<double>::infinity();
    double tropMax = -numeric_limits<double>::infinity();

    for (const auto &row : summary) {
        const double x = secondsOfDay(row.epochGps);
        ionoMean.emplace_back(x, row.meanIonoModelMeters);
        tropMean.emplace_back(x, row.meanTropMeters);
        ionoMin = min(ionoMin, row.meanIonoModelMeters);
        ionoMax = max(ionoMax, row.meanIonoModelMeters);
        tropMin = min(tropMin, row.meanTropMeters);
        tropMax = max(tropMax, row.meanTropMeters);
        if (finiteValue(row.meanDualIonoMeters)) {
            dualMean.emplace_back(x, row.meanDualIonoMeters);
            ionoMin = min(ionoMin, row.meanDualIonoMeters);
            ionoMax = max(ionoMax, row.meanDualIonoMeters);
        }
    }

    padRange(ionoMin, ionoMax, 0.12);
    padRange(tropMin, tropMax, 0.12);

    const int width = 1380;
    const int height = 980;
    const int marginL = 78;
    const int marginR = 34;
    const int marginT = 86;
    const int marginB = 56;
    const int panelGap = 86;
    const int panelH = (height - marginT - marginB - panelGap) / 2;
    const int plotW = width - marginL - marginR;
    const int topY = marginT;
    const int bottomY = marginT + panelH + panelGap;
    const double xMin = 0.0;
    const double xMax = 86400.0;

    SvgSeries ionoMeanSeries{"Model ionosphere mean", "#1d4ed8", "", 2.4, 0.0, 1.0, ionoMean};
    SvgSeries dualMeanSeries{"Dual-frequency leveled mean", "#c2410c", "", 2.4, 0.0, 1.0, dualMean};
    SvgSeries tropMeanSeries{"Troposphere mean", "#15803d", "", 2.4, 0.0, 1.0, tropMean};

    auto drawAxes = [&](ostringstream &oss, double x0, double y0, double yMin, double yMax, const string &yLabel) {
        const double y1 = y0 + panelH;
        oss << "<rect x=\"" << x0 << "\" y=\"" << y0 << "\" width=\"" << plotW
            << "\" height=\"" << panelH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";
        for (int i = 0; i <= 6; ++i) {
            const double frac = static_cast<double>(i) / 6.0;
            const double y = y1 - frac * panelH;
            const double val = yMin + frac * (yMax - yMin);
            oss << "<line x1=\"" << x0 << "\" y1=\"" << y << "\" x2=\"" << (x0 + plotW)
                << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\"/>\n";
            oss << "<text x=\"" << (x0 - 10) << "\" y=\"" << (y + 4)
                << "\" text-anchor=\"end\" class=\"tick\">" << fixed << setprecision(1)
                << val << "</text>\n";
        }
        for (int i = 0; i <= 8; ++i) {
            const double frac = static_cast<double>(i) / 8.0;
            const double x = x0 + frac * plotW;
            const int totalMinutes = static_cast<int>((xMin + frac * (xMax - xMin)) / 60.0);
            const int hh = (totalMinutes / 60) % 24;
            const int mm = totalMinutes % 60;
            ostringstream hhmm;
            hhmm << setfill('0') << setw(2) << hh << ":" << setw(2) << mm;
            oss << "<line x1=\"" << x << "\" y1=\"" << y0 << "\" x2=\"" << x
                << "\" y2=\"" << y1 << "\" stroke=\"#eef2f7\"/>\n";
            oss << "<text x=\"" << x << "\" y=\"" << (y1 + 22)
                << "\" text-anchor=\"middle\" class=\"tick\">" << hhmm.str() << "</text>\n";
        }
        oss << "<line x1=\"" << x0 << "\" y1=\"" << y1 << "\" x2=\"" << (x0 + plotW)
            << "\" y2=\"" << y1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
        oss << "<line x1=\"" << x0 << "\" y1=\"" << y0 << "\" x2=\"" << x0
            << "\" y2=\"" << y1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
        oss << "<text x=\"" << (x0 - 58) << "\" y=\"" << (y0 + panelH / 2)
            << "\" transform=\"rotate(-90 " << (x0 - 58) << " " << (y0 + panelH / 2)
            << ")\" text-anchor=\"middle\" class=\"axis\">" << yLabel << "</text>\n";
        oss << "<text x=\"" << (x0 + plotW / 2) << "\" y=\"" << (y1 + 42)
            << "\" text-anchor=\"middle\" class=\"axis\">time of day (GPST)</text>\n";
    };

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    svg << "<style>"
        << ".title{font:700 24px 'Arial','Microsoft YaHei',sans-serif;fill:#111827;}"
        << ".subtitle{font:13px 'Arial','Microsoft YaHei',sans-serif;fill:#4b5563;}"
        << ".panel{font:700 16px 'Arial','Microsoft YaHei',sans-serif;fill:#1f2937;}"
        << ".axis{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << ".tick{font:12px 'Arial','Microsoft YaHei',sans-serif;fill:#6b7280;}"
        << ".legend{font:13px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << "</style>\n";

    svg << "<text x=\"" << marginL << "\" y=\"36\" class=\"title\">Exercise 5.8: satellite state, ionosphere and troposphere over one GPS day</text>\n";
    svg << "<text x=\"" << marginL << "\" y=\"58\" class=\"subtitle\">"
        << "Rows=" << rows.size()
        << ", ionosphere mean=" << fixed << setprecision(3) << ionoStats.mean << " m"
        << ", troposphere mean=" << tropStats.mean << " m";
    if (dualStats.n > 0) {
        svg << ", dual-minus-model mean=" << dualStats.mean << " m";
    }
    svg << "</text>\n";

    ostringstream body;
    drawAxes(body, marginL, topY, ionoMin, ionoMax, "ionospheric delay / m");
    drawAxes(body, marginL, bottomY, tropMin, tropMax, "tropospheric delay / m");

    body << "<text x=\"" << (marginL + 12) << "\" y=\"" << (topY + 24)
         << "\" class=\"panel\">GPS L1 ionosphere: epoch-mean model and dual-frequency leveled mean</text>\n";
    body << "<text x=\"" << (marginL + 12) << "\" y=\"" << (bottomY + 24)
         << "\" class=\"panel\">Saastamoinen slant troposphere: epoch mean</text>\n";

    body << svgPolyline(ionoMeanSeries, xMin, xMax, ionoMin, ionoMax, marginL, topY, plotW, panelH);
    body << svgPolyline(dualMeanSeries, xMin, xMax, ionoMin, ionoMax, marginL, topY, plotW, panelH);

    body << svgPolyline(tropMeanSeries, xMin, xMax, tropMin, tropMax, marginL, bottomY, plotW, panelH);

    const double legendX = marginL + plotW - 330;
    const double legendY = topY + 20;
    const vector<pair<string, string>> legendItems = {
        {"#1d4ed8", "Model ionosphere mean"},
        {"#c2410c", "Dual-frequency leveled mean"},
        {"#15803d", "Troposphere mean"}
    };
    for (size_t i = 0; i < legendItems.size(); ++i) {
        const double y = legendY + i * 20.0;
        body << "<line x1=\"" << legendX << "\" y1=\"" << y << "\" x2=\"" << (legendX + 22)
             << "\" y2=\"" << y << "\" stroke=\"" << legendItems[i].first
             << "\" stroke-width=\"3\"/>\n";
        body << "<text x=\"" << (legendX + 30) << "\" y=\"" << (y + 4)
             << "\" class=\"legend\">" << legendItems[i].second << "</text>\n";
    }

    svg << body.str();
    svg << "</svg>\n";
}

#ifdef _WIN32
static int getPngEncoderClsid(CLSID *pClsid) {
    UINT num = 0;
    UINT size = 0;
    if (GetImageEncodersSize(&num, &size) != Ok || size == 0) return -1;
    vector<BYTE> buffer(size);
    auto *pImageCodecInfo = reinterpret_cast<ImageCodecInfo *>(buffer.data());
    if (GetImageEncoders(num, size, pImageCodecInfo) != Ok) return -1;
    for (UINT j = 0; j < num; ++j) {
        if (wcscmp(pImageCodecInfo[j].MimeType, L"image/png") == 0) {
            *pClsid = pImageCodecInfo[j].Clsid;
            return static_cast<int>(j);
        }
    }
    return -1;
}

struct GdiplusSession {
    ULONG_PTR token{0};
    GdiplusSession() {
        GdiplusStartupInput input;
        input.GdiplusVersion = 1;
        if (GdiplusStartup(&token, &input, nullptr) != Ok) {
            throw runtime_error("GDI+ startup failed.");
        }
    }
    ~GdiplusSession() {
        if (token != 0) {
            GdiplusShutdown(token);
        }
    }
};

static wstring widenPath(const string &path) {
    if (path.empty()) return wstring();
    int sizeNeeded = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, nullptr, 0);
    if (sizeNeeded <= 0) {
        sizeNeeded = MultiByteToWideChar(CP_ACP, 0, path.c_str(), -1, nullptr, 0);
    }
    if (sizeNeeded <= 0) {
        return wstring(path.begin(), path.end());
    }
    wstring wide(static_cast<size_t>(sizeNeeded - 1), L'\0');
    if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, &wide[0], sizeNeeded) <= 0) {
        MultiByteToWideChar(CP_ACP, 0, path.c_str(), -1, &wide[0], sizeNeeded);
    }
    return wide;
}

static wstring widenText(const string &text) {
    return widenPath(text);
}

static Color parseColor(const string &hex, double opacity = 1.0) {
    if (hex.size() != 7 || hex[0] != '#') {
        return Color(static_cast<BYTE>(255 * opacity), 0, 0, 0);
    }
    auto fromHex = [](char ch) -> int {
        if (ch >= '0' && ch <= '9') return ch - '0';
        if (ch >= 'a' && ch <= 'f') return 10 + ch - 'a';
        if (ch >= 'A' && ch <= 'F') return 10 + ch - 'A';
        return 0;
    };
    const int r = fromHex(hex[1]) * 16 + fromHex(hex[2]);
    const int g = fromHex(hex[3]) * 16 + fromHex(hex[4]);
    const int b = fromHex(hex[5]) * 16 + fromHex(hex[6]);
    return Color(static_cast<BYTE>(max(0.0, min(255.0, opacity * 255.0))), r, g, b);
}

static void drawSeries(Graphics &graphics,
                       const SvgSeries &series,
                       double xMin,
                       double xMax,
                       double yMin,
                       double yMax,
                       REAL x0,
                       REAL y0,
                       REAL width,
                       REAL height) {
    auto mapX = [&](double x) -> REAL {
        if (xMax <= xMin) return x0;
        return x0 + static_cast<REAL>((x - xMin) / (xMax - xMin) * width);
    };
    auto mapY = [&](double y) -> REAL {
        if (yMax <= yMin) return y0 + height * 0.5f;
        return y0 + height - static_cast<REAL>((y - yMin) / (yMax - yMin) * height);
    };

    const Color color = parseColor(series.color, series.opacity);
    Pen pen(color, static_cast<REAL>(max(1.0, series.lineWidth)));
    pen.SetLineJoin(LineJoinRound);
    pen.SetStartCap(LineCapRound);
    pen.SetEndCap(LineCapRound);

    vector<PointF> segment;
    const double maxGapSec = 180.0;
    for (size_t i = 0; i < series.points.size(); ++i) {
        const auto &pt = series.points[i];
        if (!finiteValue(pt.second)) {
            if (segment.size() >= 2 && series.lineWidth > 0.0) {
                graphics.DrawLines(&pen, segment.data(), static_cast<INT>(segment.size()));
            }
            segment.clear();
            continue;
        }
        if (!segment.empty()) {
            const double gap = pt.first - series.points[i - 1].first;
            if (gap > maxGapSec) {
                if (segment.size() >= 2 && series.lineWidth > 0.0) {
                    graphics.DrawLines(&pen, segment.data(), static_cast<INT>(segment.size()));
                }
                segment.clear();
            }
        }
        segment.emplace_back(mapX(pt.first), mapY(pt.second));
    }
    if (segment.size() >= 2 && series.lineWidth > 0.0) {
        graphics.DrawLines(&pen, segment.data(), static_cast<INT>(segment.size()));
    }

    if (series.pointRadius > 0.0) {
        SolidBrush brush(parseColor(series.fillColor.empty() ? series.color : series.fillColor, series.opacity));
        const REAL r = static_cast<REAL>(series.pointRadius);
        for (const auto &pt : series.points) {
            if (!finiteValue(pt.second)) continue;
            const REAL cx = mapX(pt.first);
            const REAL cy = mapY(pt.second);
            graphics.FillEllipse(&brush, cx - r, cy - r, 2 * r, 2 * r);
        }
    }
}

static void writePanelsPng(const string &pngFile,
                           const string &title,
                           const string &subtitle,
                           const vector<ChartPanel> &panels) {
    if (panels.empty()) return;

    GdiplusSession gdiplusSession;
    CLSID pngClsid;
    if (getPngEncoderClsid(&pngClsid) < 0) {
        throw runtime_error("PNG encoder not found in GDI+.");
    }

    const int width = 1520;
    const int height = 280 + static_cast<int>(panels.size()) * 420;
    const REAL marginL = 118.0f;
    const REAL marginR = 48.0f;
    const REAL marginT = 110.0f;
    const REAL marginB = 70.0f;
    const REAL panelGap = 72.0f;
    const REAL plotW = static_cast<REAL>(width) - marginL - marginR;
    const REAL panelH = (static_cast<REAL>(height) - marginT - marginB
                         - panelGap * static_cast<REAL>(panels.size() - 1))
                        / static_cast<REAL>(panels.size());

    Bitmap bitmap(width, height, PixelFormat32bppARGB);
    Graphics graphics(&bitmap);
    graphics.SetSmoothingMode(SmoothingModeHighQuality);
    graphics.Clear(Color(255, 255, 255, 255));

    SolidBrush textBrush(Color(255, 17, 24, 39));
    SolidBrush subBrush(Color(255, 75, 85, 99));
    Pen gridPen(Color(255, 229, 231, 235), 1.0f);
    Pen axisPen(Color(255, 71, 85, 105), 1.5f);
    Pen framePen(Color(255, 203, 213, 225), 1.0f);

    Font titleFont(L"Arial", 24.0f, FontStyleBold, UnitPixel);
    Font panelFont(L"Arial", 16.0f, FontStyleBold, UnitPixel);
    Font axisFont(L"Arial", 14.0f, FontStyleRegular, UnitPixel);
    Font tickFont(L"Arial", 12.0f, FontStyleRegular, UnitPixel);
    Font legendFont(L"Arial", 13.0f, FontStyleRegular, UnitPixel);

    graphics.DrawString(widenText(title).c_str(), -1, &titleFont,
                        PointF(marginL, 24.0f), &textBrush);
    graphics.DrawString(widenText(subtitle).c_str(), -1, &axisFont,
                        PointF(marginL, 58.0f), &subBrush);

    const double xMin = 0.0;
    const double xMax = 86400.0;

    for (size_t p = 0; p < panels.size(); ++p) {
        const REAL panelY = marginT + static_cast<REAL>(p) * (panelH + panelGap);
        const REAL plotY = panelY + 26.0f;
        const REAL innerH = panelH - 52.0f;

        graphics.DrawRectangle(&framePen, marginL, plotY, plotW, innerH);
        graphics.DrawString(widenText(panels[p].title).c_str(), -1, &panelFont,
                            PointF(marginL + 12.0f, panelY), &textBrush);

        for (int i = 0; i <= 6; ++i) {
            const double frac = static_cast<double>(i) / 6.0;
            const REAL y = plotY + innerH - static_cast<REAL>(frac * innerH);
            const double val = panels[p].yMin + frac * (panels[p].yMax - panels[p].yMin);
            graphics.DrawLine(&gridPen, marginL, y, marginL + plotW, y);

            ostringstream tick;
            tick << fixed << setprecision(1) << val;
            graphics.DrawString(widenText(tick.str()).c_str(), -1, &tickFont,
                                PointF(26.0f, y - 8.0f), &subBrush);
        }

        for (int i = 0; i <= 8; ++i) {
            const double frac = static_cast<double>(i) / 8.0;
            const REAL x = marginL + static_cast<REAL>(frac * plotW);
            graphics.DrawLine(&gridPen, x, plotY, x, plotY + innerH);

            const int totalMinutes = static_cast<int>((xMin + frac * (xMax - xMin)) / 60.0);
            const int hh = (totalMinutes / 60) % 24;
            const int mm = totalMinutes % 60;
            ostringstream hhmm;
            hhmm << setfill('0') << setw(2) << hh << ":" << setw(2) << mm;
            graphics.DrawString(widenText(hhmm.str()).c_str(), -1, &tickFont,
                                PointF(x - 18.0f, plotY + innerH + 8.0f), &subBrush);
        }

        graphics.DrawLine(&axisPen, marginL, plotY + innerH, marginL + plotW, plotY + innerH);
        graphics.DrawLine(&axisPen, marginL, plotY, marginL, plotY + innerH);

        for (const auto &series : panels[p].seriesList) {
            drawSeries(graphics, series, xMin, xMax, panels[p].yMin, panels[p].yMax,
                       marginL, plotY, plotW, innerH);
        }

        graphics.DrawString(widenText("time of day (GPST)").c_str(), -1, &axisFont,
                            PointF(marginL + plotW * 0.5f - 64.0f, plotY + innerH + 36.0f), &textBrush);

        GraphicsState state = graphics.Save();
        graphics.TranslateTransform(32.0f, plotY + innerH * 0.5f);
        graphics.RotateTransform(-90.0f);
        graphics.DrawString(widenText(panels[p].yLabel).c_str(), -1, &axisFont,
                            PointF(0.0f, 0.0f), &textBrush);
        graphics.Restore(state);

        const REAL legendX = marginL + plotW - 290.0f;
        const REAL legendY = panelY + 2.0f;
        for (size_t i = 0; i < panels[p].seriesList.size(); ++i) {
            const auto &series = panels[p].seriesList[i];
            const REAL y = legendY + static_cast<REAL>(i) * 18.0f;
            Pen pen(parseColor(series.color), 3.0f);
            graphics.DrawLine(&pen, legendX, y + 8.0f, legendX + 22.0f, y + 8.0f);
            graphics.DrawString(widenText(series.label).c_str(), -1, &legendFont,
                                PointF(legendX + 30.0f, y), &subBrush);
        }
    }

    const wstring wpath = widenPath(pngFile);
    if (bitmap.Save(wpath.c_str(), &pngClsid, nullptr) != Ok) {
        throw runtime_error("Failed to save PNG: " + pngFile);
    }
}
#endif

static string formatTimeOnly(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    ostringstream oss;
    oss << setw(2) << setfill('0') << civ.hour << ":"
        << setw(2) << civ.minute << ":"
        << setw(2) << static_cast<int>(round(civ.second));
    return oss.str();
}

static vector<ChartPanel> buildOverviewPanels(const vector<EpochSummary> &summary) {
    vector<pair<double, double>> ionoMean;
    vector<pair<double, double>> dualMean;
    vector<pair<double, double>> tropMean;

    for (const auto &row : summary) {
        const double x = secondsOfDay(row.epochGps);
        ionoMean.emplace_back(x, row.meanIonoModelMeters);
        tropMean.emplace_back(x, row.meanTropMeters);
        if (finiteValue(row.meanDualIonoMeters)) {
            dualMean.emplace_back(x, row.meanDualIonoMeters);
        }
    }

    ChartPanel ionoPanel;
    ionoPanel.title = "Epoch-mean ionosphere";
    ionoPanel.yLabel = "ionospheric delay / m";
    ionoPanel.seriesList = {
        SvgSeries{"Model ionosphere mean", "#1d4ed8", "", 2.4, 0.0, 1.0, ionoMean},
        SvgSeries{"Dual-frequency leveled mean", "#c2410c", "", 2.4, 0.0, 1.0, dualMean}
    };
    {
        auto range = computeSeriesRange(ionoPanel.seriesList);
        ionoPanel.yMin = range.first;
        ionoPanel.yMax = range.second;
    }

    ChartPanel tropPanel;
    tropPanel.title = "Epoch-mean troposphere";
    tropPanel.yLabel = "tropospheric delay / m";
    tropPanel.seriesList = {
        SvgSeries{"Troposphere mean", "#15803d", "", 2.4, 0.0, 1.0, tropMean}
    };
    {
        auto range = computeSeriesRange(tropPanel.seriesList);
        tropPanel.yMin = range.first;
        tropPanel.yMax = range.second;
    }

    return {ionoPanel, tropPanel};
}

static vector<ChartPanel> buildG01Panels(const vector<ObservationRow> &g01Rows) {
    vector<pair<double, double>> modelIono;
    vector<pair<double, double>> dualIono;
    vector<pair<double, double>> trop;

    for (const auto &row : g01Rows) {
        const double x = secondsOfDay(row.epochGps);
        modelIono.emplace_back(x, row.modelIonoMeters);
        trop.emplace_back(x, row.modelTropMeters);
        if (row.hasDualLeveled) {
            dualIono.emplace_back(x, row.dualLeveledIonoMeters);
        }
    }

    ChartPanel ionoPanel;
    ionoPanel.title = "G01 ionosphere";
    ionoPanel.yLabel = "delay / m";
    ionoPanel.seriesList = {
        SvgSeries{"Model ionosphere", "#1d4ed8", "", 2.2, 0.0, 1.0, modelIono},
        SvgSeries{"Dual-frequency leveled ionosphere", "#c2410c", "", 2.2, 0.0, 1.0, dualIono}
    };
    {
        auto range = computeSeriesRange(ionoPanel.seriesList);
        ionoPanel.yMin = range.first;
        ionoPanel.yMax = range.second;
    }

    ChartPanel tropPanel;
    tropPanel.title = "G01 troposphere";
    tropPanel.yLabel = "delay / m";
    tropPanel.seriesList = {
        SvgSeries{"Saastamoinen troposphere", "#15803d", "", 2.2, 0.0, 1.0, trop}
    };
    {
        auto range = computeSeriesRange(tropPanel.seriesList);
        tropPanel.yMin = range.first;
        tropPanel.yMax = range.second;
    }

    return {ionoPanel, tropPanel};
}

} // namespace

int main() {
    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\";
    const string obsFile = dataDir + "WUH200CHN_R_20250010000_01D_30S_MO.rnx";
    const string navFile = dataDir + "BRDC00IGS_R_20250010000_01D_MN.rnx";
    const string detailCsvFile = dataDir + "exam-5.8-tropo_iono_satpos.csv";
    const string summaryCsvFile = dataDir + "exam-5.8-tropo_iono_satpos_summary.csv";
    const string svgFile = dataDir + "exam-5.8-tropo_iono_satpos.svg";
    const string overviewPngFile = dataDir + "exam-5.8-tropo_iono_satpos.png";
    const string g01PngFile = dataDir + "exam-5.8-tropo_iono_satpos_g01.png";

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
    const KlobucharCoeffs ionoCoeffs = findGpsKlobucharCoeffs(navHeader);

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

    map<SatID, vector<ObservationRow>> rowsBySat;
    size_t totalEpochs = 0;
    size_t totalGpsSeen = 0;
    size_t totalAccepted = 0;
    size_t rejectedNoCode = 0;
    size_t rejectedNav = 0;
    size_t rejectedLowElev = 0;
    size_t dualRawCount = 0;

    cout << "Exercise 5.8: satellite position + system error corrections" << endl;
    cout << "Observation file: " << obsFile << endl;
    cout << "Navigation file : " << navFile << endl;
    cout << "Station         : " << obsHeader.markerName << endl;
    cout << "Receiver XYZ    : (" << fixed << setprecision(3)
         << recXYZ.x() << ", " << recXYZ.y() << ", " << recXYZ.z() << ") m" << endl;
    cout << "Receiver BLH    : lat=" << setprecision(8) << recBLH.latRad * RAD_TO_DEG
         << " deg, lon=" << recBLH.lonRad * RAD_TO_DEG
         << " deg, h=" << setprecision(3) << recBLH.hMeters << " m" << endl;
    cout << "GPS nav records : " << gpsRecordCount << endl;
    cout << "Elevation cutoff: " << kElevationCutoffDeg << " deg" << endl;
    cout << "Relative humidity used in Saastamoinen: " << kRelativeHumidity << endl << endl;

    while (true) {
        ObsEpochData epochData;
        try {
            epochData = obsReader.readNextEpoch();
        } catch (const EndOfObsFile &) {
            break;
        }

        ++totalEpochs;
        const CommonTime epochGps = toGpsTime(epochData.time);
        const double sod = secondsOfDay(epochGps);

        map<SatID, ObservationRow> epochRows;
        map<SatID, Xvt> satXvtTransTime;

        ObsData tgdObs;
        tgdObs.station = obsHeader.markerName;
        tgdObs.epoch = epochGps;
        tgdObs.antennaPosition = recXYZ;

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

            TransmitSolution tx;
            try {
                tx = computeTransmitSolution(navStore, sat, epochGps, tv.at(primaryCodeType));
            } catch (...) {
                ++rejectedNav;
                continue;
            }

            const Xvt xvtRx = rotateSatelliteToReceive(recXYZ, tx.xvtTx);
            const AzEl azel = computeAzEl(recXYZ, recBLH, xvtRx.x);
            const double elevDeg = azel.elRad * RAD_TO_DEG;
            if (elevDeg < kElevationCutoffDeg) {
                ++rejectedLowElev;
                continue;
            }

            ObservationRow row;
            row.epochGps = epochGps;
            row.sat = sat;
            row.primaryCodeType = primaryCodeType;
            row.primaryCodeRawMeters = tv.at(primaryCodeType);
            auto c1cIt = tv.find("C1C");
            auto c1wIt = tv.find("C1W");
            auto c2wIt = tv.find("C2W");
            auto l1wIt = tv.find("L1W");
            auto l2wIt = tv.find("L2W");
            row.c1cMeters = (c1cIt != tv.end()) ? c1cIt->second : numeric_limits<double>::quiet_NaN();
            row.c1wMeters = (c1wIt != tv.end()) ? c1wIt->second : numeric_limits<double>::quiet_NaN();
            row.c2wMeters = (c2wIt != tv.end()) ? c2wIt->second : numeric_limits<double>::quiet_NaN();
            row.l1wCycles = (l1wIt != tv.end()) ? l1wIt->second : numeric_limits<double>::quiet_NaN();
            row.l2wCycles = (l2wIt != tv.end()) ? l2wIt->second : numeric_limits<double>::quiet_NaN();
            row.satPosTx = tx.xvtTx.x;
            row.satVelTx = tx.xvtTx.v;
            row.satPosRx = xvtRx.x;
            row.satVelRx = xvtRx.v;
            row.transmitOffsetSec = tx.recvMinusTxSec;
            row.satClockBiasSec = tx.xvtTx.clkbias;
            row.relativitySec = tx.xvtTx.relcorr;
            auto tgdIt = tx.xvtTx.typeTGDData.find("TGD");
            if (tgdIt != tx.xvtTx.typeTGDData.end()) {
                row.tgdSec = tgdIt->second;
            }
            row.elevationDeg = elevDeg;
            row.azimuthDeg = azel.azRad * RAD_TO_DEG;
            row.modelIonoMeters = computeKlobuchar(recBLH, azel.azRad, azel.elRad, sod, ionoCoeffs).delayMeters;
            row.modelTropMeters = computeSaastamoinen(recBLH, azel.elRad).slantDelayMeters;

            if (finitePositive(row.c1wMeters) &&
                finitePositive(row.c2wMeters) &&
                finiteValue(row.l1wCycles) &&
                finiteValue(row.l2wCycles)) {
                row.codeIonoMeters = dualFreqL1IonoMeters(row.c1wMeters, row.c2wMeters);
                row.phaseIonoMeters = dualFreqL1PhaseIonoMeters(row.l1wCycles, row.l2wCycles);
                row.hasDualRaw = finiteValue(row.codeIonoMeters) &&
                                 finiteValue(row.phaseIonoMeters) &&
                                 abs(row.codeIonoMeters) < 100.0 &&
                                 abs(row.phaseIonoMeters) < 1.0e8;
                if (row.hasDualRaw) {
                    ++dualRawCount;
                }
            }

            satXvtTransTime[sat] = tx.xvtTx;
            addIfFiniteCode(tgdObs, sat, tv, "C1C");
            addIfFiniteCode(tgdObs, sat, tv, "C1W");
            addIfFiniteCode(tgdObs, sat, tv, "C2W");
            epochRows[sat] = row;
            ++totalAccepted;
        }

        if (!satXvtTransTime.empty() && !tgdObs.satTypeValueData.empty()) {
            vector<ObsDelayCorrectionRecord> tgdReport;
            correctTGD(tgdObs, satXvtTransTime, nullptr, &tgdReport);

            map<pair<SatID, string>, ObsDelayCorrectionRecord> reportMap;
            for (const auto &rec : tgdReport) {
                reportMap[{rec.sat, rec.obsType}] = rec;
            }

            for (auto &entry : epochRows) {
                auto &row = entry.second;
                auto satIt = tgdObs.satTypeValueData.find(entry.first);
                if (satIt != tgdObs.satTypeValueData.end()) {
                    auto codeIt = satIt->second.find(row.primaryCodeType);
                    if (codeIt != satIt->second.end() && finitePositive(codeIt->second)) {
                        row.primaryCodeTgdCorrectedMeters = codeIt->second;
                    } else {
                        row.primaryCodeTgdCorrectedMeters = row.primaryCodeRawMeters;
                    }
                } else {
                    row.primaryCodeTgdCorrectedMeters = row.primaryCodeRawMeters;
                }

                auto reportIt = reportMap.find({entry.first, row.primaryCodeType});
                if (reportIt != reportMap.end()) {
                    row.tgdCorrectionMeters = reportIt->second.correctionMeters;
                }
                rowsBySat[entry.first].push_back(row);
            }
        } else {
            for (auto &entry : epochRows) {
                entry.second.primaryCodeTgdCorrectedMeters = entry.second.primaryCodeRawMeters;
                rowsBySat[entry.first].push_back(entry.second);
            }
        }
    }

    applyDualIonoLeveling(rowsBySat);
    const vector<ObservationRow> rows = flattenRows(rowsBySat);
    const vector<EpochSummary> summary = summarizeByEpoch(rows);

    const MetricStats ionoStats = computeMetricStats(
        rows,
        [](const ObservationRow &) { return true; },
        [](const ObservationRow &row) { return row.modelIonoMeters; });

    const MetricStats tropStats = computeMetricStats(
        rows,
        [](const ObservationRow &) { return true; },
        [](const ObservationRow &row) { return row.modelTropMeters; });

    const MetricStats dualStats = computeMetricStats(
        rows,
        [](const ObservationRow &row) { return row.hasDualLeveled; },
        [](const ObservationRow &row) { return row.dualMinusModelMeters; });

    const vector<ChartPanel> overviewPanels = buildOverviewPanels(summary);
    const SatID g01("G01");
    const vector<ObservationRow> g01Rows = filterRowsBySat(rows, g01);
    const vector<ChartPanel> g01Panels = buildG01Panels(g01Rows);

    writeDetailCsv(detailCsvFile, rows);
    writeSummaryCsv(summaryCsvFile, summary);
    writeOverviewSvg(svgFile, rows, summary, ionoStats, tropStats, dualStats);
#ifdef _WIN32
    writePanelsPng(overviewPngFile,
                   "Exercise 5.8: epoch-mean ionosphere and troposphere",
                   "Mean curves only, derived from all accepted GPS samples.",
                   overviewPanels);
    if (!g01Rows.empty()) {
        writePanelsPng(g01PngFile,
                       "Exercise 5.8: G01 ionosphere and troposphere",
                       "Model ionosphere, dual-frequency leveled ionosphere, and Saastamoinen troposphere.",
                       g01Panels);
    }
#endif

    cout << "--- Daily processing summary ---" << endl;
    cout << "Epochs read                : " << totalEpochs << endl;
    cout << "GPS observations seen      : " << totalGpsSeen << endl;
    cout << "Accepted GPS samples       : " << rows.size() << endl;
    cout << "Rejected no primary code   : " << rejectedNoCode << endl;
    cout << "Rejected no nav solution   : " << rejectedNav << endl;
    cout << "Rejected low elevation     : " << rejectedLowElev << endl;
    cout << "Dual-frequency raw samples : " << dualRawCount << endl;
    cout << "Dual-frequency leveled     : "
         << count_if(rows.begin(), rows.end(), [](const ObservationRow &row) { return row.hasDualLeveled; })
         << endl;
    cout << endl;

    cout << fixed << setprecision(3);
    cout << "Model ionosphere mean/min/max : " << ionoStats.mean << " / " << ionoStats.minVal
         << " / " << ionoStats.maxVal << " m" << endl;
    cout << "  min epoch: " << formatCivilTime(ionoStats.minEpoch)
         << ", max epoch: " << formatCivilTime(ionoStats.maxEpoch) << endl;
    cout << "Troposphere mean/min/max      : " << tropStats.mean << " / " << tropStats.minVal
         << " / " << tropStats.maxVal << " m" << endl;
    cout << "  min epoch: " << formatCivilTime(tropStats.minEpoch)
         << ", max epoch: " << formatCivilTime(tropStats.maxEpoch) << endl;
    if (dualStats.n > 0) {
        cout << "Dual-model ionosphere residual: mean=" << dualStats.mean
             << " min=" << dualStats.minVal
             << " max=" << dualStats.maxVal << " m" << endl;
        cout << "  min epoch: " << formatCivilTime(dualStats.minEpoch)
             << ", max epoch: " << formatCivilTime(dualStats.maxEpoch) << endl;
    }
    cout << endl;

    cout << "Interpretation" << endl;
    cout << "  Ionosphere varies strongly within the day and also tracks satellite geometry." << endl;
    cout << "  Troposphere is smoother, with most variation driven by elevation-dependent mapping." << endl;
    if (dualStats.n > 0) {
        cout << "  The dual-frequency leveled ionosphere offers a practical check on the broadcast Klobuchar model." << endl;
    }
    cout << endl;

    cout << "Detail CSV  : " << detailCsvFile << endl;
    cout << "Summary CSV : " << summaryCsvFile << endl;
    cout << "Overview SVG: " << svgFile << endl;
#ifdef _WIN32
    cout << "Overview PNG: " << overviewPngFile << endl;
    if (!g01Rows.empty()) {
        cout << "G01 PNG     : " << g01PngFile << endl;
    }
#endif

    obsStream.close();
    navStream.close();
    return 0;
}
