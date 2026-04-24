/**
 * Exercise 5.6 ionosphere.
 *
 * Reproduce textbook example 5-2 with the GPS Klobuchar model, then compare
 * the L1 ionospheric delay estimated from dual-frequency GPS pseudoranges
 * against the broadcast ionosphere model over one day.
 */

#include <algorithm>
#include <array>
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
#include "GnssStruct.h"
#include "NavEphGPS.hpp"
#include "RinexNavStore.hpp"
#include "TimeConvert.h"

#include "navreadheader.h"
#include "navreader.h"
#include "obsreadheader.h"
#include "obsreader.h"

using namespace std;

namespace {

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
    double psiRad{0.0};
    double ippLatRad{0.0};
    double ippLonRad{0.0};
    double ippGeomagLatRad{0.0};
    double localTimeSec{0.0};
    double amplitudeSec{0.0};
    double periodSec{0.0};
    double phaseRad{0.0};
    double mapping{0.0};
    double delaySec{0.0};
    double delayMeters{0.0};
};

struct IonoRow {
    CommonTime epochGps;
    SatID sat;
    double elevationDeg{0.0};
    double azimuthDeg{0.0};
    double c1Meters{0.0};
    double c2Meters{0.0};
    double codeIonoMeters{0.0};
    double phaseIonoMeters{0.0};
    double dualIonoMeters{0.0};
    double modelIonoMeters{0.0};
    double diffMeters{0.0};
};

struct Stats {
    size_t n{0};
    double mean{0.0};
    double stddev{0.0};
    double minVal{0.0};
    double maxVal{0.0};
    double maxAbs{0.0};
    size_t over15{0};
};

struct EpochSummary {
    CommonTime epochGps;
    size_t satCount{0};
    size_t dualCount{0};
    double meanElevationDeg{0.0};
    double meanModelIonoMeters{0.0};
    double meanDualIonoMeters{std::numeric_limits<double>::quiet_NaN()};
    double meanDiffMeters{std::numeric_limits<double>::quiet_NaN()};
};

static double clampUnit(double x) {
    return std::max(-1.0, std::min(1.0, x));
}

static double wrapTwoPi(double x) {
    const double twoPi = 2.0 * PI;
    x = std::fmod(x, twoPi);
    if (x < 0.0) x += twoPi;
    return x;
}

static double wrapPi(double x) {
    x = wrapTwoPi(x);
    if (x > PI) x -= 2.0 * PI;
    return x;
}

static double wrapSecondsOfDay(double sec) {
    sec = std::fmod(sec, static_cast<double>(SEC_PER_DAY));
    if (sec < 0.0) sec += SEC_PER_DAY;
    return sec;
}

static bool finitePositive(double v) {
    return std::isfinite(v) && v > 0.0;
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
        << setw(2) << setfill('0') << civ.day << " "
        << setw(2) << setfill('0') << civ.hour << ":"
        << setw(2) << setfill('0') << civ.minute << ":"
        << fixed << setprecision(3) << setw(6) << setfill('0') << civ.second;
    return oss.str();
}

static double secondsOfDay(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    return civ.hour * 3600.0 + civ.minute * 60.0 + civ.second;
}

static Geodetic ecefToGeodetic(const Eigen::Vector3d &xyz) {
    const double a = RadiusEarth;
    const double f = 1.0 / 298.257223563;
    const double e2 = 2.0 * f - f * f;
    const double p = std::hypot(xyz.x(), xyz.y());
    double lat = std::atan2(xyz.z(), p * (1.0 - e2));
    double h = 0.0;

    for (int i = 0; i < 20; ++i) {
        const double sinLat = std::sin(lat);
        const double n = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
        h = p / std::cos(lat) - n;
        const double nextLat = std::atan2(xyz.z(), p * (1.0 - e2 * n / (n + h)));
        if (std::abs(nextLat - lat) < 1.0e-14) {
            lat = nextLat;
            break;
        }
        lat = nextLat;
    }

    Geodetic out;
    out.latRad = lat;
    out.lonRad = std::atan2(xyz.y(), xyz.x());
    out.hMeters = h;
    return out;
}

static AzEl computeAzEl(const Eigen::Vector3d &recXYZ,
                        const Geodetic &recBLH,
                        const Eigen::Vector3d &satXYZ) {
    const Eigen::Vector3d los = (satXYZ - recXYZ).normalized();
    const double lat = recBLH.latRad;
    const double lon = recBLH.lonRad;

    const Eigen::Vector3d east(-std::sin(lon), std::cos(lon), 0.0);
    const Eigen::Vector3d north(-std::sin(lat) * std::cos(lon),
                                -std::sin(lat) * std::sin(lon),
                                std::cos(lat));
    const Eigen::Vector3d up(std::cos(lat) * std::cos(lon),
                             std::cos(lat) * std::sin(lon),
                             std::sin(lat));

    AzEl out;
    out.azRad = wrapTwoPi(std::atan2(los.dot(east), los.dot(north)));
    out.elRad = std::asin(clampUnit(los.dot(up)));
    return out;
}

static Eigen::Vector3d rotateSatelliteToReceive(const Eigen::Vector3d &recXYZ,
                                                const Xvt &xvtTx) {
    const double travelSec = (xvtTx.x - recXYZ).norm() / C_MPS;
    const double wt = OMEGA_EARTH * travelSec;
    Eigen::Vector3d rot;
    rot.x() = std::cos(wt) * xvtTx.x.x() + std::sin(wt) * xvtTx.x.y();
    rot.y() = -std::sin(wt) * xvtTx.x.x() + std::cos(wt) * xvtTx.x.y();
    rot.z() = xvtTx.x.z();
    return rot;
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

static KlobucharResult computeKlobuchar(const Geodetic &recBLH,
                                        double azRad,
                                        double elevRad,
                                        double gpsSecondsOfDay,
                                        const KlobucharCoeffs &coeffs,
                                        double signalFreqHz = L1_FREQ_GPS) {
    const double geomagPoleLat = 78.3 * DEG_TO_RAD;
    const double geomagPoleLon = 283.0 * DEG_TO_RAD;

    KlobucharResult out;
    const double elevSemiCircle = elevRad / PI;
    out.psiRad = (0.0137 / (elevSemiCircle + 0.11) - 0.022) * PI;

    out.ippLatRad = std::asin(clampUnit(std::sin(recBLH.latRad) * std::cos(out.psiRad)
                      + std::cos(recBLH.latRad) * std::sin(out.psiRad) * std::cos(azRad)));
    out.ippLonRad = wrapPi(recBLH.lonRad + out.psiRad * std::sin(azRad) / std::cos(out.ippLatRad));
    out.ippGeomagLatRad = std::asin(clampUnit(std::sin(out.ippLatRad) * std::sin(geomagPoleLat)
                              + std::cos(out.ippLatRad) * std::cos(geomagPoleLat)
                              * std::cos(out.ippLonRad - geomagPoleLon)));
    out.localTimeSec = wrapSecondsOfDay(43200.0 * out.ippLonRad / PI + gpsSecondsOfDay);

    const double phi = out.ippGeomagLatRad / PI;
    double powPhi = 1.0;
    for (int i = 0; i < 4; ++i) {
        out.amplitudeSec += coeffs.alpha[i] * powPhi;
        out.periodSec += coeffs.beta[i] * powPhi;
        powPhi *= phi;
    }
    if (out.amplitudeSec < 0.0) out.amplitudeSec = 0.0;
    if (out.periodSec < 72000.0) out.periodSec = 72000.0;

    out.phaseRad = 2.0 * PI * (out.localTimeSec - 50400.0) / out.periodSec;
    out.mapping = 1.0 + 16.0 * std::pow(0.53 - elevSemiCircle, 3.0);

    if (std::abs(out.phaseRad) < PI / 2.0) {
        out.delaySec = (5.0e-9 + out.amplitudeSec * std::cos(out.phaseRad)) * out.mapping;
    } else {
        out.delaySec = 5.0e-9 * out.mapping;
    }

    const double scale = (L1_FREQ_GPS / signalFreqHz) * (L1_FREQ_GPS / signalFreqHz);
    out.delaySec *= scale;
    out.delayMeters = out.delaySec * C_MPS;
    return out;
}

static KlobucharCoeffs findGpsKlobucharCoeffs(const NavHeaderAll &navHeader) {
    KlobucharCoeffs coeffs;
    bool hasAlpha = false;
    bool hasBeta = false;

    for (const auto &rec : navHeader.ionoCorrs) {
        if (rec.type == "GPSA" && rec.coeffs.size() >= 4) {
            std::copy_n(rec.coeffs.begin(), 4, coeffs.alpha.begin());
            hasAlpha = true;
        } else if (rec.type == "GPSB" && rec.coeffs.size() >= 4) {
            std::copy_n(rec.coeffs.begin(), 4, coeffs.beta.begin());
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
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    double med = values[mid];
    if (values.size() % 2 == 0) {
        std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
        med = 0.5 * (med + values[mid - 1]);
    }
    return med;
}

static void appendLeveledArc(const vector<IonoRow> &arc,
                             vector<IonoRow> &leveledRows,
                             size_t &rejectedShortArc,
                             size_t &rejectedOutlier) {
    if (arc.size() < 10) {
        rejectedShortArc += arc.size();
        return;
    }

    vector<double> leveling;
    leveling.reserve(arc.size());
    for (const auto &row : arc) {
        leveling.push_back(row.codeIonoMeters - row.phaseIonoMeters);
    }
    const double bias = median(leveling);

    for (auto row : arc) {
        row.dualIonoMeters = row.phaseIonoMeters + bias;
        row.diffMeters = row.dualIonoMeters - row.modelIonoMeters;
        if (!std::isfinite(row.diffMeters) || std::abs(row.diffMeters) > 15.0) {
            ++rejectedOutlier;
            continue;
        }
        leveledRows.push_back(row);
    }
}

static vector<IonoRow> levelPhaseByCode(const map<SatID, vector<IonoRow>> &rawRowsBySat,
                                        size_t &rejectedShortArc,
                                        size_t &rejectedOutlier) {
    vector<IonoRow> leveledRows;
    rejectedShortArc = 0;
    rejectedOutlier = 0;

    for (const auto &satRows : rawRowsBySat) {
        vector<IonoRow> arc;
        CommonTime prevEpoch;
        double prevPhaseIono = 0.0;
        bool hasPrev = false;

        for (const auto &row : satRows.second) {
            bool newArc = false;
            if (hasPrev) {
                const double dt = row.epochGps - prevEpoch;
                if (dt > 120.0 || std::abs(row.phaseIonoMeters - prevPhaseIono) > 5.0) {
                    newArc = true;
                }
            }

            if (newArc) {
                appendLeveledArc(arc, leveledRows, rejectedShortArc, rejectedOutlier);
                arc.clear();
            }

            arc.push_back(row);
            prevEpoch = row.epochGps;
            prevPhaseIono = row.phaseIonoMeters;
            hasPrev = true;
        }

        appendLeveledArc(arc, leveledRows, rejectedShortArc, rejectedOutlier);
    }

    std::sort(leveledRows.begin(), leveledRows.end(), [](const IonoRow &a, const IonoRow &b) {
        const double dt = a.epochGps - b.epochGps;
        if (std::abs(dt) > 1.0e-9) return dt < 0.0;
        return a.sat < b.sat;
    });
    return leveledRows;
}

static Stats computeStats(const vector<double> &values) {
    Stats st;
    if (values.empty()) return st;

    st.n = values.size();
    st.mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(st.n);
    double ss = 0.0;
    st.minVal = values.front();
    st.maxVal = values.front();
    for (double v : values) {
        const double d = v - st.mean;
        ss += d * d;
        st.minVal = std::min(st.minVal, v);
        st.maxVal = std::max(st.maxVal, v);
        st.maxAbs = std::max(st.maxAbs, std::abs(v));
        if (std::abs(v) > 15.0) ++st.over15;
    }
    st.stddev = (st.n > 1) ? std::sqrt(ss / static_cast<double>(st.n - 1)) : 0.0;
    return st;
}

static vector<EpochSummary> summarizeByEpoch(const vector<IonoRow> &rows) {
    struct Accumulator {
        size_t satCount{0};
        size_t dualCount{0};
        double elevationSum{0.0};
        double modelSum{0.0};
        double dualSum{0.0};
        double diffSum{0.0};
    };

    map<CommonTime, Accumulator> acc;
    for (const auto &row : rows) {
        auto &item = acc[row.epochGps];
        ++item.satCount;
        item.elevationSum += row.elevationDeg;
        item.modelSum += row.modelIonoMeters;
        if (std::isfinite(row.dualIonoMeters)) {
            ++item.dualCount;
            item.dualSum += row.dualIonoMeters;
            item.diffSum += row.diffMeters;
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
        row.meanModelIonoMeters = entry.second.modelSum / static_cast<double>(entry.second.satCount);
        if (entry.second.dualCount > 0) {
            row.meanDualIonoMeters = entry.second.dualSum / static_cast<double>(entry.second.dualCount);
            row.meanDiffMeters = entry.second.diffSum / static_cast<double>(entry.second.dualCount);
        }
        out.push_back(row);
    }
    return out;
}

static void padRange(double &mn, double &mx, double ratio) {
    if (!std::isfinite(mn) || !std::isfinite(mx)) {
        mn = -1.0;
        mx = 1.0;
        return;
    }
    if (mx <= mn) {
        const double delta = std::max(1.0, std::abs(mx) * ratio);
        mn -= delta;
        mx += delta;
        return;
    }
    const double span = mx - mn;
    mn -= span * ratio;
    mx += span * ratio;
}

static void writeDistributionSvg(const string &svgFile,
                                 const vector<IonoRow> &rows,
                                 const Stats &stats) {
    if (rows.empty()) return;

    ofstream svg(svgFile);
    if (!svg) throw runtime_error("failed to open svg file: " + svgFile);

    vector<double> values;
    values.reserve(rows.size());
    double timeMin = numeric_limits<double>::infinity();
    double timeMax = -numeric_limits<double>::infinity();
    double diffMin = numeric_limits<double>::infinity();
    double diffMax = -numeric_limits<double>::infinity();
    for (const auto &row : rows) {
        const double t = secondsOfDay(row.epochGps);
        values.push_back(row.diffMeters);
        timeMin = std::min(timeMin, t);
        timeMax = std::max(timeMax, t);
        diffMin = std::min(diffMin, row.diffMeters);
        diffMax = std::max(diffMax, row.diffMeters);
    }
    diffMin = std::min(diffMin, -15.0);
    diffMax = std::max(diffMax, 15.0);
    padRange(diffMin, diffMax, 0.08);

    const int width = 1280;
    const int height = 760;
    const int marginL = 76;
    const int marginR = 42;
    const int marginT = 76;
    const int marginB = 70;
    const int gap = 70;
    const int panelW = (width - marginL - marginR - gap) / 2;
    const int plotH = height - marginT - marginB;
    const int histX0 = marginL;
    const int scatterX0 = marginL + panelW + gap;
    const int plotY0 = marginT;
    const int plotY1 = marginT + plotH;

    const int binCount = 40;
    vector<int> bins(binCount, 0);
    const double binWidth = (diffMax - diffMin) / binCount;
    for (double v : values) {
        int idx = static_cast<int>((v - diffMin) / binWidth);
        if (idx < 0) idx = 0;
        if (idx >= binCount) idx = binCount - 1;
        ++bins[idx];
    }
    const int maxBin = *std::max_element(bins.begin(), bins.end());
    double yHistMax = std::max(1.0, static_cast<double>(maxBin) * 1.18);
    if (stats.stddev > 0.0) {
        const double normalPeak = values.size() * binWidth / (stats.stddev * std::sqrt(2.0 * PI));
        yHistMax = std::max(yHistMax, normalPeak * 1.18);
    }

    auto histX = [&](double x) {
        return histX0 + (x - diffMin) / (diffMax - diffMin) * panelW;
    };
    auto histY = [&](double y) {
        return plotY1 - y / yHistMax * plotH;
    };
    auto timeX = [&](double t) {
        if (timeMax <= timeMin) return static_cast<double>(scatterX0);
        return scatterX0 + (t - timeMin) / (timeMax - timeMin) * panelW;
    };
    auto diffY = [&](double y) {
        return plotY1 - (y - diffMin) / (diffMax - diffMin) * plotH;
    };

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    svg << "<style>"
        << ".title{font:700 24px 'Arial','Microsoft YaHei',sans-serif;fill:#111827;}"
        << ".axis{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << ".tick{font:12px 'Arial','Microsoft YaHei',sans-serif;fill:#6b7280;}"
        << ".note{font:13px 'Arial','Microsoft YaHei',sans-serif;fill:#4b5563;}"
        << ".panel{font:700 16px 'Arial','Microsoft YaHei',sans-serif;fill:#1f2937;}"
        << "</style>\n";

    svg << "<text x=\"" << marginL << "\" y=\"38\" class=\"title\">Exercise 5.6: dual-frequency ionosphere minus Klobuchar model</text>\n";
    svg << "<text x=\"" << marginL << "\" y=\"60\" class=\"note\">"
        << "N=" << stats.n
        << ", mean=" << fixed << setprecision(3) << stats.mean << " m"
        << ", std=" << stats.stddev << " m"
        << ", min/max=" << stats.minVal << "/" << stats.maxVal << " m"
        << ", |diff|&gt;15 m: " << stats.over15
        << "</text>\n";

    svg << "<rect x=\"" << histX0 << "\" y=\"" << plotY0 << "\" width=\"" << panelW
        << "\" height=\"" << plotH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";
    svg << "<rect x=\"" << scatterX0 << "\" y=\"" << plotY0 << "\" width=\"" << panelW
        << "\" height=\"" << plotH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";
    svg << "<text x=\"" << (histX0 + 12) << "\" y=\"" << (plotY0 + 24)
        << "\" class=\"panel\">Histogram with fitted normal curve</text>\n";
    svg << "<text x=\"" << (scatterX0 + 12) << "\" y=\"" << (plotY0 + 24)
        << "\" class=\"panel\">Residuals by epoch</text>\n";

    for (int i = 0; i <= 6; ++i) {
        const double frac = static_cast<double>(i) / 6.0;
        const double y = plotY1 - frac * plotH;
        const double histVal = frac * yHistMax;
        const double diffVal = diffMin + frac * (diffMax - diffMin);
        svg << "<line x1=\"" << histX0 << "\" y1=\"" << y << "\" x2=\""
            << (histX0 + panelW) << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\"/>\n";
        svg << "<line x1=\"" << scatterX0 << "\" y1=\"" << y << "\" x2=\""
            << (scatterX0 + panelW) << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\"/>\n";
        svg << "<text x=\"" << (histX0 - 10) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\" class=\"tick\">" << fixed << setprecision(0) << histVal << "</text>\n";
        svg << "<text x=\"" << (scatterX0 - 10) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\" class=\"tick\">" << fixed << setprecision(1) << diffVal << "</text>\n";
    }

    for (int i = 0; i <= 8; ++i) {
        const double frac = static_cast<double>(i) / 8.0;
        const double xh = histX0 + frac * panelW;
        const double xs = scatterX0 + frac * panelW;
        const double diffVal = diffMin + frac * (diffMax - diffMin);
        const double sec = timeMin + frac * (timeMax - timeMin);
        int hour = static_cast<int>(sec / 3600.0);
        int minute = static_cast<int>(std::fmod(sec, 3600.0) / 60.0);
        ostringstream hhmm;
        hhmm << setfill('0') << setw(2) << hour << ":" << setw(2) << minute;

        svg << "<line x1=\"" << xh << "\" y1=\"" << plotY0 << "\" x2=\"" << xh
            << "\" y2=\"" << plotY1 << "\" stroke=\"#eef2f7\"/>\n";
        svg << "<line x1=\"" << xs << "\" y1=\"" << plotY0 << "\" x2=\"" << xs
            << "\" y2=\"" << plotY1 << "\" stroke=\"#eef2f7\"/>\n";
        svg << "<text x=\"" << xh << "\" y=\"" << (plotY1 + 22)
            << "\" text-anchor=\"middle\" class=\"tick\">" << fixed << setprecision(1) << diffVal << "</text>\n";
        svg << "<text x=\"" << xs << "\" y=\"" << (plotY1 + 22)
            << "\" text-anchor=\"middle\" class=\"tick\">" << hhmm.str() << "</text>\n";
    }

    for (int i = 0; i < binCount; ++i) {
        const double x0 = histX(diffMin + i * binWidth);
        const double x1 = histX(diffMin + (i + 1) * binWidth);
        const double y = histY(bins[i]);
        svg << "<rect x=\"" << (x0 + 1.0) << "\" y=\"" << y
            << "\" width=\"" << std::max(0.0, x1 - x0 - 2.0)
            << "\" height=\"" << (plotY1 - y)
            << "\" fill=\"#70a6ff\" opacity=\"0.78\"/>\n";
    }

    if (stats.stddev > 0.0) {
        svg << "<polyline fill=\"none\" stroke=\"#1f2937\" stroke-width=\"2.5\" points=\"";
        const int curveN = 240;
        for (int i = 0; i <= curveN; ++i) {
            const double x = diffMin + (diffMax - diffMin) * i / curveN;
            const double z = (x - stats.mean) / stats.stddev;
            const double y = values.size() * binWidth
                             * std::exp(-0.5 * z * z)
                             / (stats.stddev * std::sqrt(2.0 * PI));
            svg << histX(x) << "," << histY(y) << " ";
        }
        svg << "\"/>\n";
    }

    for (double ref : {-15.0, 15.0}) {
        if (ref >= diffMin && ref <= diffMax) {
            svg << "<line x1=\"" << histX(ref) << "\" y1=\"" << plotY0
                << "\" x2=\"" << histX(ref) << "\" y2=\"" << plotY1
                << "\" stroke=\"#ef4444\" stroke-dasharray=\"8 6\" stroke-width=\"1.6\"/>\n";
            svg << "<line x1=\"" << scatterX0 << "\" y1=\"" << diffY(ref)
                << "\" x2=\"" << (scatterX0 + panelW) << "\" y2=\"" << diffY(ref)
                << "\" stroke=\"#ef4444\" stroke-dasharray=\"8 6\" stroke-width=\"1.6\"/>\n";
        }
    }

    const size_t pointStep = std::max<size_t>(1, rows.size() / 5000);
    for (size_t i = 0; i < rows.size(); i += pointStep) {
        const auto &row = rows[i];
        svg << "<circle cx=\"" << timeX(secondsOfDay(row.epochGps))
            << "\" cy=\"" << diffY(row.diffMeters)
            << "\" r=\"1.7\" fill=\"#2563eb\" opacity=\"0.55\"/>\n";
    }

    svg << "<line x1=\"" << histX0 << "\" y1=\"" << plotY1 << "\" x2=\""
        << (histX0 + panelW) << "\" y2=\"" << plotY1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << histX0 << "\" y1=\"" << plotY0 << "\" x2=\""
        << histX0 << "\" y2=\"" << plotY1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << scatterX0 << "\" y1=\"" << plotY1 << "\" x2=\""
        << (scatterX0 + panelW) << "\" y2=\"" << plotY1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << scatterX0 << "\" y1=\"" << plotY0 << "\" x2=\""
        << scatterX0 << "\" y2=\"" << plotY1 << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";

    svg << "<text x=\"" << (histX0 + panelW / 2) << "\" y=\"" << (height - 18)
        << "\" text-anchor=\"middle\" class=\"axis\">dual-frequency minus model / m</text>\n";
    svg << "<text x=\"" << (scatterX0 + panelW / 2) << "\" y=\"" << (height - 18)
        << "\" text-anchor=\"middle\" class=\"axis\">epoch time of day</text>\n";
    svg << "<text x=\"24\" y=\"" << (plotY0 + plotH / 2)
        << "\" transform=\"rotate(-90 24 " << (plotY0 + plotH / 2)
        << ")\" text-anchor=\"middle\" class=\"axis\">count</text>\n";
    svg << "<text x=\"" << (scatterX0 - 52) << "\" y=\"" << (plotY0 + plotH / 2)
        << "\" transform=\"rotate(-90 " << (scatterX0 - 52) << " " << (plotY0 + plotH / 2)
        << ")\" text-anchor=\"middle\" class=\"axis\">difference / m</text>\n";

    svg << "</svg>\n";
}

static void writeCsv(const string &csvFile, const vector<IonoRow> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open csv file: " + csvFile);

    csv << "epoch,sat,elevation_deg,azimuth_deg,c1w_m,c2w_m,code_l1_iono_m,phase_l1_iono_m,"
        << "leveled_dual_l1_iono_m,klobuchar_l1_iono_m,diff_m\n";
    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.sat << ","
            << fixed << setprecision(6)
            << row.elevationDeg << ","
            << row.azimuthDeg << ","
            << setprecision(3)
            << row.c1Meters << ","
            << row.c2Meters << ","
            << row.codeIonoMeters << ","
            << row.phaseIonoMeters << ","
            << row.dualIonoMeters << ","
            << row.modelIonoMeters << ","
            << row.diffMeters << "\n";
    }
}

static void writeSummaryCsv(const string &csvFile, const vector<EpochSummary> &rows) {
    ofstream csv(csvFile);
    if (!csv) throw runtime_error("failed to open summary csv file: " + csvFile);

    csv << "epoch,sat_count,dual_count,mean_elevation_deg,mean_model_l1_iono_m,"
        << "mean_dual_l1_iono_m,mean_diff_m\n";
    for (const auto &row : rows) {
        csv << formatCivilTime(row.epochGps) << ","
            << row.satCount << ","
            << row.dualCount << ","
            << fixed << setprecision(6)
            << row.meanElevationDeg << ","
            << row.meanModelIonoMeters << ",";
        if (std::isfinite(row.meanDualIonoMeters)) {
            csv << row.meanDualIonoMeters;
        }
        csv << ",";
        if (std::isfinite(row.meanDiffMeters)) {
            csv << row.meanDiffMeters;
        }
        csv << "\n";
    }
}

static void printExample52(const KlobucharResult &res) {
    struct RefItem {
        string name;
        double calc;
        double ref;
        string unit;
        int precision;
    };

    const vector<RefItem> items = {
        {"psi", res.psiRad * RAD_TO_DEG, 0.357594, "deg", 6},
        {"IPP latitude", res.ippLatRad * RAD_TO_DEG, 30.504584, "deg", 6},
        {"IPP longitude", res.ippLonRad * RAD_TO_DEG, 113.943409, "deg", 6},
        {"IPP geomagnetic latitude", res.ippGeomagLatRad * RAD_TO_DEG, 18.999169, "deg", 6},
        {"IPP local time", res.localTimeSec, 27346.41818396, "s", 6},
        {"Amplitude", res.amplitudeSec, 0.0000000307, "s", 12},
        {"Period", res.periodSec, 133055.103903080, "s", 6},
        {"Phase", res.phaseRad, -1.0886461495, "rad", 10},
        {"Mapping", res.mapping, 1.005216221, "", 9},
        {"Ionospheric delay", res.delaySec, 0.0000000194, "s", 12},
    };

    cout << "[Example 5-2 reproduction]" << endl;
    cout << left << setw(30) << "Parameter"
         << right << setw(22) << "Computed"
         << setw(22) << "Reference"
         << setw(12) << "Diff"
         << "  Unit" << endl;
    cout << string(96, '-') << endl;
    for (const auto &item : items) {
        cout << left << setw(30) << item.name
             << right << setw(22) << fixed << setprecision(item.precision) << item.calc
             << setw(22) << fixed << setprecision(item.precision) << item.ref
             << setw(12) << scientific << setprecision(3) << (item.calc - item.ref)
             << "  " << item.unit << endl;
    }
    cout << fixed << setprecision(3)
         << "Ionospheric delay on L1: " << res.delayMeters << " m" << endl
         << endl;
}

} // namespace

int main() {
    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\";
    const string obsFile = dataDir + "WUH200CHN_R_20250010000_01D_30S_MO.rnx";
    const string navFile = dataDir + "BRDC00IGS_R_20250010000_01D_MN.rnx";
    const string detailCsvFile = dataDir + "exam-5.6-ionosphere.csv";
    const string summaryCsvFile = dataDir + "exam-5.6-ionosphere_summary.csv";
    const string svgFile = dataDir + "exam-5.6-ionosphere.svg";

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
    const KlobucharCoeffs coeffs = findGpsKlobucharCoeffs(navHeader);

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

    cout << "Exercise 5.6 ionosphere: Klobuchar model and dual-frequency check" << endl;
    cout << "Observation file: " << obsFile << endl;
    cout << "Navigation file : " << navFile << endl;
    cout << "Station         : " << obsHeader.markerName << endl;
    cout << "Receiver XYZ    : (" << fixed << setprecision(3)
         << recXYZ.x() << ", " << recXYZ.y() << ", " << recXYZ.z() << ") m" << endl;
    cout << "Receiver BLH    : lat=" << setprecision(8) << recBLH.latRad * RAD_TO_DEG
         << " deg, lon=" << recBLH.lonRad * RAD_TO_DEG
         << " deg, h=" << setprecision(3) << recBLH.hMeters << " m" << endl;
    cout << "GPS nav records : " << gpsRecordCount << endl;
    cout << "GPSA            : " << scientific << setprecision(6)
         << coeffs.alpha[0] << ", " << coeffs.alpha[1] << ", "
         << coeffs.alpha[2] << ", " << coeffs.alpha[3] << endl;
    cout << "GPSB            : "
         << coeffs.beta[0] << ", " << coeffs.beta[1] << ", "
         << coeffs.beta[2] << ", " << coeffs.beta[3] << endl << endl;

    const double bookAzDeg = 265.65936;
    const double bookElDeg = 83.007251;
    const double bookTowSec = 0.0;
    const KlobucharResult book = computeKlobuchar(recBLH,
                                                  bookAzDeg * DEG_TO_RAD,
                                                  bookElDeg * DEG_TO_RAD,
                                                  bookTowSec,
                                                  coeffs);
    printExample52(book);

    map<SatID, vector<IonoRow>> rawRowsBySat;
    size_t seenGpsWithCodes = 0;
    size_t rejectedLowElev = 0;
    size_t rejectedNav = 0;
    size_t rejectedBadIono = 0;

    while (true) {
        ObsEpochData epochData;
        try {
            epochData = obsReader.readNextEpoch();
        } catch (const EndOfObsFile &) {
            break;
        }

        const CommonTime epochGps = toGpsTime(epochData.time);
        const double sod = secondsOfDay(epochGps);

        for (const auto &satEntry : epochData.satTypeValue) {
            const SatID sat(satEntry.first);
            if (sat.system != "G") continue;

            const auto &tv = satEntry.second;
            auto itC1W = tv.find("C1W");
            auto itC2W = tv.find("C2W");
            auto itC1C = tv.find("C1C");
            auto itL1W = tv.find("L1W");
            auto itL2W = tv.find("L2W");
            if (itC1W == tv.end() || itC2W == tv.end()
                || itL1W == tv.end() || itL2W == tv.end()
                || !finitePositive(itC1W->second) || !finitePositive(itC2W->second)
                || !std::isfinite(itL1W->second) || !std::isfinite(itL2W->second)) {
                continue;
            }
            ++seenGpsWithCodes;

            const double rangeForTx = (itC1C != tv.end() && finitePositive(itC1C->second))
                                      ? itC1C->second
                                      : itC1W->second;

            Xvt xvtTx;
            try {
                xvtTx = computeTransmitXvt(navStore, sat, epochGps, rangeForTx);
            } catch (...) {
                ++rejectedNav;
                continue;
            }

            const Eigen::Vector3d satXYZRec = rotateSatelliteToReceive(recXYZ, xvtTx);
            const AzEl azel = computeAzEl(recXYZ, recBLH, satXYZRec);
            const double elevDeg = azel.elRad * RAD_TO_DEG;
            if (elevDeg < 15.0) {
                ++rejectedLowElev;
                continue;
            }

            const KlobucharResult model = computeKlobuchar(recBLH, azel.azRad, azel.elRad, sod, coeffs);
            const double codeIono = dualFreqL1IonoMeters(itC1W->second, itC2W->second);
            const double phaseIono = dualFreqL1PhaseIonoMeters(itL1W->second, itL2W->second);
            if (!std::isfinite(codeIono) || !std::isfinite(phaseIono)
                || std::abs(codeIono) > 100.0 || std::abs(phaseIono) > 1.0e8) {
                ++rejectedBadIono;
                continue;
            }

            IonoRow row;
            row.epochGps = epochGps;
            row.sat = sat;
            row.elevationDeg = elevDeg;
            row.azimuthDeg = azel.azRad * RAD_TO_DEG;
            row.c1Meters = itC1W->second;
            row.c2Meters = itC2W->second;
            row.codeIonoMeters = codeIono;
            row.phaseIonoMeters = phaseIono;
            row.modelIonoMeters = model.delayMeters;
            rawRowsBySat[sat].push_back(row);
        }
    }

    size_t rejectedShortArc = 0;
    size_t rejectedOutlier = 0;
    vector<IonoRow> rows = levelPhaseByCode(rawRowsBySat, rejectedShortArc, rejectedOutlier);
    const vector<EpochSummary> summaryRows = summarizeByEpoch(rows);

    vector<double> diffs;
    diffs.reserve(rows.size());
    for (const auto &row : rows) diffs.push_back(row.diffMeters);
    const Stats stats = computeStats(diffs);

    writeCsv(detailCsvFile, rows);
    writeSummaryCsv(summaryCsvFile, summaryRows);
    writeDistributionSvg(svgFile, rows, stats);

    cout << "--- Dual-frequency minus Klobuchar summary ---" << endl;
    cout << "GPS samples with C1W/C2W : " << seenGpsWithCodes << endl;
    cout << "Accepted samples         : " << stats.n << endl;
    cout << "Rejected low elevation   : " << rejectedLowElev << endl;
    cout << "Rejected missing nav      : " << rejectedNav << endl;
    cout << "Rejected bad ionosphere   : " << rejectedBadIono << endl;
    cout << "Rejected short arcs       : " << rejectedShortArc << endl;
    cout << "Rejected |diff| > 15 m    : " << rejectedOutlier << endl;
    cout << fixed << setprecision(3)
         << "Difference mean/std      : " << stats.mean << " / " << stats.stddev << " m" << endl
         << "Difference min/max       : " << stats.minVal << " / " << stats.maxVal << " m" << endl
         << "Max |difference|         : " << stats.maxAbs << " m" << endl
         << "|difference| > 15 m      : " << stats.over15 << endl;
    cout << "Detail CSV written to    : " << detailCsvFile << endl;
    cout << "Summary CSV written to   : " << summaryCsvFile << endl;
    cout << "SVG written to           : " << svgFile << endl;

    obsStream.close();
    navStream.close();
    return 0;
}
