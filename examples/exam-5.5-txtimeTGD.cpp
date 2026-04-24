/**
 * Exercise 5.5 example.
 *
 * Reproduce the GPS C1 TGD correction analysis:
 * - read RINEX headers with NavReadHeader / ObsReadHeader
 * - read observation epochs with ObsReader
 * - read GPS broadcast ephemeris with NavReader and store it in RinexNavStore
 * - apply broadcast TGD to GPS C1 using correctTGD
 * - emit a CSV file and two SVG figures:
 *   1) back-substitution residuals
 *   2) G10 C1 TGD correction time series
 */

#include <cmath>
#include <cerrno>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Const.h"
#include "GnssFunc.h"
#include "NavEphGPS.hpp"
#include "RinexNavStore.hpp"
#include "TimeConvert.h"

#include "navreadheader.h"
#include "navreader.h"
#include "obsreadheader.h"
#include "obsreader.h"

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

struct TxRowSummary {
    CommonTime epochGps;
    SatID sat;
    double c1Raw{0.0};
    double c1Tgd{0.0};
    double c2{0.0};
    double ifPr{0.0};
    double tgdSeconds{0.0};
    double tgdMeters{0.0};
    CommonTime txRaw;
    CommonTime txTgd;
    CommonTime txIf;
    double txRawOffsetSec{0.0};
    double txTgdOffsetSec{0.0};
    double txIfOffsetSec{0.0};
    double dtTgdNs{0.0};
    double dtIfNs{0.0};
    double dtRawIfNs{0.0};
    double dc1Meters{0.0};
    double c1RawBack{0.0};
    double c1TgdBack{0.0};
    double ifBack{0.0};
    double c1RawResidual{0.0};
    double c1TgdResidual{0.0};
    double ifResidual{0.0};
};

struct SvgSeries {
    std::string name;
    std::string color;
    std::vector<std::pair<double, double>> points;
};

struct TxSolution {
    CommonTime tx;
    double recvMinusTxSec{0.0};
};

static CommonTime toGpsTime(const RnxTime &rt) {
    CivilTime civ(rt.year, rt.month, rt.day, rt.hour, rt.minute, rt.second, TimeSystem::GPS);
    CommonTime ct = CivilTime2CommonTime(civ);
    ct.setTimeSystem(TimeSystem::GPS);
    return ct;
}

static CommonTime toGpsTime(const CommonTime &ct) {
    CommonTime gpsTime = convertTimeSystem(ct, TimeSystem::GPS);
    gpsTime.setTimeSystem(TimeSystem::GPS);
    return gpsTime;
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

static double rawGpsClockBiasSec(const NavEphGPS &eph, const CommonTime &t) {
    const double dt = t - eph.ctToc;
    return eph.af0 + dt * (eph.af1 + dt * eph.af2);
}

static TxSolution computeTransmitSolution(const NavEphGPS &eph,
                                          const CommonTime &recvGpsTime,
                                          double pseudorangeMeters,
                                          bool applyTgd) {
    const double baseTravelSec = pseudorangeMeters / C_MPS;
    double recvMinusTxSec = baseTravelSec;
    CommonTime tt = recvGpsTime;
    tt -= recvMinusTxSec;

    for (int i = 0; i < 2; ++i) {
        double clkBias = rawGpsClockBiasSec(eph, tt);
        if (applyTgd) {
            clkBias -= eph.TGD;
        }
        double relCorr = eph.svRelativity(tt);
        recvMinusTxSec = baseTravelSec + clkBias + relCorr;
        tt = recvGpsTime;
        tt -= recvMinusTxSec;
    }

    return {tt, recvMinusTxSec};
}

static CommonTime computeTransmitTime(const NavEphGPS &eph,
                                      const CommonTime &recvGpsTime,
                                      double pseudorangeMeters,
                                      bool applyTgd) {
    return computeTransmitSolution(eph, recvGpsTime, pseudorangeMeters, applyTgd).tx;
}

static double computeIfPseudorange(const TypeValueMap &tv, const string &sys) {
    auto itC1 = tv.find("C1C");
    auto itC2 = tv.find("C2W");
    if (itC1 == tv.end() || itC2 == tv.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double f1 = getFreq(sys, "C1");
    const double f2 = getFreq(sys, "C2");
    if (f1 == 0.0 || f2 == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return (f1 * f1 * itC1->second - f2 * f2 * itC2->second) / (f1 * f1 - f2 * f2);
}

static double computeClockRelSeconds(const NavEphGPS &eph,
                                     const CommonTime &tx,
                                     bool applyTgd) {
    double clkBias = rawGpsClockBiasSec(eph, tx);
    if (applyTgd) {
        clkBias -= eph.TGD;
    }
    return clkBias + eph.svRelativity(tx);
}

static double backSubstituteRangeMeters(const NavEphGPS &eph,
                                        const CommonTime &tx,
                                        double recvMinusTxSec,
                                        bool applyTgd) {
    return C_MPS * (recvMinusTxSec - computeClockRelSeconds(eph, tx, applyTgd));
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

static TxRowSummary buildRow(const ObsEpochData &epochData,
                             const SatID &sat,
                             const TypeValueMap &tv,
                             const NavEphGPS &eph) {
    TxRowSummary row;
    row.epochGps = toGpsTime(epochData.time);
    row.sat = sat;
    row.c1Raw = tv.at("C1C");
    row.c2 = tv.at("C2W");
    row.ifPr = computeIfPseudorange(tv, sat.system);

    const TxSolution rawTx = computeTransmitSolution(eph, row.epochGps, row.c1Raw, false);
    row.txRaw = rawTx.tx;
    row.txRawOffsetSec = -rawTx.recvMinusTxSec;

    ObsData delayObs;
    delayObs.epoch = row.epochGps;
    delayObs.antennaPosition = Eigen::Vector3d::Zero();
    delayObs.satTypeValueData[sat]["C1C"] = row.c1Raw;

    std::map<SatID, Xvt> satXvtTransTime;
    satXvtTransTime[sat] = eph.svXvt(row.txRaw);

    std::vector<ObsDelayCorrectionRecord> delayReport;
    correctTGD(delayObs, satXvtTransTime, nullptr, &delayReport);

    row.c1Tgd = delayObs.satTypeValueData.at(sat).at("C1C");
    row.dc1Meters = row.c1Raw - row.c1Tgd;
    if (!delayReport.empty()) {
        row.tgdSeconds = delayReport.front().tgdSeconds;
        row.tgdMeters = delayReport.front().correctionMeters;
    } else {
        row.tgdSeconds = eph.TGD;
        row.tgdMeters = 0.0;
    }

    const TxSolution tgdTx = computeTransmitSolution(eph, row.epochGps, row.c1Tgd, false);
    const TxSolution ifTx = computeTransmitSolution(eph, row.epochGps, row.ifPr, true);
    row.txTgd = tgdTx.tx;
    row.txIf = ifTx.tx;
    row.txTgdOffsetSec = -tgdTx.recvMinusTxSec;
    row.txIfOffsetSec = -ifTx.recvMinusTxSec;
    row.dtTgdNs = (rawTx.recvMinusTxSec - tgdTx.recvMinusTxSec) * 1.0e9;
    row.dtIfNs = (tgdTx.recvMinusTxSec - ifTx.recvMinusTxSec) * 1.0e9;
    row.dtRawIfNs = (rawTx.recvMinusTxSec - ifTx.recvMinusTxSec) * 1.0e9;

    row.c1RawBack = backSubstituteRangeMeters(eph, row.txRaw, rawTx.recvMinusTxSec, false);
    row.c1TgdBack = backSubstituteRangeMeters(eph, row.txTgd, tgdTx.recvMinusTxSec, false);
    row.ifBack = backSubstituteRangeMeters(eph, row.txIf, ifTx.recvMinusTxSec, true);
    row.c1RawResidual = row.c1RawBack - row.c1Raw;
    row.c1TgdResidual = row.c1TgdBack - row.c1Tgd;
    row.ifResidual = row.ifBack - row.ifPr;

    return row;
}

static void printCsvHeader(ofstream &csv) {
    csv << "epoch,sat,c1_raw_m,c1_tgd_m,c2_m,if_m,tgd_s,tgd_m,"
        << "tx_raw_offset_s,tx_tgd_offset_s,tx_if_offset_s,dt_tgd_ns,dt_if_ns,dt_raw_if_ns,dc1_m,"
        << "c1_raw_back_m,c1_tgd_back_m,if_back_m,c1_raw_residual_m,c1_tgd_residual_m,if_residual_m\n";
}

static void writeCsvRow(ofstream &csv, const TxRowSummary &row) {
    csv << formatCivilTime(row.epochGps) << ","
        << row.sat << ","
        << fixed << setprecision(3)
        << row.c1Raw << ","
        << row.c1Tgd << ","
        << row.c2 << ","
        << row.ifPr << ","
        << scientific << setprecision(12)
        << row.tgdSeconds << ","
        << row.tgdMeters << ","
        << row.txRawOffsetSec << ","
        << row.txTgdOffsetSec << ","
        << row.txIfOffsetSec << ","
        << row.dtTgdNs << ","
        << row.dtIfNs << ","
        << row.dtRawIfNs << ","
        << row.dc1Meters << ","
        << row.c1RawBack << ","
        << row.c1TgdBack << ","
        << row.ifBack << ","
        << row.c1RawResidual << ","
        << row.c1TgdResidual << ","
        << row.ifResidual << "\n";
}

static double secondsOfDay(const CommonTime &ct) {
    CivilTime civ = CommonTime2CivilTime(ct);
    return civ.hour * 3600.0 + civ.minute * 60.0 + civ.second;
}

static string svgEscape(const string &text) {
    string out;
    out.reserve(text.size());
    for (char c : text) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out.push_back(c); break;
        }
    }
    return out;
}

static void writeLineChartSvg(const string &svgFile,
                              const string &title,
                              const vector<SvgSeries> &topSeries,
                              const vector<SvgSeries> &bottomSeries,
                              double xMin,
                              double xMax,
                              double topYMin,
                              double topYMax,
                              double bottomYMin,
                              double bottomYMax) {
    ofstream svg(svgFile);
    if (!svg) {
        throw runtime_error("failed to open svg file: " + svgFile);
    }

    const int width = 1280;
    const int height = 900;
    const int marginL = 90;
    const int marginR = 40;
    const int marginT = 70;
    const int marginB = 70;
    const int gap = 40;
    const int panelH = (height - marginT - marginB - gap) / 2;
    const int panelW = width - marginL - marginR;

    auto xToPx = [&](double x) {
        if (xMax <= xMin) return static_cast<double>(marginL);
        return marginL + (x - xMin) / (xMax - xMin) * panelW;
    };

    auto yToPx = [&](double y, double top, double bottom, double yMin, double yMax) {
        if (yMax <= yMin) return bottom;
        return bottom - (y - yMin) / (yMax - yMin) * (bottom - top);
    };

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width
        << " " << height << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    svg << "<style>"
        << ".title{font:700 24px 'Arial','Microsoft YaHei',sans-serif;fill:#1f2937;}"
        << ".axis{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << ".tick{font:12px 'Arial','Microsoft YaHei',sans-serif;fill:#6b7280;}"
        << ".legend{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#111827;}"
        << "</style>\n";
    svg << "<text x=\"" << marginL << "\" y=\"" << 36
        << "\" class=\"title\">" << svgEscape(title) << "</text>\n";

    const double panelTop1 = marginT;
    const double panelBottom1 = marginT + panelH;
    const double panelTop2 = marginT + panelH + gap;
    const double panelBottom2 = panelTop2 + panelH;
    const vector<pair<double, double>> panels = {
        {panelTop1, panelBottom1},
        {panelTop2, panelBottom2}
    };

    const vector<string> panelTitles = {
        "Transmit-time shift (ns)",
        "C1 TGD correction (m)"
    };

    for (size_t pi = 0; pi < panels.size(); ++pi) {
        const double top = panels[pi].first;
        const double bottom = panels[pi].second;
        const double plotH = bottom - top;

        svg << "<rect x=\"" << marginL << "\" y=\"" << top << "\" width=\"" << panelW
            << "\" height=\"" << plotH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";
        svg << "<text x=\"" << (marginL + 12) << "\" y=\"" << (top + 24)
            << "\" class=\"axis\">" << svgEscape(panelTitles[pi]) << "</text>\n";

        const int vGrid = 6;
        const int hGrid = 8;
        const double yMin = (pi == 0 ? topYMin : bottomYMin);
        const double yMax = (pi == 0 ? topYMax : bottomYMax);
        const vector<SvgSeries> &curSeries = (pi == 0 ? topSeries : bottomSeries);
        for (int i = 0; i <= vGrid; ++i) {
            double frac = static_cast<double>(i) / vGrid;
            double y = bottom - frac * plotH;
            svg << "<line x1=\"" << marginL << "\" y1=\"" << y << "\" x2=\"" << (marginL + panelW)
                << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n";
            double val = yMin + frac * (yMax - yMin);
            ostringstream oss;
            oss << fixed << setprecision(1) << val;
            svg << "<text x=\"" << (marginL - 12) << "\" y=\"" << (y + 4)
                << "\" text-anchor=\"end\" class=\"tick\">" << oss.str() << "</text>\n";
        }

        for (int i = 0; i <= hGrid; ++i) {
            double frac = static_cast<double>(i) / hGrid;
            double x = marginL + frac * panelW;
            svg << "<line x1=\"" << x << "\" y1=\"" << top << "\" x2=\"" << x
                << "\" y2=\"" << bottom << "\" stroke=\"#eef2f7\" stroke-width=\"1\"/>\n";
            double sec = xMin + frac * (xMax - xMin);
            int hour = static_cast<int>(sec / 3600.0);
            int minute = static_cast<int>(std::fmod(sec, 3600.0) / 60.0);
            ostringstream oss;
            oss << setfill('0') << setw(2) << hour << ":" << setw(2) << minute;
            svg << "<text x=\"" << x << "\" y=\"" << (bottom + 20)
                << "\" text-anchor=\"middle\" class=\"tick\">" << oss.str() << "</text>\n";
        }

        svg << "<line x1=\"" << marginL << "\" y1=\"" << bottom << "\" x2=\""
            << (marginL + panelW) << "\" y2=\"" << bottom << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
        svg << "<line x1=\"" << marginL << "\" y1=\"" << top << "\" x2=\""
            << marginL << "\" y2=\"" << bottom << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
        svg << "<text x=\"" << (marginL + panelW / 2) << "\" y=\"" << (bottom + 46)
            << "\" text-anchor=\"middle\" class=\"axis\">Epoch time of day</text>\n";
        svg << "<text x=\"" << 26 << "\" y=\"" << (top + plotH / 2)
            << "\" transform=\"rotate(-90 26 " << (top + plotH / 2) << ")\" text-anchor=\"middle\" class=\"axis\">"
            << svgEscape(pi == 0 ? "Shift / ns" : "Correction / m") << "</text>\n";

        for (const auto &s : curSeries) {
            if (s.points.empty()) continue;
            svg << "<polyline fill=\"none\" stroke=\"" << s.color
                << "\" stroke-width=\"2\" points=\"";
            for (const auto &pt : s.points) {
                const double x = xToPx(pt.first);
                const double y = yToPx(pt.second, top + 4.0, bottom - 18.0, yMin, yMax);
                svg << x << "," << y << " ";
            }
            svg << "\"/>\n";
        }
    }

    const int legendX = width - marginR - 240;
    const int legendY = 42;
    svg << "<rect x=\"" << legendX << "\" y=\"" << legendY
        << "\" width=\"220\" height=\"110\" rx=\"8\" ry=\"8\" fill=\"#ffffff\" stroke=\"#cbd5e1\"/>\n";
    int ly = legendY + 28;
    for (const auto &s : topSeries) {
        svg << "<line x1=\"" << (legendX + 16) << "\" y1=\"" << ly - 5 << "\" x2=\""
            << (legendX + 48) << "\" y2=\"" << ly - 5
            << "\" stroke=\"" << s.color << "\" stroke-width=\"3\"/>\n";
        svg << "<text x=\"" << (legendX + 58) << "\" y=\"" << ly
            << "\" class=\"legend\">" << svgEscape(s.name) << "</text>\n";
        ly += 26;
    }
    for (const auto &s : bottomSeries) {
        svg << "<line x1=\"" << (legendX + 16) << "\" y1=\"" << ly - 5 << "\" x2=\""
            << (legendX + 48) << "\" y2=\"" << ly - 5
            << "\" stroke=\"" << s.color << "\" stroke-width=\"3\"/>\n";
        svg << "<text x=\"" << (legendX + 58) << "\" y=\"" << ly
            << "\" class=\"legend\">" << svgEscape(s.name) << "</text>\n";
        ly += 26;
    }

    svg << "</svg>\n";
}

static void writeSingleChartSvg(const string &svgFile,
                                const string &title,
                                const SvgSeries &series,
                                double xMin,
                                double xMax,
                                double yMin,
                                double yMax) {
    ofstream svg(svgFile);
    if (!svg) {
        throw runtime_error("failed to open svg file: " + svgFile);
    }

    const int width = 1280;
    const int height = 560;
    const int marginL = 90;
    const int marginR = 40;
    const int marginT = 70;
    const int marginB = 70;
    const int plotW = width - marginL - marginR;
    const int plotH = height - marginT - marginB;

    auto xToPx = [&](double x) {
        if (xMax <= xMin) return static_cast<double>(marginL);
        return marginL + (x - xMin) / (xMax - xMin) * plotW;
    };
    auto yToPx = [&](double y) {
        if (yMax <= yMin) return static_cast<double>(marginT + plotH);
        return marginT + plotH - (y - yMin) / (yMax - yMin) * plotH;
    };

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width
        << " " << height << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    svg << "<style>"
        << ".title{font:700 24px 'Arial','Microsoft YaHei',sans-serif;fill:#1f2937;}"
        << ".axis{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << ".tick{font:12px 'Arial','Microsoft YaHei',sans-serif;fill:#6b7280;}"
        << "</style>\n";
    svg << "<text x=\"" << marginL << "\" y=\"" << 36
        << "\" class=\"title\">" << svgEscape(title) << "</text>\n";

    svg << "<rect x=\"" << marginL << "\" y=\"" << marginT << "\" width=\"" << plotW
        << "\" height=\"" << plotH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";

    const int vGrid = 6;
    const int hGrid = 8;
    for (int i = 0; i <= vGrid; ++i) {
        double frac = static_cast<double>(i) / vGrid;
        double y = marginT + plotH - frac * plotH;
        svg << "<line x1=\"" << marginL << "\" y1=\"" << y << "\" x2=\"" << (marginL + plotW)
            << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n";
        double val = yMin + frac * (yMax - yMin);
        ostringstream oss;
        oss << fixed << setprecision(2) << val;
        svg << "<text x=\"" << (marginL - 12) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\" class=\"tick\">" << oss.str() << "</text>\n";
    }

    for (int i = 0; i <= hGrid; ++i) {
        double frac = static_cast<double>(i) / hGrid;
        double x = marginL + frac * plotW;
        svg << "<line x1=\"" << x << "\" y1=\"" << marginT << "\" x2=\"" << x
            << "\" y2=\"" << (marginT + plotH) << "\" stroke=\"#eef2f7\" stroke-width=\"1\"/>\n";
        double sec = xMin + frac * (xMax - xMin);
        int hour = static_cast<int>(sec / 3600.0);
        int minute = static_cast<int>(std::fmod(sec, 3600.0) / 60.0);
        ostringstream oss;
        oss << setfill('0') << setw(2) << hour << ":" << setw(2) << minute;
        svg << "<text x=\"" << x << "\" y=\"" << (marginT + plotH + 20)
            << "\" text-anchor=\"middle\" class=\"tick\">" << oss.str() << "</text>\n";
    }

    svg << "<line x1=\"" << marginL << "\" y1=\"" << (marginT + plotH) << "\" x2=\""
        << (marginL + plotW) << "\" y2=\"" << (marginT + plotH)
        << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << marginL << "\" y1=\"" << marginT << "\" x2=\""
        << marginL << "\" y2=\"" << (marginT + plotH) << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<text x=\"" << (marginL + plotW / 2) << "\" y=\"" << (height - 20)
        << "\" text-anchor=\"middle\" class=\"axis\">Epoch time of day</text>\n";
    svg << "<text x=\"" << 26 << "\" y=\"" << (marginT + plotH / 2)
        << "\" transform=\"rotate(-90 26 " << (marginT + plotH / 2) << ")\" text-anchor=\"middle\" class=\"axis\">"
        << "TGD correction / m</text>\n";

    if (!series.points.empty()) {
        svg << "<polyline fill=\"none\" stroke=\"" << series.color
            << "\" stroke-width=\"2\" points=\"";
        for (const auto &pt : series.points) {
            svg << xToPx(pt.first) << "," << yToPx(pt.second) << " ";
        }
        svg << "\"/>\n";
    }

    svg << "<rect x=\"" << (width - marginR - 260) << "\" y=\"" << 42
        << "\" width=\"240\" height=\"54\" rx=\"8\" ry=\"8\" fill=\"#ffffff\" stroke=\"#cbd5e1\"/>\n";
    svg << "<line x1=\"" << (width - marginR - 240) << "\" y1=\"" << 69 << "\" x2=\""
        << (width - marginR - 208) << "\" y2=\"" << 69
        << "\" stroke=\"" << series.color << "\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"" << (width - marginR - 198) << "\" y=\"" << 74
        << "\" class=\"axis\">" << svgEscape(series.name) << "</text>\n";

    svg << "</svg>\n";
}

static void writeMultiSeriesChartSvg(const string &svgFile,
                                     const string &title,
                                     const vector<SvgSeries> &seriesList,
                                     double xMin,
                                     double xMax,
                                     double yMin,
                                     double yMax) {
    ofstream svg(svgFile);
    if (!svg) {
        throw runtime_error("failed to open svg file: " + svgFile);
    }

    const int width = 1280;
    const int height = 560;
    const int marginL = 90;
    const int marginR = 40;
    const int marginT = 70;
    const int marginB = 70;
    const int plotW = width - marginL - marginR;
    const int plotH = height - marginT - marginB;
    const double maxAbsY = std::max(std::abs(yMin), std::abs(yMax));
    double yScale = 1.0;
    string yUnit = "m";
    if (maxAbsY > 0.0 && maxAbsY < 1.0e-6) {
        yScale = 1.0e9;
        yUnit = "nm";
    } else if (maxAbsY < 1.0e-2) {
        yScale = 1.0e3;
        yUnit = "mm";
    } else if (maxAbsY < 1.0e-3) {
        yScale = 1.0e6;
        yUnit = "um";
    }
    const double yMinPlot = yMin * yScale;
    const double yMaxPlot = yMax * yScale;

    auto xToPx = [&](double x) {
        if (xMax <= xMin) return static_cast<double>(marginL);
        return marginL + (x - xMin) / (xMax - xMin) * plotW;
    };
    auto yToPx = [&](double y) {
        if (yMaxPlot <= yMinPlot) return static_cast<double>(marginT + plotH);
        const double yPlot = y * yScale;
        return marginT + plotH - (yPlot - yMinPlot) / (yMaxPlot - yMinPlot) * plotH;
    };

    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width
        << " " << height << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    svg << "<style>"
        << ".title{font:700 24px 'Arial','Microsoft YaHei',sans-serif;fill:#1f2937;}"
        << ".axis{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#374151;}"
        << ".tick{font:12px 'Arial','Microsoft YaHei',sans-serif;fill:#6b7280;}"
        << ".legend{font:14px 'Arial','Microsoft YaHei',sans-serif;fill:#111827;}"
        << "</style>\n";
    svg << "<text x=\"" << marginL << "\" y=\"" << 36
        << "\" class=\"title\">" << svgEscape(title) << "</text>\n";

    svg << "<rect x=\"" << marginL << "\" y=\"" << marginT << "\" width=\"" << plotW
        << "\" height=\"" << plotH << "\" fill=\"#fbfdff\" stroke=\"#cbd5e1\"/>\n";

    const int vGrid = 6;
    const int hGrid = 8;
    for (int i = 0; i <= vGrid; ++i) {
        double frac = static_cast<double>(i) / vGrid;
        double y = marginT + plotH - frac * plotH;
        svg << "<line x1=\"" << marginL << "\" y1=\"" << y << "\" x2=\"" << (marginL + plotW)
            << "\" y2=\"" << y << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n";
        double val = yMinPlot + frac * (yMaxPlot - yMinPlot);
        ostringstream oss;
        oss << fixed << setprecision(yUnit == "m" ? 6 : 3) << val;
        svg << "<text x=\"" << (marginL - 12) << "\" y=\"" << (y + 4)
            << "\" text-anchor=\"end\" class=\"tick\">" << oss.str() << "</text>\n";
    }

    for (int i = 0; i <= hGrid; ++i) {
        double frac = static_cast<double>(i) / hGrid;
        double x = marginL + frac * plotW;
        svg << "<line x1=\"" << x << "\" y1=\"" << marginT << "\" x2=\"" << x
            << "\" y2=\"" << (marginT + plotH) << "\" stroke=\"#eef2f7\" stroke-width=\"1\"/>\n";
        double sec = xMin + frac * (xMax - xMin);
        int hour = static_cast<int>(sec / 3600.0);
        int minute = static_cast<int>(std::fmod(sec, 3600.0) / 60.0);
        ostringstream oss;
        oss << setfill('0') << setw(2) << hour << ":" << setw(2) << minute;
        svg << "<text x=\"" << x << "\" y=\"" << (marginT + plotH + 20)
            << "\" text-anchor=\"middle\" class=\"tick\">" << oss.str() << "</text>\n";
    }

    svg << "<line x1=\"" << marginL << "\" y1=\"" << (marginT + plotH) << "\" x2=\""
        << (marginL + plotW) << "\" y2=\"" << (marginT + plotH)
        << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << marginL << "\" y1=\"" << marginT << "\" x2=\""
        << marginL << "\" y2=\"" << (marginT + plotH) << "\" stroke=\"#475569\" stroke-width=\"1.5\"/>\n";
    svg << "<text x=\"" << (marginL + plotW / 2) << "\" y=\"" << (height - 20)
        << "\" text-anchor=\"middle\" class=\"axis\">Epoch time of day</text>\n";
    svg << "<text x=\"" << 26 << "\" y=\"" << (marginT + plotH / 2)
        << "\" transform=\"rotate(-90 26 " << (marginT + plotH / 2) << ")\" text-anchor=\"middle\" class=\"axis\">"
        << "Residual / " << yUnit << "</text>\n";

    const double maxLineGapSec = 120.0;
    for (const auto &series : seriesList) {
        if (series.points.empty()) continue;
        size_t segmentStart = 0;
        for (size_t i = 1; i <= series.points.size(); ++i) {
            const bool atEnd = (i == series.points.size());
            const bool hasGap = !atEnd && (series.points[i].first - series.points[i - 1].first > maxLineGapSec);
            if (!atEnd && !hasGap) {
                continue;
            }

            if (i > segmentStart) {
                svg << "<polyline fill=\"none\" stroke=\"" << series.color
                    << "\" stroke-width=\"2\" points=\"";
                for (size_t j = segmentStart; j < i; ++j) {
                    const auto &pt = series.points[j];
                    svg << xToPx(pt.first) << "," << yToPx(pt.second) << " ";
                }
                svg << "\"/>\n";
            }
            segmentStart = i;
        }
    }

    const int legendX = width - marginR - 260;
    const int legendY = 42;
    const int legendH = 24 * static_cast<int>(seriesList.size()) + 18;
    svg << "<rect x=\"" << legendX << "\" y=\"" << legendY
        << "\" width=\"240\" height=\"" << legendH
        << "\" rx=\"8\" ry=\"8\" fill=\"#ffffff\" stroke=\"#cbd5e1\"/>\n";
    int ly = legendY + 28;
    for (const auto &series : seriesList) {
        svg << "<line x1=\"" << (legendX + 16) << "\" y1=\"" << ly - 5 << "\" x2=\""
            << (legendX + 48) << "\" y2=\"" << ly - 5
            << "\" stroke=\"" << series.color << "\" stroke-width=\"3\"/>\n";
        svg << "<text x=\"" << (legendX + 58) << "\" y=\"" << ly
            << "\" class=\"legend\">" << svgEscape(series.name) << "</text>\n";
        ly += 26;
    }

    svg << "</svg>\n";
}

} // namespace

int main() {
    const int plotGpsPrn = 10;
    const string plotSatLabel = "G10";
    const string dataDir = "D:\\GNSSLAB\\gnssLab-2.4\\data\\";
    const string obsFile = dataDir + "WUH200CHN_R_20250010000_01D_30S_MO.rnx";
    const string navFile = dataDir + "BRDC00IGS_R_20250010000_01D_MN.rnx";
    const string csvFile = dataDir + "exam-5.5-txtimeTGD.csv";
    const string residualSvgFile = dataDir + "exam-5.5-txtimeTGD_residual.svg";
    const string tgdSvgFile = dataDir + "exam-5.5-txtimeTGD_tgd.svg";

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

    ObsReader obsReader;
    obsReader.setFileStream(&obsStream);
    obsReader.setHeader(&obsHeader);

    NavReader navReader;
    navReader.setFileStream(&navStream);
    navReader.setHeader(&navHeader);

    RinexNavStore navStore;
    size_t gpsRecordCount = 0;
    loadGpsNavIntoStore(navReader, navStore, gpsRecordCount);

    ofstream csv(csvFile);
    if (!csv) {
        cerr << "csv file open error: " << strerror(errno) << endl;
        return -1;
    }
    printCsvHeader(csv);

    TxStats tgdStats;
    TxStats ifStats;
    TxStats rawIfStats;
    TxStats c1Stats;
    TxStats rawResidualStats;
    TxStats tgdResidualStats;
    TxStats ifResidualStats;
    vector<TxRowSummary> plotRows;
    size_t printedRows = 0;
    size_t totalRows = 0;

    cout << "Exercise 5.5 txtimeTGD: transmit-time shift with TGD correction" << endl;
    cout << "Observation file: " << obsFile << endl;
    cout << "Navigation file : " << navFile << endl;
    cout << "Obs header marker: " << obsHeader.markerName << endl;
    cout << "Obs approx XYZ   : (" << obsHeader.approxX << ", " << obsHeader.approxY
         << ", " << obsHeader.approxZ << ")" << endl;
    cout << "Nav version/type : " << navHeader.version << " / " << navHeader.fileType << endl;
    cout << "GPS records loaded into RinexNavStore: " << gpsRecordCount << endl;
    cout << endl;

    while (true) {
        ObsEpochData epochData;
        try {
            epochData = obsReader.readNextEpoch();
        } catch (const EndOfObsFile &) {
            break;
        }

        CommonTime epochGps = toGpsTime(epochData.time);

        for (const auto &satEntry : epochData.satTypeValue) {
            const SatID sat(satEntry.first);
            if (sat.system != "G") {
                continue;
            }

            const TypeValueMap &tv = satEntry.second;
            auto itC1 = tv.find("C1C");
            auto itC2 = tv.find("C2W");
            if (itC1 == tv.end() || itC2 == tv.end()) {
                continue;
            }

            NavEphGPS eph = navStore.findGPSEph(sat, epochGps);
            if (!eph.isValid(epochGps)) {
                continue;
            }

            TxRowSummary row = buildRow(epochData, sat, tv, eph);
            if (!std::isfinite(row.ifPr)) {
                continue;
            }

            writeCsvRow(csv, row);
            ++totalRows;
            tgdStats.add(row.dtTgdNs);
            ifStats.add(row.dtIfNs);
            rawIfStats.add(row.dtRawIfNs);
            c1Stats.add(row.dc1Meters);
            rawResidualStats.add(row.c1RawResidual);
            tgdResidualStats.add(row.c1TgdResidual);
            ifResidualStats.add(row.ifResidual);
            if (row.sat.system == "G" && row.sat.id == plotGpsPrn) {
                plotRows.push_back(row);
            }

            if (printedRows < 12) {
                cout << "epoch=" << formatCivilTime(row.epochGps)
                     << " sat=" << row.sat
                     << " C1=" << fixed << setprecision(3) << row.c1Raw
                     << " C1_TGD=" << row.c1Tgd
                     << " IF=" << row.ifPr
                     << " TGD(s)=" << scientific << setprecision(12) << row.tgdSeconds
                     << " dC1(m)=" << row.dc1Meters
                     << " dTx_TGD(ns)=" << row.dtTgdNs
                     << " dTx_IF(ns)=" << row.dtIfNs
                     << " dTx_RAW_IF(ns)=" << row.dtRawIfNs
                     << " backC1=" << fixed << setprecision(3) << row.c1RawBack
                     << " backC1_TGD=" << row.c1TgdBack
                     << " backIF=" << row.ifBack
                     << " resC1(m)=" << row.c1RawResidual
                     << " resC1_TGD(m)=" << row.c1TgdResidual
                     << " resIF(m)=" << row.ifResidual
                     << endl;
                ++printedRows;
            }
        }
    }

    cout << "\n--- Summary ---" << endl;
    cout << "Rows written: " << totalRows << endl;
    cout << "C1 TGD correction(m): mean=" << c1Stats.mean()
         << " rms=" << c1Stats.rms()
         << " min=" << c1Stats.minVal
         << " max=" << c1Stats.maxVal << endl;
    cout << "TGD transmit-time shift(ns): mean=" << tgdStats.mean()
         << " rms=" << tgdStats.rms()
         << " min=" << tgdStats.minVal
         << " max=" << tgdStats.maxVal << endl;
    cout << "IF transmit-time shift(ns): mean=" << ifStats.mean()
         << " rms=" << ifStats.rms()
         << " min=" << ifStats.minVal
         << " max=" << ifStats.maxVal << endl;
    cout << "Raw vs IF transmit-time shift(ns): mean=" << rawIfStats.mean()
         << " rms=" << rawIfStats.rms()
         << " min=" << rawIfStats.minVal
         << " max=" << rawIfStats.maxVal << endl;
    cout << "Back-substitution residual(m):" << endl;
    cout << "  C1 raw  : mean=" << rawResidualStats.mean()
         << " rms=" << rawResidualStats.rms()
         << " min=" << rawResidualStats.minVal
         << " max=" << rawResidualStats.maxVal << endl;
    cout << "  C1 TGD  : mean=" << tgdResidualStats.mean()
         << " rms=" << tgdResidualStats.rms()
         << " min=" << tgdResidualStats.minVal
         << " max=" << tgdResidualStats.maxVal << endl;
    cout << "  IF      : mean=" << ifResidualStats.mean()
         << " rms=" << ifResidualStats.rms()
         << " min=" << ifResidualStats.minVal
         << " max=" << ifResidualStats.maxVal << endl;
    cout << "CSV written to: " << csvFile << endl;

    if (!plotRows.empty()) {
        SvgSeries resC1Raw{plotSatLabel + " raw residual", "#d62728", {}};
        SvgSeries resC1Tgd{plotSatLabel + " TGD residual", "#ff7f0e", {}};
        SvgSeries resIf{plotSatLabel + " IF residual", "#1f77b4", {}};
        SvgSeries c1TgdSeries{plotSatLabel + " C1 TGD correction", "#2ca02c", {}};

        double xMin = std::numeric_limits<double>::infinity();
        double xMax = -std::numeric_limits<double>::infinity();
        double residualYMin = std::numeric_limits<double>::infinity();
        double residualYMax = -std::numeric_limits<double>::infinity();
        double tgdYMin = std::numeric_limits<double>::infinity();
        double tgdYMax = -std::numeric_limits<double>::infinity();

        for (const auto &row : plotRows) {
            const double x = secondsOfDay(row.epochGps);
            xMin = std::min(xMin, x);
            xMax = std::max(xMax, x);
            resC1Raw.points.emplace_back(x, row.c1RawResidual);
            resC1Tgd.points.emplace_back(x, row.c1TgdResidual);
            resIf.points.emplace_back(x, row.ifResidual);
            c1TgdSeries.points.emplace_back(x, row.dc1Meters);

            residualYMin = std::min(residualYMin, row.c1RawResidual);
            residualYMax = std::max(residualYMax, row.c1RawResidual);
            residualYMin = std::min(residualYMin, row.c1TgdResidual);
            residualYMax = std::max(residualYMax, row.c1TgdResidual);
            residualYMin = std::min(residualYMin, row.ifResidual);
            residualYMax = std::max(residualYMax, row.ifResidual);

            tgdYMin = std::min(tgdYMin, row.dc1Meters);
            tgdYMax = std::max(tgdYMax, row.dc1Meters);
        }

        auto padRange = [](double &mn, double &mx, double ratio) {
            if (!std::isfinite(mn) || !std::isfinite(mx)) {
                mn = -1.0;
                mx = 1.0;
                return;
            }
            if (mx <= mn) {
                const double delta = (std::abs(mx) > 1.0 ? std::abs(mx) * ratio : 1.0);
                mn -= delta;
                mx += delta;
                return;
            }
            const double span = mx - mn;
            const double pad = span * ratio;
            mn -= pad;
            mx += pad;
        };

        padRange(tgdYMin, tgdYMax, 0.12);

        // Fix the residual plot to a symmetric +/-10 nm window so the
        // back-substitution check is easy to read.
        residualYMin = -10.0e-9;
        residualYMax = 10.0e-9;

        writeMultiSeriesChartSvg(residualSvgFile,
                                 "习题5.5：" + plotSatLabel + " 回代残差图",
                                 {resC1Raw, resC1Tgd, resIf},
                                 xMin, xMax,
                                 residualYMin, residualYMax);

        writeSingleChartSvg(tgdSvgFile,
                            "Exercise 5.5: " + plotSatLabel + " C1 TGD correction time series",
                            c1TgdSeries,
                            xMin, xMax,
                            tgdYMin, tgdYMax);
        cout << "Residual SVG written to: " << residualSvgFile << endl;
        cout << "TGD SVG written to: " << tgdSvgFile << endl;
    } else {
        cout << "SVG not written: no " << plotSatLabel << " samples collected for plotting." << endl;
    }

    obsStream.close();
    navStream.close();
    csv.close();

    return 0;
}
