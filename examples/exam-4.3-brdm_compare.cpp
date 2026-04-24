//
// Compare broadcast navigation ephemeris with precise SP3 products over one day.
// Each satellite is plotted as its own line, and BDS GEO vs IGSO/MEO statistics
// are summarized separately.
//

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <direct.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "NavEphBDS.hpp"
#include "NavEphGPS.hpp"
#include "SP3Store.hpp"
#include "TimeConvert.h"
#include "navreader.h"
#include "navreadheader.h"

#ifdef _WIN32
#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#endif

using namespace std;
using namespace Eigen;

namespace {

constexpr double kEpochStepSec = 30.0;
constexpr double kDayStartSec = 0.0;
constexpr double kDayEndSec = 86400.0 - kEpochStepSec;
constexpr int kGpsYear = 2025;
constexpr int kGpsMonth = 1;
constexpr int kGpsDay = 1;

struct SamplePoint {
    double secOfDay{0.0};
    Vector3d posDiff{Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())};
    Vector3d velDiff{Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())};
    bool valid{false};
};

struct SatSeries {
    SatID sat;
    bool geo{false};
    std::vector<SamplePoint> samples;
    size_t validCount{0};
    size_t missingCount{0};
};

struct VectorStats {
    size_t n{0};
    Vector3d sum{Vector3d::Zero()};
    Vector3d sumSq{Vector3d::Zero()};
    double normSum{0.0};
    double normSq{0.0};

    void add(const Vector3d& v) {
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2])) {
            return;
        }
        ++n;
        sum += v;
        sumSq += v.array().square().matrix();
        const double norm = v.norm();
        normSum += norm;
        normSq += norm * norm;
    }

    Vector3d mean() const {
        if (n == 0) return Vector3d::Zero();
        return sum / static_cast<double>(n);
    }

    Vector3d variance() const {
        if (n == 0) return Vector3d::Zero();
        const Vector3d m = mean();
        Vector3d var = sumSq / static_cast<double>(n) - m.array().square().matrix();
        for (int i = 0; i < 3; ++i) {
            if (var[i] < 0.0 && std::fabs(var[i]) < 1e-12) {
                var[i] = 0.0;
            }
        }
        return var;
    }

    double meanNorm() const {
        if (n == 0) return 0.0;
        return normSum / static_cast<double>(n);
    }

    double varianceNorm() const {
        if (n == 0) return 0.0;
        const double m = meanNorm();
        double var = normSq / static_cast<double>(n) - m * m;
        if (var < 0.0 && std::fabs(var) < 1e-12) var = 0.0;
        return var;
    }
};

struct GroupStats {
    VectorStats pos;
    VectorStats vel;
};

struct FigureLayout {
    int width{0};
    int height{0};
    int left{0};
    int right{0};
    int top{0};
    int bottom{0};
    int gap{0};
    int legendWidth{0};
    int plotW{0};
    int plotH{0};
    int legendX{0};
    int legendY{0};
    int legendH{0};
    int legendInnerW{0};
    int rowsPerCol{0};
};

static double safeComponentValue(const SamplePoint& pt, int component, bool velocity);

static string fmtVec(const Vector3d& v, int precision = 6) {
    ostringstream oss;
    oss << fixed << setprecision(precision)
        << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return oss.str();
}

static string fmtDouble(double v, int precision = 6) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << v;
    return oss.str();
}

static string fmtDoubleSci(double v, int precision = 12) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << v;
    return oss.str();
}

static string resolveExistingPath(const std::vector<std::string>& candidates) {
    for (const auto& path : candidates) {
        std::ifstream test(path.c_str(), std::ios::in);
        if (test.good()) {
            return path;
        }
    }
    throw std::runtime_error("Unable to locate any of the candidate data files.");
}

static void ensureDirectoryExists(const string& dir) {
    if (dir.empty()) return;
    string current;
    for (char ch : dir) {
        current.push_back(ch);
        if (ch == '\\' || ch == '/') {
            if (current.size() > 1 && current.back() != ':') {
                _mkdir(current.c_str());
            }
        }
    }
    _mkdir(dir.c_str());
}

static string joinPath(const string& base, const string& name) {
    if (base.empty()) return name;
    if (base.back() == '\\' || base.back() == '/') return base + name;
    return base + "\\" + name;
}

static string stripExtension(const string& path) {
    const size_t sep = path.find_last_of("/\\");
    const size_t dot = path.find_last_of('.');
    if (dot == string::npos || (sep != string::npos && dot < sep)) {
        return path;
    }
    return path.substr(0, dot);
}

static string openCsvWithFallback(ofstream& csv, const string& basePath) {
    const string stem = stripExtension(basePath);
    const string ext = ".csv";

    auto tryOpen = [&](const string& path) -> bool {
        csv.open(path.c_str(), ios::out | ios::trunc);
        return csv.is_open();
    };

    if (tryOpen(basePath)) {
        return basePath;
    }

    csv.clear();
    for (int i = 1; i <= 1000; ++i) {
        ostringstream oss;
        oss << stem << "_" << i << ext;
        const string candidate = oss.str();
        if (tryOpen(candidate)) {
            return candidate;
        }
        csv.clear();
    }

    throw runtime_error("Unable to open CSV output: " + basePath);
}

template<typename EphType>
static const EphType* selectBestEph(const std::map<CommonTime, EphType>& ephMap,
                                    const CommonTime& epoch) {
    if (ephMap.empty()) return nullptr;

    const EphType* bestValid = nullptr;
    double bestValidDt = std::numeric_limits<double>::infinity();
    const EphType* bestAny = nullptr;
    double bestAnyDt = std::numeric_limits<double>::infinity();

    for (const auto& kv : ephMap) {
        const EphType& eph = kv.second;
        const double dt = std::fabs(epoch - eph.ctToe);

        if (dt < bestAnyDt) {
            bestAnyDt = dt;
            bestAny = &eph;
        }
        if (eph.isValid(epoch) && dt < bestValidDt) {
            bestValidDt = dt;
            bestValid = &eph;
        }
    }

    return bestValid ? bestValid : bestAny;
}

static map<SatID, map<CommonTime, NavEphGPS> > loadGpsNavMap(const string& navPath) {
    fstream strm(navPath.c_str(), ios::in);
    if (!strm.is_open()) {
        throw runtime_error("Unable to open GPS navigation file: " + navPath);
    }

    NavReadHeader headerReader;
    headerReader.setFileStream(&strm);
    headerReader.parseHeader();

    NavReader reader;
    reader.setFileStream(&strm);
    reader.setHeader(&headerReader.getHeader());

    map<SatID, map<CommonTime, NavEphGPS> > navMap;
    while (true) {
        try {
            NavGpsRecord rec = reader.readNextGpsRecord();
            navMap[rec.sat][rec.eph.ctToe] = rec.eph;
        } catch (const EndOfNavFile&) {
            break;
        }
    }

    return navMap;
}

static map<SatID, map<CommonTime, NavEphBDS> > loadBdsNavMap(const string& navPath) {
    fstream strm(navPath.c_str(), ios::in);
    if (!strm.is_open()) {
        throw runtime_error("Unable to open BDS navigation file: " + navPath);
    }

    NavReadHeader headerReader;
    headerReader.setFileStream(&strm);
    headerReader.parseHeader();

    NavReader reader;
    reader.setFileStream(&strm);
    reader.setHeader(&headerReader.getHeader());

    map<SatID, map<CommonTime, NavEphBDS> > navMap;
    while (true) {
        try {
            NavBdsRecord rec = reader.readNextBdsRecord();
            navMap[rec.sat][rec.eph.ctToe] = rec.eph;
        } catch (const EndOfNavFile&) {
            break;
        }
    }

    return navMap;
}

static vector<SatID> intersectSats(const map<SatID, map<CommonTime, NavEphGPS> >& navMap,
                                   const SatIDSet& sp3Sats) {
    vector<SatID> sats;
    for (const auto& kv : navMap) {
        if (sp3Sats.count(kv.first) != 0) {
            sats.push_back(kv.first);
        }
    }
    sort(sats.begin(), sats.end());
    return sats;
}

static vector<SatID> intersectSats(const map<SatID, map<CommonTime, NavEphBDS> >& navMap,
                                   const SatIDSet& sp3Sats) {
    vector<SatID> sats;
    for (const auto& kv : navMap) {
        if (sp3Sats.count(kv.first) != 0) {
            sats.push_back(kv.first);
        }
    }
    sort(sats.begin(), sats.end());
    return sats;
}

static string satColor(size_t idx) {
    static const vector<string> palette = {
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600",
            "#2f4b7c", "#f95d6a", "#a05195", "#665191", "#3f8a89"
    };
    return palette[idx % palette.size()];
}

static double safeComponentValue(const SamplePoint& pt, int component, bool velocity) {
    const Vector3d& vec = velocity ? pt.velDiff : pt.posDiff;
    return vec[component];
}

static std::string componentLatin(int idx) {
    switch (idx) {
        case 0: return "X";
        case 1: return "Y";
        default: return "Z";
    }
}

static std::string componentSuffix() {
    return u8"分量";
}

static std::string legendLabel(const SatSeries& series) {
    return series.sat.toString();
}

static std::string vectorUnit(bool velocity) {
    return velocity ? u8"米/秒" : u8"米";
}

static std::string vectorStatUnit(bool velocity) {
    return velocity ? "(m/s)^2" : "m^2";
}

static FigureLayout makeFigureLayout(size_t seriesCount) {
    FigureLayout layout;
    layout.width = 2600;
    layout.height = 1800;
    layout.left = 220;
    layout.right = 80;
    layout.top = 170;
    layout.bottom = 140;
    layout.gap = 90;
    layout.legendWidth = static_cast<int>(std::round(layout.width * 0.20));
    layout.plotW = layout.width - layout.left - layout.right - layout.legendWidth - 25;
    layout.plotH = (layout.height - layout.top - layout.bottom - 2 * layout.gap) / 3;
    layout.legendX = layout.left + layout.plotW + 25;
    layout.legendY = layout.top;
    layout.legendInnerW = layout.legendWidth - 2;
    constexpr int legendCols = 2;
    layout.rowsPerCol = static_cast<int>((seriesCount + legendCols - 1) / legendCols);
    if (layout.rowsPerCol < 1) layout.rowsPerCol = 1;
    layout.legendH = 170 + layout.rowsPerCol * 50;
    return layout;
}

static std::string formatAxisLabel(double value) {
    ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << setprecision(6) << value;
    std::string s = oss.str();
    while (s.size() > 1 && s.find('.') != string::npos && s.back() == '0') {
        s.pop_back();
    }
    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }
    if (s == "-0") s = "0";
    return s;
}

static void computeRange(const vector<const SatSeries*>& seriesList,
                         bool velocity,
                         int comp,
                         double& yMin,
                         double& yMax) {
    yMin = std::numeric_limits<double>::infinity();
    yMax = -std::numeric_limits<double>::infinity();
    for (const SatSeries* series : seriesList) {
        for (const SamplePoint& pt : series->samples) {
            if (!pt.valid) continue;
            const double v = safeComponentValue(pt, comp, velocity);
            if (!std::isfinite(v)) continue;
            yMin = std::min(yMin, v);
            yMax = std::max(yMax, v);
        }
    }
    if (!std::isfinite(yMin) || !std::isfinite(yMax)) {
        yMin = -1.0;
        yMax = 1.0;
    }
    if (std::fabs(yMax - yMin) < 1e-9) {
        const double pad = std::max(1e-6, std::fabs(yMin) * 0.1 + 1e-3);
        yMin -= pad;
        yMax += pad;
    } else {
        const double pad = (yMax - yMin) * 0.08;
        yMin -= pad;
        yMax += pad;
    }
}

static double chooseNiceVelocityStep(double maxAbs) {
    if (maxAbs <= 0.015) return 0.005;
    if (maxAbs <= 0.03) return 0.01;
    if (maxAbs <= 0.08) return 0.02;
    if (maxAbs <= 0.2) return 0.05;
    if (maxAbs <= 0.5) return 0.1;
    return 0.2;
}

static void computeVelocityRange(const vector<const SatSeries*>& seriesList,
                                 double& yMin,
                                 double& yMax,
                                 double fixedAbsMax) {
    if (fixedAbsMax > 0.0) {
        yMin = -fixedAbsMax;
        yMax = fixedAbsMax;
        return;
    }

    double minVal = std::numeric_limits<double>::infinity();
    double maxVal = -std::numeric_limits<double>::infinity();
    for (const SatSeries* series : seriesList) {
        for (const SamplePoint& pt : series->samples) {
            if (!pt.valid) continue;
            for (int comp = 0; comp < 3; ++comp) {
                const double v = safeComponentValue(pt, comp, true);
                if (!std::isfinite(v)) continue;
                minVal = std::min(minVal, v);
                maxVal = std::max(maxVal, v);
            }
        }
    }

    if (!std::isfinite(minVal) || !std::isfinite(maxVal)) {
        yMin = -0.01;
        yMax = 0.01;
        return;
    }

    const double maxAbs = std::max(std::fabs(minVal), std::fabs(maxVal));
    const double step = chooseNiceVelocityStep(maxAbs);
    double upper = std::ceil(maxAbs / step) * step;
    if (upper < step) upper = step;
    yMin = -upper;
    yMax = upper;
}

static void writeSvgFigure(const string& outPath,
                           const string& title,
                           const string& subtitle,
                           const vector<const SatSeries*>& seriesList,
                           bool velocity,
                           double yMin,
                           double yMax) {
    const FigureLayout layout = makeFigureLayout(seriesList.size());
    const int width = layout.width;
    const int height = layout.height;
    const int left = layout.left;
    const int top = layout.top;
    const int bottom = layout.bottom;
    const int gap = layout.gap;
    const int plotH = layout.plotH;
    const int plotW = layout.plotW;
    const int legendX = layout.legendX;
    const int legendY = layout.legendY;
    const int legendW = layout.legendWidth;
    const int legendH = layout.legendH;
    const int rowsPerCol = layout.rowsPerCol;
    const int legendColW = (legendW - 40) / 2;
    const int legendRowStep = 50;

    auto sampleToX = [&](double sec) {
        return left + (sec / 86400.0) * plotW;
    };

    ofstream os(outPath.c_str(), ios::out | ios::trunc);
    if (!os.is_open()) {
        throw runtime_error("Unable to write SVG file: " + outPath);
    }

    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
       << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";
    os << "<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    if (title.rfind("GPS", 0) == 0) {
        os << "<text x=\"" << (width / 2) << "\" y=\"64\" text-anchor=\"middle\" font-size=\"49\" font-weight=\"700\" fill=\"#000\">"
           << "<tspan font-family=\"Times New Roman, serif\">GPS</tspan>"
           << "<tspan font-family=\"SimSun, serif\">" << title.substr(3) << "</tspan></text>\n";
    } else {
        os << "<text x=\"" << (width / 2) << "\" y=\"64\" text-anchor=\"middle\" font-family=\"SimSun, Times New Roman, serif\""
           << " font-size=\"49\" font-weight=\"700\" fill=\"#000\">" << title << "</text>\n";
    }
    if (!subtitle.empty()) {
        os << "<text x=\"" << (width / 2) << "\" y=\"112\" text-anchor=\"middle\" font-family=\"SimSun, Times New Roman, serif\""
           << " font-size=\"37\" fill=\"#000\">" << subtitle << "</text>\n";
    }
    os << "<rect x=\"" << legendX << "\" y=\"" << legendY << "\" width=\"" << legendW
       << "\" height=\"" << legendH << "\" rx=\"10\" ry=\"10\" fill=\"#fffdf9\" stroke=\"#d0d0d0\"/>\n";
    os << "<text x=\"" << (legendX + 14) << "\" y=\"" << (legendY + 24)
       << "\" font-family=\"SimSun, Times New Roman, serif\" font-size=\"36\" font-weight=\"700\" fill=\"#000\">图例</text>\n";

    for (int comp = 0; comp < 3; ++comp) {
        const int y0 = top + comp * (plotH + gap);
        const int plotY = y0 + 12;
        const int innerH = plotH - 12;

        const double plotTop = plotY;
        const double plotBottom = plotY + innerH;
        const double plotHeight = plotBottom - plotTop;

        os << "<rect x=\"" << left << "\" y=\"" << plotY << "\" width=\"" << plotW
           << "\" height=\"" << innerH << "\" fill=\"#ffffff\" stroke=\"#d0d0d0\"/>\n";

        os << "<text x=\"" << left + 6 << "\" y=\"" << (plotY + 34)
           << "\" font-family=\"Times New Roman, SimSun, serif\" font-size=\"36\" font-weight=\"600\" fill=\"#000\">"
           << "<tspan font-family=\"Times New Roman, serif\">" << componentLatin(comp) << "</tspan>"
           << "<tspan font-family=\"SimSun, serif\">" << componentSuffix() << "</tspan></text>\n";

        // Grid lines and y tick labels.
        for (int i = 0; i <= 4; ++i) {
            const double frac = static_cast<double>(i) / 4.0;
            const double y = plotBottom - frac * plotHeight;
            const double val = yMin + frac * (yMax - yMin);
            os << "<line x1=\"" << left << "\" y1=\"" << y << "\" x2=\"" << (left + plotW)
               << "\" y2=\"" << y << "\" stroke=\"#ececec\" stroke-width=\"1\"/>\n";
        os << "<text x=\"" << (left - 10) << "\" y=\"" << (y + 10)
               << "\" font-family=\"Times New Roman, SimSun, serif\" font-size=\"30\" fill=\"#000\" text-anchor=\"end\">"
               << formatAxisLabel(val) << "</text>\n";
        }

        for (int i = 0; i <= 4; ++i) {
            const double frac = static_cast<double>(i) / 4.0;
            const double x = left + frac * plotW;
            const double sec = frac * 86400.0;
            os << "<line x1=\"" << x << "\" y1=\"" << plotY << "\" x2=\"" << x
               << "\" y2=\"" << plotBottom << "\" stroke=\"#f0f0f0\" stroke-width=\"1\"/>\n";
            ostringstream tick;
            tick.setf(std::ios::fixed);
            tick << setprecision(0) << (sec / 3600.0);
            os << "<text x=\"" << x << "\" y=\"" << (plotBottom + 34)
               << "\" font-family=\"Times New Roman, SimSun, serif\" font-size=\"30\" fill=\"#000\" text-anchor=\"middle\">"
               << tick.str() << "</text>\n";
        }

        os << "<text x=\"" << (left + plotW / 2) << "\" y=\"" << (plotBottom + 78)
           << "\" font-family=\"SimSun, Times New Roman, serif\" font-size=\"32\" fill=\"#000\" text-anchor=\"middle\">"
           << u8"时间（小时）</text>\n";
        os << "<text x=\"" << 118 << "\" y=\"" << (plotTop + plotHeight / 2)
           << "\" transform=\"rotate(-90 118 " << (plotTop + plotHeight / 2) << ")\""
           << " font-family=\"SimSun, Times New Roman, serif\" font-size=\"32\" fill=\"#000\" text-anchor=\"middle\">"
           << u8"差值（" << vectorUnit(velocity) << u8"）</text>\n";

        // Draw each satellite line.
        for (size_t sidx = 0; sidx < seriesList.size(); ++sidx) {
            const SatSeries& series = *seriesList[sidx];
            string color = satColor(sidx);
            ostringstream path;
            path << fixed << setprecision(2);
            bool started = false;
            for (const SamplePoint& pt : series.samples) {
                if (!pt.valid) {
                    started = false;
                    continue;
                }
                const double val = safeComponentValue(pt, comp, velocity);
                if (!std::isfinite(val)) {
                    started = false;
                    continue;
                }
                const double x = sampleToX(pt.secOfDay);
                const double y = plotBottom - (val - yMin) / (yMax - yMin) * plotHeight;
                if (!started) {
                    path << "M " << x << " " << y << " ";
                    started = true;
                } else {
                    path << "L " << x << " " << y << " ";
                }
            }
            os << "<path d=\"" << path.str() << "\" fill=\"none\" stroke=\"" << color
               << "\" stroke-width=\"" << (series.geo ? 2.4 : 1.6) << "\" stroke-opacity=\"0.92\"";
            if (series.geo) {
                os << " stroke-dasharray=\"16,10\" stroke-linecap=\"round\"";
            }
            os << "/>\n";
        }

        // Panel frame on top of traces.
        os << "<rect x=\"" << left << "\" y=\"" << plotY << "\" width=\"" << plotW
           << "\" height=\"" << innerH << "\" fill=\"none\" stroke=\"#bbbbbb\" stroke-width=\"1\"/>\n";
    }

    for (size_t sidx = 0; sidx < seriesList.size(); ++sidx) {
        const SatSeries& series = *seriesList[sidx];
        const size_t col = sidx / static_cast<size_t>(rowsPerCol);
        const size_t row = sidx % static_cast<size_t>(rowsPerCol);
        const int lx = legendX + 20 + static_cast<int>(col) * legendColW;
        const int ly = legendY + 72 + static_cast<int>(row) * legendRowStep;
        const string color = satColor(sidx);
        os << "<line x1=\"" << lx << "\" y1=\"" << ly - 8 << "\" x2=\"" << (lx + 48)
           << "\" y2=\"" << ly - 8 << "\" stroke=\"" << color
           << "\" stroke-width=\"4\" stroke-opacity=\"0.95\"";
        if (series.geo) {
            os << " stroke-dasharray=\"18,12\" stroke-linecap=\"round\"";
        }
        os << "/>\n";
        os << "<text x=\"" << (lx + 80) << "\" y=\"" << ly
           << "\" font-family=\"Times New Roman, SimSun, serif\" font-size=\"30\" fill=\"#000\">"
           << legendLabel(series) << "</text>\n";
    }

    os << "</svg>\n";
}

#ifdef _WIN32
static int getPngEncoderClsid(CLSID* pClsid) {
    UINT num = 0;
    UINT size = 0;
    if (GetImageEncodersSize(&num, &size) != Ok || size == 0) return -1;
    std::vector<BYTE> buffer(size);
    ImageCodecInfo* pImageCodecInfo = reinterpret_cast<ImageCodecInfo*>(buffer.data());
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

static std::wstring widenPath(const std::string& path) {
    if (path.empty()) return std::wstring();
    int sizeNeeded = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, nullptr, 0);
    if (sizeNeeded <= 0) {
        sizeNeeded = MultiByteToWideChar(CP_ACP, 0, path.c_str(), -1, nullptr, 0);
    }
    if (sizeNeeded <= 0) {
        return std::wstring(path.begin(), path.end());
    }
    std::wstring wide(static_cast<size_t>(sizeNeeded - 1), L'\0');
    if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, &wide[0], sizeNeeded) <= 0) {
        MultiByteToWideChar(CP_ACP, 0, path.c_str(), -1, &wide[0], sizeNeeded);
    }
    return wide;
}

static void drawText(Graphics& g, const std::wstring& text, const Font& font,
                     const Brush& brush, float x, float y, StringFormat* format = nullptr) {
    RectF rect(x, y, 1000.0f, 100.0f);
    g.DrawString(text.c_str(), static_cast<INT>(text.size()), &font, rect, format, &brush);
}

static void drawMixedComponentTitle(Graphics& g,
                                    const std::wstring& latin,
                                    const std::wstring& cjk,
                                    const Font& latinFont,
                                    const Font& cjkFont,
                                    const Brush& brush,
                                    float x,
                                    float y) {
    g.DrawString(latin.c_str(), static_cast<INT>(latin.size()), &latinFont, PointF(x, y), &brush);
    g.DrawString(cjk.c_str(), static_cast<INT>(cjk.size()), &cjkFont, PointF(x + 28.0f, y), &brush);
}

static void drawDashedLine(Graphics& g,
                           const Color& color,
                           float width,
                           const PointF& a,
                           const PointF& b,
                           float dashLen = 14.0f,
                           float gapLen = 10.0f) {
    const float dx = b.X - a.X;
    const float dy = b.Y - a.Y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len <= 1e-4f) {
        return;
    }

    Pen pen(color, width);
    pen.SetLineJoin(LineJoinRound);

    const float step = dashLen + gapLen;
    for (float pos = 0.0f; pos < len; pos += step) {
        const float segStart = pos;
        const float segEnd = std::min(pos + dashLen, len);
        const float t0 = segStart / len;
        const float t1 = segEnd / len;
        PointF p0(a.X + dx * t0, a.Y + dy * t0);
        PointF p1(a.X + dx * t1, a.Y + dy * t1);
        g.DrawLine(&pen, p0, p1);
    }
}

static float measureTextWidth(Graphics& g, const std::wstring& text, const Font& font) {
    RectF bounds;
    g.MeasureString(text.c_str(), static_cast<INT>(text.size()), &font, PointF(0.0f, 0.0f), &bounds);
    return bounds.Width;
}

static void writePngFigure(const string& outPath,
                           const string& title,
                           const string& subtitle,
                           const vector<const SatSeries*>& seriesList,
                           bool velocity,
                           double yMin,
                           double yMax) {
    const FigureLayout layout = makeFigureLayout(seriesList.size());
    const int width = layout.width;
    const int height = layout.height;
    const int left = layout.left;
    const int top = layout.top;
    const int bottom = layout.bottom;
    const int gap = layout.gap;
    const int plotH = layout.plotH;
    const int plotW = layout.plotW;
    const int legendX = layout.legendX;
    const int legendY = layout.legendY;
    const int legendW = layout.legendWidth;
    const int legendH = layout.legendH;
    const int rowsPerCol = layout.rowsPerCol;
    const int legendColW = (legendW - 40) / 2;
    const int legendRowStep = 50;

    Bitmap bitmap(width, height, PixelFormat32bppARGB);
    Graphics graphics(&bitmap);
    graphics.SetSmoothingMode(SmoothingModeHighQuality);
    graphics.SetTextRenderingHint(TextRenderingHintClearTypeGridFit);

    SolidBrush white(Color(255, 255, 255, 255));
    SolidBrush black(Color(255, 17, 17, 17));
    SolidBrush gray(Color(255, 85, 85, 85));
    SolidBrush panelFill(Color(255, 255, 253, 249));
    Pen borderPen(Color(255, 208, 208, 208), 1.0f);
    Pen gridPen(Color(255, 236, 236, 236), 1.0f);
    Pen gridPen2(Color(255, 240, 240, 240), 1.0f);
    Pen panelPen(Color(255, 187, 187, 187), 1.0f);
    Pen legendPen(Color(255, 208, 208, 208), 1.0f);

    FontFamily times(L"Times New Roman");
    FontFamily simsun(L"SimSun");
    Font titleLatinFont(&times, 49.0f, FontStyleBold, UnitPixel);
    Font titleCjkFont(&simsun, 49.0f, FontStyleBold, UnitPixel);
    Font subFont(&simsun, 37.0f, FontStyleRegular, UnitPixel);
    Font panelLatinFont(&times, 41.0f, FontStyleRegular, UnitPixel);
    Font panelCjkFont(&simsun, 41.0f, FontStyleBold, UnitPixel);
    Font legendTitleFont(&simsun, 41.0f, FontStyleBold, UnitPixel);
    Font legendFont(&times, 35.0f, FontStyleRegular, UnitPixel);
    Font tickFont(&times, 35.0f, FontStyleRegular, UnitPixel);
    StringFormat centerFmt;
    centerFmt.SetAlignment(StringAlignmentCenter);
    centerFmt.SetLineAlignment(StringAlignmentCenter);
    StringFormat rightFmt;
    rightFmt.SetAlignment(StringAlignmentFar);
    rightFmt.SetLineAlignment(StringAlignmentCenter);

    graphics.FillRectangle(&white, 0, 0, width, height);
    const std::wstring wTitle = widenPath(title);
    if (title.rfind("GPS", 0) == 0) {
        const std::wstring wGps = L"GPS";
        const std::wstring wRest = widenPath(title.substr(3));
        const float gpsWidth = measureTextWidth(graphics, wGps, titleLatinFont);
        const float restWidth = measureTextWidth(graphics, wRest, titleCjkFont);
        const float startX = static_cast<float>((width - (gpsWidth + restWidth)) / 2.0);
        graphics.DrawString(wGps.c_str(), static_cast<INT>(wGps.size()), &titleLatinFont,
                            PointF(startX, 12.0f), &black);
        graphics.DrawString(wRest.c_str(), static_cast<INT>(wRest.size()), &titleCjkFont,
                            PointF(startX + gpsWidth, 12.0f), &black);
    } else {
        RectF titleRect(0.0f, 12.0f, static_cast<REAL>(width), 70.0f);
        graphics.DrawString(wTitle.c_str(), static_cast<INT>(wTitle.size()), &titleCjkFont,
                            titleRect, &centerFmt, &black);
    }
    if (!subtitle.empty()) {
        const std::wstring wSubtitle = widenPath(subtitle);
        RectF subtitleRect(0.0f, 72.0f, static_cast<REAL>(width), 48.0f);
        graphics.DrawString(wSubtitle.c_str(), static_cast<INT>(wSubtitle.size()), &subFont,
                            subtitleRect, &centerFmt, &black);
    }
    graphics.FillRectangle(&panelFill, legendX, legendY, legendW, legendH);
    graphics.DrawRectangle(&legendPen, legendX, legendY, legendW, legendH);
    graphics.DrawString(L"图例", -1, &legendTitleFont, PointF(static_cast<REAL>(legendX + 14), static_cast<REAL>(legendY + 10)), &black);

    for (int comp = 0; comp < 3; ++comp) {
        const int y0 = top + comp * (plotH + gap);
        const int plotY = y0 + 12;
        const int innerH = plotH - 12;
        const int plotBottom = plotY + innerH;
        const int plotHeight = plotBottom - plotY;

        graphics.FillRectangle(&white, left, plotY, plotW, innerH);
        graphics.DrawRectangle(&borderPen, left, plotY, plotW, innerH);

        drawMixedComponentTitle(graphics,
                                widenPath(componentLatin(comp)),
                                widenPath(componentSuffix()),
                                panelLatinFont,
                                panelCjkFont,
                                black,
                                static_cast<REAL>(left + 6),
                                static_cast<REAL>(plotY + 4));

        for (int i = 0; i <= 4; ++i) {
            const float frac = static_cast<float>(i) / 4.0f;
            const float y = static_cast<float>(plotBottom - frac * plotHeight);
            const double val = yMin + frac * (yMax - yMin);
            graphics.DrawLine(&gridPen, static_cast<REAL>(left), y, static_cast<REAL>(left + plotW), y);
            const std::wstring ws = widenPath(formatAxisLabel(val));
            graphics.DrawString(ws.c_str(), static_cast<INT>(ws.size()), &tickFont,
                                RectF(static_cast<REAL>(left - 170), y - 18.0f, 160.0f, 40.0f), &rightFmt, &black);
        }

        for (int i = 0; i <= 4; ++i) {
            const float frac = static_cast<float>(i) / 4.0f;
            const float x = static_cast<float>(left + frac * plotW);
            const double sec = frac * 86400.0;
            graphics.DrawLine(&gridPen2, x, static_cast<REAL>(plotY), x, static_cast<REAL>(plotBottom));
            std::wostringstream oss;
            oss.setf(std::ios::fixed);
            oss << setprecision(0) << (sec / 3600.0);
            const std::wstring ws = oss.str();
            graphics.DrawString(ws.c_str(), static_cast<INT>(ws.size()), &tickFont,
                                RectF(x - 52.0f, static_cast<REAL>(plotBottom + 8), 104.0f, 38.0f), &centerFmt, &black);
        }

        graphics.DrawString(L"时间（小时）", -1, &subFont,
                            PointF(static_cast<REAL>(left + plotW / 2 - 90), static_cast<REAL>(plotBottom + 58)), &black);

        std::wstring yLabel = L"差值（";
        yLabel += widenPath(vectorUnit(velocity));
        yLabel += L"）";
        GraphicsState state = graphics.Save();
        graphics.TranslateTransform(72.0f, static_cast<REAL>(plotY + plotHeight / 2 + 102));
        graphics.RotateTransform(-90.0f);
        graphics.DrawString(yLabel.c_str(), static_cast<INT>(yLabel.size()), &subFont,
                            PointF(0.0f, 0.0f), &black);
        graphics.Restore(state);

        for (size_t sidx = 0; sidx < seriesList.size(); ++sidx) {
            const SatSeries& series = *seriesList[sidx];
            const Color color = Color(255,
                                      static_cast<BYTE>(std::stoi(satColor(sidx).substr(1, 2), nullptr, 16)),
                                      static_cast<BYTE>(std::stoi(satColor(sidx).substr(3, 2), nullptr, 16)),
                                      static_cast<BYTE>(std::stoi(satColor(sidx).substr(5, 2), nullptr, 16)));
            const float lineWidth = series.geo ? 2.6f : 1.8f;
            std::vector<PointF> runPoints;
            for (const SamplePoint& pt : series.samples) {
                if (!pt.valid) {
                    if (runPoints.size() >= 2) {
                        if (series.geo) {
                            GraphicsPath path;
                            path.AddLines(runPoints.data(), static_cast<INT>(runPoints.size()));
                            Pen linePen(color, lineWidth);
                            linePen.SetLineJoin(LineJoinRound);
                            linePen.SetDashStyle(DashStyleDash);
                            graphics.DrawPath(&linePen, &path);
                        } else {
                            Pen linePen(color, lineWidth);
                            linePen.SetLineJoin(LineJoinRound);
                            graphics.DrawLines(&linePen, runPoints.data(), static_cast<INT>(runPoints.size()));
                        }
                    }
                    runPoints.clear();
                    continue;
                }
                const double val = safeComponentValue(pt, comp, velocity);
                if (!std::isfinite(val)) {
                    if (runPoints.size() >= 2) {
                        if (series.geo) {
                            GraphicsPath path;
                            path.AddLines(runPoints.data(), static_cast<INT>(runPoints.size()));
                            Pen linePen(color, lineWidth);
                            linePen.SetLineJoin(LineJoinRound);
                            linePen.SetDashStyle(DashStyleDash);
                            graphics.DrawPath(&linePen, &path);
                        } else {
                            Pen linePen(color, lineWidth);
                            linePen.SetLineJoin(LineJoinRound);
                            graphics.DrawLines(&linePen, runPoints.data(), static_cast<INT>(runPoints.size()));
                        }
                    }
                    runPoints.clear();
                    continue;
                }
                const float x = static_cast<float>(left + (pt.secOfDay / 86400.0) * plotW);
                const float y = static_cast<float>(plotBottom - (val - yMin) / (yMax - yMin) * plotHeight);
                runPoints.emplace_back(x, y);
            }
            if (runPoints.size() >= 2) {
                if (series.geo) {
                    GraphicsPath path;
                    path.AddLines(runPoints.data(), static_cast<INT>(runPoints.size()));
                    Pen linePen(color, lineWidth);
                    linePen.SetLineJoin(LineJoinRound);
                    linePen.SetDashStyle(DashStyleDash);
                    graphics.DrawPath(&linePen, &path);
                } else {
                    Pen linePen(color, lineWidth);
                    linePen.SetLineJoin(LineJoinRound);
                    graphics.DrawLines(&linePen, runPoints.data(), static_cast<INT>(runPoints.size()));
                }
            }
        }

    }

    for (size_t sidx = 0; sidx < seriesList.size(); ++sidx) {
        const SatSeries& series = *seriesList[sidx];
        const size_t col = sidx / static_cast<size_t>(rowsPerCol);
        const size_t row = sidx % static_cast<size_t>(rowsPerCol);
        const int lx = legendX + 14 + static_cast<int>(col) * legendColW;
        const int ly = legendY + 72 + static_cast<int>(row) * legendRowStep;
        const string colorHex = satColor(sidx);
        const Color color = Color(255,
                                  static_cast<BYTE>(std::stoi(colorHex.substr(1, 2), nullptr, 16)),
                                  static_cast<BYTE>(std::stoi(colorHex.substr(3, 2), nullptr, 16)),
                                  static_cast<BYTE>(std::stoi(colorHex.substr(5, 2), nullptr, 16)));
        const float lineWidth = 5.0f;
        if (series.geo) {
            drawDashedLine(graphics, color, lineWidth, PointF(static_cast<REAL>(lx), static_cast<REAL>(ly - 8)),
                           PointF(static_cast<REAL>(lx + 48), static_cast<REAL>(ly - 8)), 16.0f, 10.0f);
        } else {
            Pen linePen(color, lineWidth);
            graphics.DrawLine(&linePen, lx, ly - 8, lx + 48, ly - 8);
        }
        std::wstring ws = widenPath(legendLabel(series));
        graphics.DrawString(ws.c_str(), static_cast<INT>(ws.size()), &legendFont,
                            PointF(static_cast<REAL>(lx + 82), static_cast<REAL>(ly - 20)), &black);
    }

    CLSID pngClsid{};
    if (getPngEncoderClsid(&pngClsid) < 0) {
        throw runtime_error("PNG encoder not found in GDI+.");
    }
    const std::wstring wpath = widenPath(outPath);
    if (bitmap.Save(wpath.c_str(), &pngClsid, nullptr) != Ok) {
        throw runtime_error("Failed to save PNG: " + outPath);
    }
}
#endif

static void writeCsvHeader(ofstream& csv) {
    csv << "system,satellite,class,sec_of_day,pos_dx,pos_dy,pos_dz,vel_dx,vel_dy,vel_dz\n";
}

static void appendCsvRow(ofstream& csv,
                         const string& system,
                         const SatID& sat,
                         const string& satClass,
                         const SamplePoint& pt) {
    csv << system << ","
        << sat.toString() << ","
        << satClass << ","
        << fixed << setprecision(0) << pt.secOfDay << ","
        << fixed << setprecision(12)
        << pt.posDiff[0] << ","
        << pt.posDiff[1] << ","
        << pt.posDiff[2] << ","
        << pt.velDiff[0] << ","
        << pt.velDiff[1] << ","
        << pt.velDiff[2] << "\n";
}

static void printVectorStats(const string& label, const VectorStats& stats, bool velocity) {
    cout << label << " samples: " << stats.n << "\n";
    if (velocity) {
        cout << "  mean velocity bias: " << fmtVec(stats.mean(), 12) << " m/s\n";
        cout << "  variance: " << fmtVec(stats.variance(), 12) << " (m/s)^2\n";
        cout << "  mean norm: " << fmtDoubleSci(stats.meanNorm(), 12) << "\n";
        cout << "  variance norm: " << fmtDoubleSci(stats.varianceNorm(), 12) << "\n";
    } else {
        cout << "  mean position bias: " << fmtVec(stats.mean(), 4) << " m\n";
        cout << "  variance: " << fmtVec(stats.variance(), 4) << " m^2\n";
        cout << "  mean norm: " << fmtDouble(stats.meanNorm(), 4) << "\n";
        cout << "  variance norm: " << fmtDouble(stats.varianceNorm(), 4) << "\n";
    }
}

static void printGroupStats(const string& label, const GroupStats& stats) {
    cout << "\n[" << label << "]\n";
    printVectorStats("Position", stats.pos, false);
    printVectorStats("Velocity", stats.vel, true);
}

static void evaluateGpsDay(const map<SatID, map<CommonTime, NavEphGPS> >& navMap,
                           SP3Store& sp3,
                           const CommonTime& startGps,
                           vector<SatSeries>& seriesOut,
                           ofstream& csv) {
    const SatIDSet sp3Sats = sp3.getSatSet();
    const vector<SatID> sats = intersectSats(navMap, sp3Sats);

    cout << "GPS satellites with both nav and SP3 data: " << sats.size() << "\n";

    for (const SatID& sat : sats) {
        SatSeries series;
        series.sat = sat;
        series.geo = false;
        series.samples.reserve(static_cast<size_t>((kDayEndSec - kDayStartSec) / kEpochStepSec) + 1);

        const auto navIt = navMap.find(sat);
        if (navIt == navMap.end()) continue;

        for (double sec = kDayStartSec; sec <= kDayEndSec + 0.5; sec += kEpochStepSec) {
            SamplePoint pt;
            pt.secOfDay = sec;

            CommonTime epochGps = startGps + sec;
            epochGps.setTimeSystem(TimeSystem::GPS);

            try {
                const NavEphGPS* eph = selectBestEph(navIt->second, epochGps);
                if (!eph) throw runtime_error("No GPS ephemeris found.");

                Xvt navXvt = eph->svXvt(epochGps);
                Xvt sp3Xvt = sp3.getXvt(sat, epochGps);

                pt.posDiff = navXvt.x - sp3Xvt.x;
                pt.velDiff = navXvt.v - sp3Xvt.v;
                pt.valid = true;

                ++series.validCount;
                appendCsvRow(csv, "GPS", sat, "GPS", pt);
            } catch (const exception& e) {
                pt.valid = false;
                ++series.missingCount;
                pt.posDiff = Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
                pt.velDiff = Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
            }

            series.samples.push_back(pt);
        }

        seriesOut.push_back(series);
    }
}

static void evaluateBdsDay(const map<SatID, map<CommonTime, NavEphBDS> >& navMap,
                           SP3Store& sp3,
                           const CommonTime& startGps,
                           vector<SatSeries>& seriesOut,
                           GroupStats& geoStats,
                           GroupStats& meoStats,
                           ofstream& csv) {
    const SatIDSet sp3Sats = sp3.getSatSet();
    const vector<SatID> sats = intersectSats(navMap, sp3Sats);

    cout << "BDS satellites with both nav and SP3 data: " << sats.size() << "\n";

    for (const SatID& sat : sats) {
        SatSeries series;
        series.sat = sat;
        series.geo = NavEphBDS::isGeoSatellite(sat);
        series.samples.reserve(static_cast<size_t>((kDayEndSec - kDayStartSec) / kEpochStepSec) + 1);

        const auto navIt = navMap.find(sat);
        if (navIt == navMap.end()) continue;

        for (double sec = kDayStartSec; sec <= kDayEndSec + 0.5; sec += kEpochStepSec) {
            SamplePoint pt;
            pt.secOfDay = sec;

            CommonTime epochGps = startGps + sec;
            epochGps.setTimeSystem(TimeSystem::GPS);
            CommonTime epochBdt = convertTimeSystem(epochGps, TimeSystem::BDT);
            epochBdt.setTimeSystem(TimeSystem::BDT);

            try {
                const NavEphBDS* eph = selectBestEph(navIt->second, epochBdt);
                if (!eph) throw runtime_error("No BDS ephemeris found.");

                Xvt navXvt = eph->svXvt(sat, epochBdt);
                Xvt sp3Xvt = sp3.getXvt(sat, epochGps);

                pt.posDiff = navXvt.x - sp3Xvt.x;
                pt.velDiff = navXvt.v - sp3Xvt.v;
                pt.valid = true;

                ++series.validCount;
                appendCsvRow(csv, "BDS", sat, series.geo ? "GEO" : "IGSO/MEO", pt);

                if (series.geo) {
                    geoStats.pos.add(pt.posDiff);
                    geoStats.vel.add(pt.velDiff);
                } else {
                    meoStats.pos.add(pt.posDiff);
                    meoStats.vel.add(pt.velDiff);
                }
            } catch (const exception& e) {
                pt.valid = false;
                ++series.missingCount;
                pt.posDiff = Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
                pt.velDiff = Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
            }

            series.samples.push_back(pt);
        }

        seriesOut.push_back(series);
    }
}

static void printComparisonSummary(const GroupStats& geoStats, const GroupStats& meoStats) {
    cout << "\n[BDS GEO vs IGSO/MEO comparison]\n";
    cout << "Position mean norm GEO        : " << fmtDouble(geoStats.pos.meanNorm(), 6) << "\n";
    cout << "Position mean norm IGSO/MEO   : " << fmtDouble(meoStats.pos.meanNorm(), 6) << "\n";
    cout << "Position variance norm GEO    : " << fmtDouble(geoStats.pos.varianceNorm(), 6) << "\n";
    cout << "Position variance norm IGSO/MEO: " << fmtDouble(meoStats.pos.varianceNorm(), 6) << "\n";
    cout << "Velocity mean norm GEO        : " << fmtDoubleSci(geoStats.vel.meanNorm(), 12) << "\n";
    cout << "Velocity mean norm IGSO/MEO   : " << fmtDoubleSci(meoStats.vel.meanNorm(), 12) << "\n";
    cout << "Velocity variance norm GEO    : " << fmtDoubleSci(geoStats.vel.varianceNorm(), 12) << "\n";
    cout << "Velocity variance norm IGSO/MEO: " << fmtDoubleSci(meoStats.vel.varianceNorm(), 12) << "\n";

    cout << "Interpretation: compare the component-wise means and variances above. "
         << "GEO satellites are expected to differ from IGSO/MEO mainly because of the GEO-specific "
         << "rotation treatment in the BDS broadcast model." << "\n";
}

} // namespace

int main(int argc, char* argv[]) {
    try {
#ifdef _WIN32
        GdiplusSession gdiplusSession;
#endif
        const string navPath = resolveExistingPath({
                argc > 1 ? argv[1] : string("data/BRDC00IGS_R_20250010000_01D_MN.rnx"),
                argc > 1 ? argv[1] : string("../data/BRDC00IGS_R_20250010000_01D_MN.rnx"),
                argc > 1 ? argv[1] : string("../../data/BRDC00IGS_R_20250010000_01D_MN.rnx")
        });
        const string sp3Path = resolveExistingPath({
                argc > 2 ? argv[2] : string("data/WUM0MGXFIN_20250010000_01D_05M_ORB.SP3"),
                argc > 2 ? argv[2] : string("../data/WUM0MGXFIN_20250010000_01D_05M_ORB.SP3"),
                argc > 2 ? argv[2] : string("../../data/WUM0MGXFIN_20250010000_01D_05M_ORB.SP3")
        });
        const string outDir = (argc > 4)
                              ? argv[4]
                              : string(R"(D:\GNSSLAB\gnssLab-2.4\result\nav_sp3_compare)");

        cout << "Navigation file : " << navPath << "\n";
        cout << "SP3 file        : " << sp3Path << "\n";
        cout << "Output dir      : " << outDir << "\n";

        ensureDirectoryExists(outDir);

        SP3Store sp3;
        sp3.loadSP3File(sp3Path);
        cout << "SP3 satellites  : " << sp3.getSatSet().size() << "\n";

        const auto gpsNavMap = loadGpsNavMap(navPath);
        const auto bdsNavMap = loadBdsNavMap(navPath);

        CivilTime startCivil(kGpsYear, kGpsMonth, kGpsDay, 0, 0, 0.0, TimeSystem::GPS);
        CommonTime startGps = CivilTime2CommonTime(startCivil);
        startGps.setTimeSystem(TimeSystem::GPS);

        vector<SatSeries> gpsSeries;
        vector<SatSeries> bdsSeries;
        GroupStats bdsGeoStats;
        GroupStats bdsMeoStats;

        const string csvBasePath = joinPath(outDir, "brdm_compare_20250101_30s_satellite_diff.csv");
        ofstream csv;
        const string csvPath = openCsvWithFallback(csv, csvBasePath);
        writeCsvHeader(csv);

        evaluateGpsDay(gpsNavMap, sp3, startGps, gpsSeries, csv);
        evaluateBdsDay(bdsNavMap, sp3, startGps, bdsSeries, bdsGeoStats, bdsMeoStats, csv);

        sort(gpsSeries.begin(), gpsSeries.end(), [](const SatSeries& a, const SatSeries& b) {
            return a.sat < b.sat;
        });
        sort(bdsSeries.begin(), bdsSeries.end(), [](const SatSeries& a, const SatSeries& b) {
            return a.sat < b.sat;
        });

        vector<const SatSeries*> gpsPtrs;
        for (const auto& series : gpsSeries) gpsPtrs.push_back(&series);
        vector<const SatSeries*> bdsPtrs;
        for (const auto& series : bdsSeries) bdsPtrs.push_back(&series);

        double gpsVelMin = 0.0;
        double gpsVelMax = 0.0;
        computeVelocityRange(gpsPtrs, gpsVelMin, gpsVelMax, 0.20);

        double bdsVelMin = 0.0;
        double bdsVelMax = 0.0;
        computeVelocityRange(bdsPtrs, bdsVelMin, bdsVelMax, 0.0005);

        writeSvgFigure(joinPath(outDir, "brdm_compare_20250101_30s_gps_pos_diff.svg"),
                       u8"GPS广播星历与精密星历位置差",
                       "",
                       gpsPtrs,
                       false,
                       -10.0,
                       10.0);
        writeSvgFigure(joinPath(outDir, "brdm_compare_20250101_30s_gps_vel_diff.svg"),
                       u8"GPS广播星历与精密星历速度差",
                       "",
                       gpsPtrs,
                       true,
                       gpsVelMin,
                       gpsVelMax);
        writeSvgFigure(joinPath(outDir, "brdm_compare_20250101_30s_bds_pos_diff.svg"),
                       u8"北斗广播星历与精密星历位置差",
                       "",
                       bdsPtrs,
                       false,
                       -13.0,
                       13.0);
        writeSvgFigure(joinPath(outDir, "brdm_compare_20250101_30s_bds_vel_diff.svg"),
                       u8"北斗广播星历与精密星历速度差",
                       "",
                       bdsPtrs,
                       true,
                       bdsVelMin,
                       bdsVelMax);

#ifdef _WIN32
        writePngFigure(joinPath(outDir, "brdm_compare_20250101_30s_gps_pos_diff.png"),
                       u8"GPS广播星历与精密星历位置差",
                       "",
                       gpsPtrs,
                       false,
                       -10.0,
                       10.0);
        writePngFigure(joinPath(outDir, "brdm_compare_20250101_30s_gps_vel_diff.png"),
                       u8"GPS广播星历与精密星历速度差",
                       "",
                       gpsPtrs,
                       true,
                       gpsVelMin,
                       gpsVelMax);
        writePngFigure(joinPath(outDir, "brdm_compare_20250101_30s_bds_pos_diff.png"),
                       u8"北斗广播星历与精密星历位置差",
                       "",
                       bdsPtrs,
                       false,
                       -13.0,
                       13.0);
        writePngFigure(joinPath(outDir, "brdm_compare_20250101_30s_bds_vel_diff.png"),
                       u8"北斗广播星历与精密星历速度差",
                       "",
                       bdsPtrs,
                       true,
                       bdsVelMin,
                       bdsVelMax);
#endif

        cout << "\nDaily CSV written to: " << csvPath << "\n";
        cout << "SVG figures written to:\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_gps_pos_diff.svg") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_gps_vel_diff.svg") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_bds_pos_diff.svg") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_bds_vel_diff.svg") << "\n";
        cout << "PNG figures written to:\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_gps_pos_diff.png") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_gps_vel_diff.png") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_bds_pos_diff.png") << "\n";
        cout << "  " << joinPath(outDir, "brdm_compare_20250101_30s_bds_vel_diff.png") << "\n";

        size_t gpsValid = 0, gpsMissing = 0, bdsValid = 0, bdsMissing = 0;
        for (const auto& s : gpsSeries) {
            gpsValid += s.validCount;
            gpsMissing += s.missingCount;
        }
        for (const auto& s : bdsSeries) {
            bdsValid += s.validCount;
            bdsMissing += s.missingCount;
        }
        cout << "\nCoverage summary:\n";
        cout << "  GPS valid samples  : " << gpsValid << ", missing: " << gpsMissing << "\n";
        cout << "  BDS valid samples  : " << bdsValid << ", missing: " << bdsMissing << "\n";

        printGroupStats("BDS GEO", bdsGeoStats);
        printGroupStats("BDS IGSO/MEO", bdsMeoStats);
        printComparisonSummary(bdsGeoStats, bdsMeoStats);

        if (!gpsSeries.empty()) {
            cout << "\nGPS satellites plotted individually: " << gpsSeries.size() << "\n";
        }
        if (!bdsSeries.empty()) {
            cout << "BDS satellites plotted individually: " << bdsSeries.size() << "\n";
        }

        return 0;
    } catch (const std::exception& e) {
        cerr << "brdm_compare failed: " << e.what() << endl;
        return 1;
    }
}
