/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: Shoujian Zhang，shjzhang@sgg.whu.edu.cn， 2024-10-10
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernández-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */


#include <string>
#include <algorithm> //replace 函数
#include <cmath>
#include <cctype>
#include "TimeConvert.h"
#include "GnssStruct.h"
#include "GnssFunc.h"
#include "ARLambda.hpp"
#include "CoordConvert.h"
#include "RinexNavStore.hpp"

#define debug 1
#define debugCSMW 1

void parseRinexHeader(std::fstream &rinexFileStream, RinexHeader &rinexHeader) {

    double version;
    XYZ antennaPosition;
    string satSys;
    std::map<string, std::vector<string>> mapObsTypes;
    while (true) {
        string line;
        getline(rinexFileStream, line);

        cout << "parseRinexHeader:" << line << endl;

        string label;
        if (line.size() >= 80)
            label = line.substr(60, 20);

        strip(label);

        if (label == "END OF HEADER") {
            break;
        } else if (label == "MARKER NAME") {
            string markerName = line.substr(0, 60);
            std::replace(markerName.begin(), markerName.end(), ' ', '_');
            rinexHeader.station = markerName;

        } else if (label == "RINEX VERSION / TYPE") {
            version = safeStod(line.substr(0, 20));
            if (version != 3.04) {
                cerr << "only support rinex 3.04 version!" << endl;
                exit(-1);
            }
            rinexHeader.version = version;
        } else if (label == "APPROX POSITION XYZ") {
            antennaPosition[0] = safeStod(line.substr(0, 14));
            antennaPosition[1] = safeStod(line.substr(14, 14));
            antennaPosition[2] = safeStod(line.substr(28, 14));
            rinexHeader.antennaPosition = antennaPosition;
        } else if (label == "SYS / # / OBS TYPES") {
            string sysStr;
            sysStr = line.substr(0, 1);
            strip(sysStr);

            int numObs;

            if (sysStr != "") {
                numObs = stoi(line.substr(3, 3));
                satSys = sysStr;
            }

            const int maxObsPerLine = 13;
            for (int i = 0; i < maxObsPerLine && mapObsTypes[satSys].size() < numObs; i++) {
                std::string typeStr = (line.substr(4 * i + 7, 3));
                // insert into mapObsTypes
                mapObsTypes[satSys].push_back(typeStr);
            }
            rinexHeader.mapObsTypes = mapObsTypes;
        }
    }
};

ObsData parseRinexObs(std::fstream &rinexFileStream) {
    static bool isHeaderRead = false;
    static RinexHeader rinexHeader;

    if (!isHeaderRead) {
        parseRinexHeader(rinexFileStream, rinexHeader);
        isHeaderRead = true;
    }

    // 读取观测值
    std::string line;
    getline(rinexFileStream, line);

    if (rinexFileStream.eof()) {
        EndOfFile err("EOF encountered!");
        throw err;
    }

    if (debug) {
        std::cout << "current record line is:" << std::endl;
        std::cout << line << std::endl;
    }

    // 检查并解析历元行
    // 检查历元标记 ('>') 和随后的空格。
    if (line[0] != '>' || line[1] != ' ') {
        FFStreamError e("Bad epoch line: >" + line + "<");
        throw e;
    }

    int epochFlag = stoi(line.substr(31, 1));
    if (epochFlag < 0 || epochFlag > 6) {
        FFStreamError e("Invalid epoch flag: " + std::to_string(epochFlag));
        throw e;
    }

    CommonTime currEpoch = parseTime(line);
    if (debug) {
        std::cout << " currEpoch" << currEpoch << std::endl;
    }

    int numSats = stoi(line.substr(32, 3));

    if (debug) cout << numSats << endl;

    // 读取观测：SV ID 和数据
    SatTypeValueMap stvData;
    if (epochFlag == 0 || epochFlag == 1 || epochFlag == 6) {

        std::vector<SatID> satIndex(numSats);
        for (int isv = 0; isv < numSats; ++isv) {
            getline(rinexFileStream, line); // 修改了这里的变量名以匹配上下文

            if (debug) {
                cout << "parseRinexObs:" << line << endl;
            }

            if (rinexFileStream.eof()) {
                EndOfFile err("EOF encountered!");
                throw err;
            }

            // 获取 SV ID
            try {
                satIndex[isv] = SatID(line.substr(0, 3));
            } catch (std::exception &e) {
                FFStreamError ffse(e.what());
                throw ffse;
            }

            SatID sat = SatID(satIndex[isv]);

            // 如果卫星系统不是GPS("G")也不是北斗("C")，则跳过当前循环迭代。
            if (sat.system != "G" && sat.system != "C") {
                continue;
            }

            int size = rinexHeader.mapObsTypes.at(satIndex[isv].system).size();

            // 有些文件没有观测值，后面就没有输出，这里用空格来替换，否则解析错误
            size_t minSize = 3 + 16 * size;
            if (line.size() < minSize) {
                line += std::string(minSize - line.size(), ' ');
            }

            // 获取数据 (# entries in ObsType map of maps from header)
            TypeValueMap typeObs;
            TypeValueMap typeLLI;
            TypeValueMap typeSSI;
            for (int i = 0; i < size; ++i) {
                size_t pos = 3 + 16 * i;
                std::string str = line.substr(pos, 16);

                // ObsType
                std::string obsTypeStr = rinexHeader.mapObsTypes.at(sat.system)[i];

                // 观测值
                std::string tmpStr = str.substr(0, 14);

                double data = safeStod(tmpStr);

                // 载波相位
                if (obsTypeStr[0] == 'L') {
                    double wavelength = 0.0;

                    // 获取观测值频率，比如L1C，其频率为1
                    int n;
                    if (obsTypeStr[1] == 'A') {
                        n = 1;
                    } else {
                        n = stoi(obsTypeStr.substr(1, 1));
                    }

                    wavelength = getWavelength(sat.system, n);

                    if (wavelength == 0.0) continue;

                    if (debug) {
                        std::cout << obsTypeStr << " wavelength"
                                  << std::setprecision(12)
                                  << wavelength << std::endl;
                    }

                    // 将周期转换为米
                    data = data * wavelength;
                }

                // 观测值异常
                if (std::abs(data) == 0.0) {
                    continue;
                }

                typeObs[obsTypeStr] = data;
            }

            // 插入当前卫星的数据到 stvData
            stvData[satIndex[isv]] = typeObs;

        }
    }

    ObsData obsData;
    obsData.station = rinexHeader.station;
    obsData.epoch = currEpoch;
    obsData.satTypeValueData = stvData;

    return obsData;
}


CommonTime parseTime(const string &line) {

    // check if the spaces are in the right place - an easy
    // way to check if there's corruption in the file
    if ((line[1] != ' ') || (line[6] != ' ') || (line[9] != ' ') ||
        (line[12] != ' ') || (line[15] != ' ') || (line[18] != ' ') ||
        (line[29] != ' ') || (line[30] != ' ')) {
        FFStreamError e("Invalid time format");
        throw (e);
    }

    // if there's no time, just return a bad time
    if (line.substr(2, 27) == string(27, ' '))
        return BEGINNING_OF_TIME;

    int year, month, day, hour, min;
    double sec;

    year = stoi(line.substr(2, 4));
    month = stoi(line.substr(7, 2));
    day = stoi(line.substr(10, 2));
    hour = stoi(line.substr(13, 2));
    min = stoi(line.substr(16, 2));
    sec = safeStod(line.substr(19, 11));

    // Real Rinex has epochs 'yy mm dd hr 59 60.0' surprisingly often.
    double ds = 0;
    if (sec >= 60.0) {
        ds = sec;
        sec = 0.0;
    }

    CommonTime ctime;
    CivilTime cv = CivilTime(year, month, day, hour, min, sec);
    ctime = CivilTime2CommonTime(cv);

    if (ds != 0)
        ctime = ctime + ds;

    return ctime;


}  // end parseTime

void chooseObs(ObsData &obsData, std::map<std::string, std::set<std::string>> &sysTypes) {
    SatTypeValueMap filteredSatTypeValueData;

    // Iterate over all satellite entries in satTypeValueData
    for (const auto &satEntry: obsData.satTypeValueData) {
        const auto &satId = satEntry.first;
        const auto &typeValueMap = satEntry.second;

        // Check if the satellite's system is in sysTypes
        auto itSys = sysTypes.find(satId.system);
        if (itSys != sysTypes.end()) { // If the system is found in sysTypes
            const auto &allowedTypes = itSys->second;
            TypeValueMap filteredTypeValueMap;

            // Filter the observations based on the allowed types
            for (const auto &typeValueEntry: typeValueMap) {
                if (allowedTypes.find(typeValueEntry.first) != allowedTypes.end()) {
                    filteredTypeValueMap.insert(typeValueEntry);
                }
            }

            // Only add the satellite entry if there are any remaining observations
            if (!filteredTypeValueMap.empty()) {
                filteredSatTypeValueData[satId] = std::move(filteredTypeValueMap);
            }
        }
    }

    // Replace the original data with the filtered data
    obsData.satTypeValueData.swap(filteredSatTypeValueData);
}

// L1C => L1
// L2W => L2
void convertObsType(ObsData &obsData) {
    SatTypeValueMap stvData;
    for (auto sd: obsData.satTypeValueData) {
        TypeValueMap tvData;
        for (auto td: sd.second) {
            tvData[td.first.substr(0, 2)] = td.second;
        }
        stvData[sd.first] = tvData;
    }

    // 替代
    obsData.satTypeValueData = stvData;
};

static bool isCodeLikeObs(const std::string &obsType) {
    return obsType.size() >= 2 &&
           (obsType[0] == 'C' || obsType[0] == 'P') &&
           std::isdigit(static_cast<unsigned char>(obsType[1]));
}

static char codeFreqDigit(const std::string &obsType) {
    if (!isCodeLikeObs(obsType)) {
        return '\0';
    }
    return obsType[1];
}

static std::string tgdKeyForObs(const SatID &sat, const std::string &obsType) {
    const char freq = codeFreqDigit(obsType);
    if (freq == '\0') {
        return "";
    }

    if (sat.system == "G") {
        if (freq == '1' || freq == '2') {
            return "TGD";
        }
    } else if (sat.system == "C") {
        // BDS-2 / BDS-3 textbook mapping:
        // B1I (freq digit 2) -> TGD1, B2I/B2b (freq digit 7) -> TGD2.
        if (freq == '2') {
            return "TGD1";
        }
        if (freq == '7') {
            return "TGD2";
        }
    }

    return "";
}

static double tgdScaleForObs(const SatID &sat, const std::string &obsType) {
    const char freq = codeFreqDigit(obsType);
    if (freq == '\0') {
        return 0.0;
    }

    if (sat.system == "G") {
        if (freq == '1') {
            // Convert the broadcast DCB-like parameter into the textbook P1
            // delay using gamma2 = (f2/f1)^2.
            const double gamma2 = 1.0 / getGamma("G", "C1", "C2");
            return 1.0 / (1.0 - gamma2);
        }
        if (freq == '2') {
            const double gamma2 = 1.0 / getGamma("G", "C1", "C2");
            return gamma2 / (1.0 - gamma2);
        }
        if (freq == '5') {
            // CNAV / L5-specific ISC is not populated yet.
            return 0.0;
        }
    } else if (sat.system == "C") {
        // BDS B1I/B2I use TGD1/TGD2 directly in the single-frequency model.
        if (freq == '2' || freq == '7') {
            return 1.0;
        }
    }

    return 0.0;
}

static bool hasObsFreq(const TypeValueMap &typeValueMap, char freqDigit) {
    for (const auto &entry : typeValueMap) {
        if (codeFreqDigit(entry.first) == freqDigit && std::isfinite(entry.second)) {
            return true;
        }
    }
    return false;
}

static void appendBdsIfReport(const SatID &sat,
                              const TypeValueMap &typeValueMap,
                              const Xvt &xvt,
                              std::vector<ObsDelayCorrectionRecord> *report) {
    if (!report || sat.system != "C") {
        return;
    }

    if (!hasObsFreq(typeValueMap, '2') || !hasObsFreq(typeValueMap, '7')) {
        return;
    }

    auto tgd1It = xvt.typeTGDData.find("TGD1");
    auto tgd2It = xvt.typeTGDData.find("TGD2");
    if (tgd1It == xvt.typeTGDData.end() || tgd2It == xvt.typeTGDData.end()) {
        return;
    }

    const double f1 = getFreq("C", "C2");
    const double f2 = getFreq("C", "C7");
    const double ifSeconds = (f1 * f1 * tgd1It->second - f2 * f2 * tgd2It->second)
                             / (f1 * f1 - f2 * f2);

    ObsDelayCorrectionRecord rec;
    rec.sat = sat;
    rec.obsType = "IF(C2,C7)";
    rec.tgdSeconds = ifSeconds;
    rec.correctionSeconds = ifSeconds;
    rec.correctionMeters = C_MPS * ifSeconds;
    report->push_back(rec);
}

void correctTGD(ObsData &obsData,
                const std::map<SatID, Xvt> &satXvtTransTime,
                const ObsDelayBiasTable *biasTable,
                std::vector<ObsDelayCorrectionRecord> *report) {
    if (report) {
        report->clear();
    }

    for (auto &satEntry: obsData.satTypeValueData) {
        const SatID &sat = satEntry.first;
        auto xvtIt = satXvtTransTime.find(sat);
        if (xvtIt == satXvtTransTime.end()) {
            continue;
        }

        for (auto &tv: satEntry.second) {
            const std::string &obsType = tv.first;

            double correctionSeconds = 0.0;
            double tgdSeconds = 0.0;
            bool applied = false;

            if (biasTable != nullptr) {
                auto satBiasIt = biasTable->find(sat);
                if (satBiasIt != biasTable->end()) {
                    auto obsBiasIt = satBiasIt->second.find(obsType);
                    if (obsBiasIt != satBiasIt->second.end()) {
                        const SignalDelayTerm &term = obsBiasIt->second;
                        if (term.hasReferenceDelay) {
                            correctionSeconds = term.referenceDelaySeconds;
                            tgdSeconds = term.referenceDelaySeconds;
                            applied = true;
                        }
                        if (term.hasIsc) {
                            correctionSeconds -= term.iscSeconds;
                            applied = true;
                        }
                    }
                }
            }

            if (!applied) {
                std::string tgdKey = tgdKeyForObs(sat, obsType);
                if (tgdKey.empty()) {
                    continue;
                }

                auto tgdIt = xvtIt->second.typeTGDData.find(tgdKey);
                if (tgdIt == xvtIt->second.typeTGDData.end()) {
                    continue;
                }

                double scale = tgdScaleForObs(sat, obsType);
                if (scale == 0.0) {
                    continue;
                }

                tgdSeconds = tgdIt->second;
                correctionSeconds = scale * tgdSeconds;
            }

            double correctionMeters = C_MPS * correctionSeconds;
            tv.second -= correctionMeters;

            if (report) {
                ObsDelayCorrectionRecord rec;
                rec.sat = sat;
                rec.obsType = obsType;
                rec.tgdSeconds = tgdSeconds;
                rec.correctionSeconds = correctionSeconds;
                rec.correctionMeters = correctionMeters;
                report->push_back(rec);
            }
        }

        appendBdsIfReport(sat, satEntry.second, xvtIt->second, report);
    }
}

std::map<SatID, Xvt> computeSatPos(ObsData &obsData, RinexNavStore& navStore) {
    std::map<SatID, Xvt> satXvtData;
    SatIDSet satRejectedSet;
    CommonTime time = obsData.epoch;
    // Loop through all the satellites
    for (auto stv: obsData.satTypeValueData) {
        SatID sat(stv.first);
        Xvt xvt;
        // compute satellite ephemeris at transmitting time
        // Scalar to hold temporal value
        double obs(0.0);
        string codeType;
        if (sat.system == "G") {
            codeType = "C1";
        }
            // todo
            // 请增加bds或其他系统的观测值选择
        else {
            satRejectedSet.insert(sat);
            continue;
        }

        // code obs
        try {
            obs = stv.second.at(codeType);
            if(debug)
                cout << "sat:" << sat << "obs:" << codeType << "value:" << obs << endl;
        }
        catch (...) {
            satRejectedSet.insert(sat);
            continue;
        }

        // now, compute xvt
        try {
            xvt = computeAtTransmitTime(time, obs, sat, navStore);
        }
        catch (InvalidRequest &e) {
            satRejectedSet.insert(sat);
            continue;
        }
        satXvtData[sat] = xvt;
    }

    // remove bad sat;
    for (auto sat: satRejectedSet) {
        obsData.satTypeValueData.erase(sat);
    }

    return satXvtData;

};

Xvt computeAtTransmitTime(const CommonTime &tr,
                          const double &pr,
                          const SatID &sat,
                          RinexNavStore& navStore)
noexcept(false) {
    Xvt xvt;
    CommonTime tt;
    CommonTime transmit = tr;

    transmit -= pr / C_MPS;
    tt = transmit;

    // 这里也可以用while循环来替换这里的迭代次数
    for (int i = 0; i < 2; i++) {
        xvt = navStore.getXvt(sat, tt);
        tt = transmit;
        tt -= (xvt.clkbias + xvt.relcorr);
    }
    return xvt;
};

std::map<SatID, Xvt> earthRotation(Eigen::Vector3d &xyz,
                                   std::map<SatID, Xvt> &satXvtTransTime) {

    std::map<SatID, Xvt> satXvtRecTime;
    for(auto stv: satXvtTransTime) {
        SatID sat = stv.first;
        XYZ xyzSat(stv.second.x);
        double dt = (xyzSat - xyz).norm() / C_MPS;

        double wt(0.0);
        wt = OMEGA_EARTH * dt;

        // todo:
        // Eigen中Vector3d是不是支持坐标旋转？
        // 请查询并修改

        double xSat, ySat, zSat;
        xSat = stv.second.x[0];
        ySat = stv.second.x[1];
        zSat = stv.second.x[2];

        double xSatRot(0.0), ySatRot(0.0);
        xSatRot = +std::cos(wt) * xSat + std::sin(wt) * ySat;
        ySatRot = -std::sin(wt) * xSat + std::cos(wt) * ySat;

        XYZ xyzRecTime;
        xyzRecTime[0] = xSatRot;
        xyzRecTime[1] = ySatRot;
        xyzRecTime[2] = zSat; // z轴不变

        double vxSat, vySat, vzSat;
        vxSat = stv.second.v[0];
        vySat = stv.second.v[1];
        vzSat = stv.second.v[2];

        double vxSatRot(0.0), vySatRot(0.0);
        vxSatRot = +std::cos(wt) * vxSat + std::sin(wt) * vySat;
        vySatRot = -std::sin(wt) * vxSat + std::cos(wt) * vySat;

        XYZ velRecTime;
        velRecTime[0] = vxSatRot;
        velRecTime[1] = vySatRot;
        velRecTime[2] = vzSat; // 不变

        // 替换位置和速度，得到旋转后的卫星产品
        Xvt xvtRecTime = stv.second;
        xvtRecTime.x = xyzRecTime;
        xvtRecTime.v = velRecTime;

        satXvtRecTime[sat] = xvtRecTime;
    };

    return satXvtRecTime;
};

void computeElevAzim(Eigen::Vector3d& xyz,
                     std::map<SatID,Xvt> & satXvt,
                     SatValueMap& tempElevData,
                     SatValueMap& tempAzimData
)
{

    for(auto sx: satXvt)
    {
        SatID sat = sx.first;

        XYZ satXYZ = sx.second.x;

        // elevation
        double elev(0.0);
        double azim(0.0);
        elev = elevation(xyz, satXYZ);
        azim = azimuth(xyz, satXYZ);

        tempElevData[sat] = elev;
        tempAzimData[sat] = azim;
    }
};

double wavelengthOfMW(string sys, string L1Type, string L2Type) {
    double f1 = getFreq(sys, L1Type);
    double f2 = getFreq(sys, L2Type);
    double wavelength = C_MPS / (f1 - f2);
    return wavelength;
};

double varOfMW(string, string L1Type, string L2Type) {
    double var = sqrt(2.0) / 2 * 0.3;
    return var;
};

void detectCSMW(ObsData &obsData,
                std::map<Variable, int> &csFlagData,
                SatEpochValueMap &satEpochMWData,
                SatEpochValueMap &satEpochMeanMWData,
                SatEpochValueMap &satEpochCSFlagData) {
    //==================
    // 初始化常数和static变量
    //==================
    double deltaTMax(120.0);
    double minCycles(2.0);

    // A structure used to store filter data for a SV.
    struct MWData {
        // Default constructor initializing the data in the structure
        MWData()
                : formerEpoch(BEGINNING_OF_TIME), windowSize(0), meanMW(0.0), varMW(0.0) {};

        CommonTime formerEpoch; ///< The previous epoch time stamp.
        int windowSize;         ///< Size of current window, in samples.
        double meanMW;          ///< Accumulated mean value of combination.
        double varMW;           ///< Accumulated std value of combination.
    };

    // 这个数据在下次调用时需要用到，所以定位为static变量
    static std::map<SatID, MWData> satMWData;

    //==========================
    // 逐个卫星做周跳探测
    //==========================
    // Loop through all the satellites
    CommonTime currentEpoch = obsData.epoch;
    SatIDSet badSatSet;
    for (auto stv: obsData.satTypeValueData) {
        SatID sat = (stv).first;
        string L1Type, L2Type, C1Type, C2Type;
        if (sat.system == "G") {
            L1Type = "L1";
            L2Type = "L2";
            C1Type = "C1";
            C2Type = "C2";
        } else // 请增加bds的处理
        {
            badSatSet.insert(sat);
        }

        // wavelengthMW of MW-combination, see LinearCombination
        double wavelengthMW = wavelengthOfMW(sat.system, L1Type, L2Type);
        double varianceMW = varOfMW(sat.system, L1Type, L2Type);

        double f1 = getFreq(sat.system, L1Type);
        double f2 = getFreq(sat.system, L2Type);

        if (debug) {
            cout << "f1:" << f1 << "f2:" << f2 << endl;
        }

        double L1Value, L2Value, C1Value, C2Value, mwValue;

        try {
            L1Value = stv.second.at(L1Type);
            L2Value = stv.second.at(L2Type);
            C1Value = stv.second.at(C1Type);
            C2Value = stv.second.at(C2Type);

            mwValue
                    = (f1 * L1Value - f2 * L2Value) / (f1 - f2)
                      - (f1 * C1Value + f2 * C2Value) / (f1 + f2);
        } catch (std::out_of_range) {
            // 无法构成mw，这个卫星观测值周跳无法探测，删除这个卫星
            badSatSet.insert(sat);
            continue; // 继续处理下一个卫星
        }

        satEpochMWData[sat][currentEpoch] = mwValue;

        if (debugCSMW) {
            cout << "L1Value:" << L1Value << endl;
            cout << "L2Value:" << L2Value << endl;
            cout << "C1Value:" << C1Value << endl;
            cout << "C2Value:" << C2Value << endl;
            cout << "mwValue:" << mwValue << endl;
            cout << "wavelength:" << C_MPS / (f1 - f2) << endl;
        }

        //-------------------
        double currentDeltaT(0.0);
        double currentBias(0.0);
        int csFlag(0.0);

        currentDeltaT = (currentEpoch - satMWData[sat].formerEpoch);
        satMWData[sat].formerEpoch = currentEpoch;
        if (debugCSMW) {
            cout << "currentDeltaT:" << currentDeltaT << endl;
        }
        // Difference between current value of MW and average value
        currentBias = std::abs(mwValue - satMWData[sat].meanMW);
        if (debugCSMW) {
            cout << "currentBias:" << currentBias << endl;
        }

        // Increment window size
        satMWData[sat].windowSize++;

        /**
         * cycle-slip condition
         * 1. if data interrupt for a given time gap, then cyce slip should be set
         * 2. if current bias is greater than 1 cycle and greater than 4 sigma of mean mw.
         */
        double sigLimit = 4 * std::sqrt(satMWData[sat].varMW);

        if (debugCSMW) {
            cout << "deltaTMax:" << deltaTMax << endl;
            cout << "wavelengthMW:" << wavelengthMW << endl;
            cout << "sigLimit:" << sigLimit << endl;
            cout << "minCycles:" << minCycles * wavelengthMW << endl;
        }

        // 波长有可能为负值
        if (currentDeltaT > deltaTMax ||
            currentBias > std::abs(minCycles * wavelengthMW) ||
            currentBias > sigLimit) {

            // reset the filter window size/meanMW/InitialVarofMW
            satMWData[sat].meanMW = mwValue;
            satMWData[sat].varMW = varianceMW;
            satMWData[sat].windowSize = 1;

            if (debugCSMW) {
                cout << "* CS happened!" << endl;
            }
            csFlag = 1.0;
        } else {
            // MW bias from the mean value
            double mwBias(mwValue - satMWData[sat].meanMW);
            double size(static_cast<double>(satMWData[sat].windowSize));

            // Compute average
            satMWData[sat].meanMW += mwBias / size;

            // Compute variance
            // Var(i) = Var(i-1) + [ ( mw(i) - meanMW)^2/(i)- 1*Var(i-1) ]/(i);
            satMWData[sat].varMW += (mwBias * mwBias - satMWData[sat].varMW) / size;
        }

        // for print
        satEpochMeanMWData[sat][currentEpoch] = satMWData[sat].meanMW;

        // 放大到mw数值，以方便绘图
        satEpochCSFlagData[sat][currentEpoch] = csFlag * mwValue;

        // 将周跳探测标志存到模糊度变量中
        Variable amb1(obsData.station, sat, static_cast<Parameter>(Parameter::ambiguity), ObsID(sat.system, L1Type));
        Variable amb2(obsData.station, sat, static_cast<Parameter>(Parameter::ambiguity), ObsID(sat.system, L2Type));

        csFlagData[amb1] = csFlag;
        csFlagData[amb2] = csFlag;
    }

    // 删除坏卫星
    for (auto sat: badSatSet)
        obsData.satTypeValueData.erase(sat);

};

void differenceStation(EquSys& equSysRover,
                       EquSys& equSysBase,
                       EquSys& equSysSD)
{
    // 逐个观测值类型取出类型
    std::map<EquID, EquData> obsEquDataDiff;
    VariableSet varSetDiff;
    for(auto& oe: equSysRover.obsEquData)
    {
        // 在参考站中查找当前观测值，如果没有找到就跳过;
        // 需要注意测站名是不同的， 只需要查找卫星号和观测值
        if(equSysBase.obsEquData.find(oe.first)==equSysBase.obsEquData.end())
        {
            continue;
        }

        // 取出参考站的线性化残差观测值
        EquData oeBase = equSysBase.obsEquData.at(oe.first);

        // 在参考站中查找类型的观测值，找到了就计算站间差分观测值
        double diffPrefit;
        diffPrefit = oe.second.prefit - oeBase.prefit;

        if(debug)
        {
            cout << ">>>>>>>>> differenceStation" << endl;
            cout << oe.first
            << "rover: " << oe.second.prefit
            << "base:  " << oeBase.prefit
            << "diffPrefit:" << diffPrefit
            << endl;
        }

        // 把当前卫星的tvDiff数据插入到gDataDiff;
        obsEquDataDiff[oe.first].prefit = diffPrefit;

        //------------------------------------------------------------
        // 因为对于短基线来说，可以不用估计电离层和对流程，
        // 这里为了简单起见，直接将电离层和对流层参数从站间差分未知参数表中删除
        // todo:
        // 更优雅的处理方式是在参数估计时，对电离层和电离层进行约束，
        // 并根据基线长度对约束的方差进行动态调整。
        // 比如：
        // iono = 0, sigmaIono = 0.001*0.001*baseline
        // trop = 0, sigmaTrop = 0.0001*.0001*baseline
        // 通过增加电离层、对流层约束方程，实现通用rtk定位模型
        //------------------------------------------------------------
        // 未知参数与流动站的参数是相同的。
        std::map<Variable, double> vcDataTemp;
        for(auto vc: oe.second.varCoeffData)
        {
            if(vc.first.getParaType() != Parameter::iono)
            {
                vcDataTemp[vc.first] = vc.second;
                varSetDiff.insert(vc.first);
            }
        }
        obsEquDataDiff[oe.first].varCoeffData = vcDataTemp;

        // 权函数
        double weightRover = oe.second.weight;
        double weightBase = oeBase.weight;
        double varDiff = 1.0/weightRover + 1.0/weightBase;
        obsEquDataDiff[oe.first].weight = 1.0/varDiff;
    }

    equSysSD.obsEquData = obsEquDataDiff;
    equSysSD.varSet = varSetDiff;

};



SatID findDatumSat(bool& firstEpoch,
                   SatValueMap& satElevData) {
    // 确定基准卫星，必须是上一个历元已经固定的卫星才能选作基准
    SatID datumSat;
    auto maxIt = max_element(satElevData.begin(),satElevData.end(),
                              [](const auto& a, const auto& b){ return a.second < b.second; } );
    datumSat = maxIt->first;
    return datumSat;
}

void differenceSat( SatID& datumSat,
                    EquSys& equSysSD,
                    EquSys& equSysDD ) {
    //----------------------------------------------------
    // 根据基准卫星，选择每个观测类型的观测值，并将其他的与基准卫星对应观测值求差
    // warning:
    // 因为星间单差需要消除接收机钟差和接收机端硬件延迟，
    // 因此必须为每个类型独立构建星间单差观测方程，而不能混合在一起；
    // 因此，基准观测值的方程数据，应该存在以观测类型为key键值的map中，
    // 由于观测类型我们采用了C1，C2，L1，L2作为名字，
    // 当采用GPS+BDS时，两个系统均存在L2，无法有效区分，
    // 因此，这里需要创建一个独立的数据结构ObsID来管理观测类型ID，
    // 其由两个成员构成，一个是obsType；一个是卫星系统system
    // 另一个简单的处理：
    // string obsStr = obsType + system;
    //----------------------------------------------------

    std::map<ObsID, EquData> datumEquData;
    std::map<EquID, EquData> otherEquData;

    std::map<EquID, EquData> equData;
    equData = equSysSD.obsEquData;
    if (debug)
    {
        cout << "differenceSat:" << endl;
        cout << "datumSat:" << datumSat << endl;
    }

    for(auto ed: equData)
    {
        if(ed.first.sat == datumSat)
        {
            ObsID obsID(ed.first.sat.system, ed.first.obsType);
            if(debug)
                cout << "datum obsid:" << obsID << endl;

            datumEquData[obsID] = ed.second;
        }
        else
        {
            otherEquData[ed.first] = ed.second;
        }
    }

    // dd
    // 先验残差求差；
    // dx，dy，dz的系数求差；
    // 接收机钟差进一步差分掉了；
    // 模糊度除了基准卫星，其他卫星变成双差模式，系数不变
    std::map<EquID, EquData> equDataDD;
    VariableSet varSetDD;
    for(auto ed: otherEquData)
    {
        // 先验残差
        // 需要在基准ObsID里找EquData，来构成星间差分，
        // 如果找不到就剔除这个卫星；
        // 因此需要捕获异常，来处理找不到的情况；
        cout << "differenceSat:" << "sat:" << ed.first.sat << endl;
        double prefitDatum;
        ObsID currentObsID = ObsID(ed.first.sat.system, ed.first.obsType);
        try {
            prefitDatum= datumEquData.at(currentObsID).prefit;

            // dd prefit
            double prefitDD = ed.second.prefit - prefitDatum;
            equDataDD[ed.first].prefit = prefitDD;

            // 系数与未知参数
            VariableDataMap vcDatum = datumEquData.at(currentObsID).varCoeffData;

            // 接收机钟差消除了，只保留了坐标和模糊度参数
            for(auto vc: ed.second.varCoeffData)
            {
                if( vc.first.getParaType()==Parameter::dX ||
                    vc.first.getParaType()==Parameter::dY ||
                    vc.first.getParaType()==Parameter::dZ )
                {
                    double coeffDiff;
                    coeffDiff = vc.second - vcDatum.at(vc.first);
                    equDataDD[ed.first].varCoeffData[vc.first] = coeffDiff;
                    varSetDD.insert(vc.first);
                }
                else if(vc.first.getParaType()==Parameter::ambiguity)
                {
                    equDataDD[ed.first].varCoeffData[vc.first] = vc.second;
                    varSetDD.insert(vc.first);
                }
            }

            // 双差的方差近似等于单差观测值的方差的和；
            double weightCurrent = ed.second.weight;
            double weightDatum = datumEquData.at(currentObsID).weight;
            double varDiff = 1.0/weightCurrent + 1.0/weightDatum;
            equDataDD[ed.first].weight = 1.0/varDiff;

            //
            // todo
            // 构建完整的方差协方差阵，并比较定位结果的不同
        }
        catch(...)
        {
            continue;
        }
    }
    equSysDD.obsEquData = equDataDD;
    equSysDD.varSet = varSetDD;
};

// 对于Kalman滤波来说，需要对流动站和参考站周跳进行周跳标识符的合并，
// 只要流动站和参考站对应频率模糊度有一个发生了周跳就需要对周跳进行合并；
// 可以通过重载函数来实现对现有函数功能的复用。
void differenceStation(EquSys& equSysRover, VariableDataMap& csFlagRover,
                       EquSys& equSysBase, VariableDataMap& csFlagBase,
                       EquSys& equSysSD, VariableDataMap& csFlagSD)
{
    differenceStation(equSysRover,
                      equSysBase,
                      equSysSD);

    string roverStation = equSysRover.station;

    // todo:
    // 去掉流动站和参考站的站名，否则无法查找并匹配周跳
    VariableDataMap tempFlagRover;
    for(auto vd: csFlagRover)
    {
        Variable tempVar = vd.first;
        tempVar.station = std::string("");
        tempFlagRover[tempVar] = vd.second;
    }

    // 去掉基准站名字
    VariableDataMap tempFlagBase;
    for(auto vd: csFlagBase)
    {
        Variable tempVar = vd.first;
        tempVar.station = std::string("");
        tempFlagBase[tempVar] = vd.second;
    }

    // 现在，从基准站中寻找流动站模糊度，如果找到了，就把周跳标志合并
    // 如果没找到，就跳过，说明无法形成站间差分观测
    for(auto vd: tempFlagRover) {
        double flagRover = vd.second;
        if (tempFlagBase.find(vd.first) != tempFlagBase.end())
        {
            double flagBase = tempFlagBase.at(vd.first);
            double flagSD(0.0);
            // 基准站或者流动站一个发生周跳，就标志周跳
            if(flagRover|| flagBase)
            {
                flagSD = 1.0;
            }
            Variable varSD = vd.first; // 得到流动站模糊度变量
            varSD.station = roverStation; // 把流动站名站再次赋值进来
            csFlagSD[varSD] = flagSD;
        }
    }

    // 单差周跳
    cout << "differenceStation:" << "csFlagRover:" << endl;
    for(auto cd:csFlagRover)
    {
        cout << "cs:" << cd.first << " flag:" << cd.second;
    }

    cout << "differenceStation:" << "csFlagBase:" << endl;
    for(auto cd:csFlagBase)
    {
        cout << "cs:" << cd.first << " flag:" << cd.second;
    }

    cout << "differenceStation:" << "csFlagSD:" << endl;
    for(auto cd:csFlagSD)
    {
        cout << "cs:" << cd.first << " flag:" << cd.second;
    }

};

void differenceSat( SatID& datumSat,
                    EquSys& equSysSD, VariableDataMap& csFlagSD,
                    EquSys& equSysDD, VariableDataMap& csFlagDD ) {
    //----------------------------------------------------
    // 根据基准卫星，选择每个观测类型的观测值，并将其他的与基准卫星对应观测值求差
    // warning:
    // 因为星间单差需要消除接收机钟差和接收机端硬件延迟，
    // 因此必须为每个类型独立构建星间单差观测方程，而不能混合在一起；
    // 因此，基准观测值的方程数据，应该存在以观测类型为key键值的map中，
    // 由于观测类型我们采用了C1，C2，L1，L2作为名字，
    // 当采用GPS+BDS时，两个系统均存在L2，无法有效区分，
    // 因此，这里需要创建一个独立的数据结构ObsID来管理观测类型ID，
    // 其由两个成员构成，一个是obsType；一个是卫星系统system
    // 另一个简单的处理：
    // string obsStr = obsType + system;
    //----------------------------------------------------

    std::map<ObsID, EquData> datumEquData;
    std::map<EquID, EquData> otherEquData;

    std::map<EquID, EquData> equData;
    equData = equSysSD.obsEquData;

    for(auto ed: equData)
    {
        if(ed.first.sat == datumSat)
        {
            ObsID obsID(ed.first.sat.system, ed.first.obsType);
            if(debug)
                cout << "datum obsid:" << obsID << endl;

            datumEquData[obsID] = ed.second;
        }
        else
        {
            otherEquData[ed.first] = ed.second;
        }
    }

    // dd
    // 先验残差求差；
    // dx，dy，dz的系数求差；
    // 接收机钟差进一步差分掉了；
    // 模糊度除了基准卫星，其他卫星变成双差模式，系数不变
    std::map<EquID, EquData> equDataDD;
    VariableSet varSetDD;
    for(auto ed: otherEquData)
    {
        // 先验残差
        // 需要在基准ObsID里找EquData，来构成星间差分，
        // 如果找不到就剔除这个卫星；
        // 因此需要捕获异常，来处理找不到的情况；
        cout << "differenceSat:" << "sat:" << ed.first.sat << endl;
        double prefitDatum;
        ObsID currentObsID = ObsID(ed.first.sat.system, ed.first.obsType);
        try {
            prefitDatum= datumEquData.at(currentObsID).prefit;

            // dd prefit
            double prefitDD = ed.second.prefit - prefitDatum;
            equDataDD[ed.first].prefit = prefitDD;

            // 系数与未知参数
            VariableDataMap vcDatum = datumEquData.at(currentObsID).varCoeffData;

            // 接收机钟差消除了，只保留了坐标和模糊度参数
            for(auto vc: ed.second.varCoeffData)
            {
                if( vc.first.getParaType()==Parameter::dX ||
                    vc.first.getParaType()==Parameter::dY ||
                    vc.first.getParaType()==Parameter::dZ )
                {
                    double coeffDiff;
                    coeffDiff = vc.second - vcDatum.at(vc.first);
                    equDataDD[ed.first].varCoeffData[vc.first] = coeffDiff;
                    varSetDD.insert(vc.first);
                }
                else if(vc.first.getParaType()==Parameter::ambiguity)
                {
                    std::pair<Variable, double> ambData;
                    for(auto vc2: vcDatum)
                    {
                        if(vc2.first.getParaType() == Parameter::ambiguity)
                        {
                            ambData.first = vc2.first;
                            ambData.second = vc2.second;
                        }
                    }
                    // 把基准模糊度插入到方程中，也就是估计基准模糊度，而不是合并成双差模糊度
                    // warning: 是负号
                    equDataDD[ed.first].varCoeffData[ambData.first] = -ambData.second;
                    equDataDD[ed.first].varCoeffData[vc.first] = vc.second;
                    varSetDD.insert(ambData.first);
                    varSetDD.insert(vc.first);
                }
            }

            // 双差的方差近似等于单差观测值的方差的和；
            double weightCurrent = ed.second.weight;
            double weightDatum = datumEquData.at(currentObsID).weight;
            double varDiff = 1.0/weightCurrent + 1.0/weightDatum;
            equDataDD[ed.first].weight = 1.0/varDiff;

            //
            // todo
            // 构建完整的方差协方差阵，并比较定位结果的不同
        }
        catch(...)
        {
            continue;
        }
    }
    equSysDD.obsEquData = equDataDD;
    equSysDD.varSet = varSetDD;

    // 直接把站间单差模糊度标志给双差即可，因为估计的模糊度仍然为站间单差模糊度
    csFlagDD = csFlagSD;

};

void fixSolution(VectorXd& stateVec,
                 MatrixXd& covMatrix,
                 VariableSet& varSet,
                 double& ratio,
                 Vector3d& dxyzFixed,
                 VariableDataMap& fixedAmbData)
{
    // 按照卫星把模糊度进行分类；
    VariableSet ambVarSet;
    for (auto var: varSet) {
        if(var.getParaType()==Parameter::ambiguity)
        {
            ambVarSet.insert(var);
        }
    };

    int numAmb = ambVarSet.size();
    int numXYZT = varSet.size() - numAmb;

    // 取出来星间差分模糊度ambVarSetSD的估值和方差，利用lambda方法固定
    VectorXd ambSol;
    MatrixXd ambCov;
    ambSol = stateVec.tail(numAmb);
    ambCov = covMatrix.block(numXYZT, numXYZT, numAmb, numAmb);

    if (debug) {
        cout << "ambVarSet" << endl;
        for (auto var: ambVarSet) {
            cout << var << " ";
        }
        cout << endl;

        cout << "ambSol" << endl;
        cout << ambSol.transpose() << endl;

        cout << "ambCov" << endl;
        cout << ambCov << endl;
    }

    // 如果ratio值大于3，则表明模糊度可以固定。
    ARLambda arLambda;
    VectorXd ambSolFixed = arLambda.resolve(ambSol, ambCov);

    int iamb=0;
    fixedAmbData.clear();
    for(auto var:ambVarSet)
    {
        fixedAmbData[var] = ambSolFixed(iamb);
        iamb++;
    }

    ratio = arLambda.squaredRatio;

    if (debug) {
        cout << "ratio:" << endl;
        cout << ratio << endl;
    }

    VectorXd xVecFixed;
    VectorXd xVec = VectorXd::Zero(numXYZT);

    // x/y/z
    xVec = stateVec.head(numXYZT);

    MatrixXd Qxx = covMatrix.block(0, 0, numXYZT, numXYZT);
    MatrixXd Qxb = covMatrix.block(0, numXYZT, numXYZT, numAmb);
    MatrixXd Qbb = covMatrix.block(numXYZT, numXYZT, numAmb, numAmb);

    if (debug) {
        cout << fixed << setprecision(5) << endl;

        cout << "covMatrix" << endl;
        cout << covMatrix << endl;

        cout << "Qxx" << endl;
        cout << Qxx << endl;

        cout << "Qxb" << endl;
        cout << Qxb << endl;

        cout << "Qbb" << endl;
        cout << Qbb << endl;
    }

    xVecFixed = xVec - Qxb * Qbb.inverse() * (ambSol - ambSolFixed);

    // return fixed solutions
    dxyzFixed = xVecFixed;

};



void ambiguityDatum(bool& firstEpoch,
                    SatID& datumSat,
                    VariableDataMap& fixedAmbData,
                    EquSys& equSysDD){

    // 对于第一个历元，直接将基准卫星模糊度固定为零即可。
    if(firstEpoch)
    {
        for(auto var:equSysDD.varSet)
        {
            if(var.getSat() == datumSat)
            {
                EquID equIDDatum;
                equIDDatum.sat = datumSat;
                equIDDatum.obsType = var.getParaType().toString()
                                     + var.getObsID().toString();

                EquData equDataDatum;
                equDataDatum.prefit = 0.0;
                equDataDatum.varCoeffData[var] = 1.0;
                equDataDatum.weight = 1.0E+8;

                // 将模糊度基准观测方程加入到观测系统中
                equSysDD.obsEquData[equIDDatum] = equDataDatum;
            }
        }
    }
    else
    {
        for(auto vd: fixedAmbData)
        {
            // 生成模糊度基准的观测方程
            if(vd.first.getSat()==datumSat)
            {
                EquID equIDDatum;
                equIDDatum.sat = datumSat;
                equIDDatum.obsType = vd.first.getParaType().toString()
                                     + vd.first.getObsID().toString();

                EquData equDataDatum;
                equDataDatum.prefit = vd.second;
                equDataDatum.varCoeffData[vd.first] = 1.0;
                equDataDatum.weight = 1.0E+8;

                // 将模糊度基准观测方程加入到观测系统中
                equSysDD.obsEquData[equIDDatum] = equDataDatum;
            }
        }
    }
};

// print solution to files
void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover,
                   Eigen::Vector3d& xyzRTKFloat,
                   double& ratio,
                   Eigen::Vector3d& xyzRTKFixed)
{
    YDSTime ydsTime = CommonTime2YDSTime(ctTime);
    solStream
            << ydsTime
            << fixed << setprecision(3)
            << "spp: " << xyzRover.transpose()
            << " float-rtk: " << xyzRTKFloat.transpose()
            << " ratio:" << ratio
            << " fixed-rtk:"<< xyzRTKFixed.transpose()
            << endl;
};

// print solution to files
void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover,
                   Eigen::Vector3d& xyzRTKFloat)
{
    YDSTime ydsTime = CommonTime2YDSTime(ctTime);
    solStream
    << ydsTime
    << fixed << setprecision(3)
    << "spp: " << xyzRover.transpose()
    << " rtk: " << xyzRTKFloat.transpose() << endl;
};

void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover)
{
    YDSTime ydsTime = CommonTime2YDSTime(ctTime);
    solStream
    << ydsTime
    << " "
    << fixed << setprecision(3)
    << xyzRover.transpose() << endl;
};
