/**
 * Copyright:
 * ---------
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author:
 * ---------
 * Shoujian Zhang，shjzhang@sgg.whu.edu.cn， 2024-10-10
 *
 * References:
 * ---------
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernández-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */


#include "StringUtils.h"
#include "RinexObsReader.h"
#include "TimeConvert.h"

#define debug 0
void RinexObsReader::parseRinexHeader() {

    double version;
    XYZ antennaPosition;
    string satSys;
    std::map<string, std::vector<string>> mapObsTypes;
    while (true) {
        string line;
        getline(*pFileStream, line);

        cout << "parseRinexHeader:" << line << endl;

        string label;
        label = strip(line.substr(60, 20));

        if(debug)
            cout << "label:" << label << endl;

        if (label == "END OF HEADER") {
            break;
        }
        else if (label == "MARKER NAME") {
            string markerName = strip(line.substr(0, 60));
            std::replace(markerName.begin(), markerName.end(), ' ', '_');
            rinexHeader.station = strip(markerName);
            cout << "markerName:" << markerName << endl;
        }
        else if (label == "RINEX VERSION / TYPE") {
            version = safeStod(line.substr(0, 20));
            if (version < 3.04 || version > 3.05) {
                cerr << "only support rinex 3.04/3.05 version!" << endl;
                exit(-1);
            }
            rinexHeader.version = version;
        }
        else if (label == "APPROX POSITION XYZ") {
            antennaPosition[0] = safeStod(line.substr(0, 14));
            antennaPosition[1] = safeStod(line.substr(14, 14));
            antennaPosition[2] = safeStod(line.substr(28, 14));
            rinexHeader.antennaPosition = antennaPosition;
        }
        else if (label == "SYS / # / OBS TYPES") {
            string sysStr;
            sysStr = line.substr(0, 1);
            strip(sysStr);

            int numObs;

            if (sysStr != "") {
                numObs = safeStoi(line.substr(3, 3));
                satSys = sysStr;
            }

            if(debug)
                cout << "numObs:" << numObs << endl;

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

ObsData RinexObsReader::parseRinexObs() {

    if (!isHeaderRead) {
        parseRinexHeader();
        isHeaderRead = true;
    }

    // 读取观测值
    std::string line;
    getline(*pFileStream, line);

    if ((*pFileStream).eof()) {
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

    int epochFlag = safeStoi(line.substr(31, 1));
    if (epochFlag < 0 || epochFlag > 6) {
        FFStreamError e("Invalid epoch flag: " + std::to_string(epochFlag));
        throw e;
    }

    CommonTime currEpoch = parseTime(line);
    if (debug) {
        std::cout << " currEpoch" << currEpoch << std::endl;
    }

    int numSats = safeStoi(line.substr(32, 3));

    if (debug) cout << numSats << endl;

    // 读取观测：SV ID 和数据
    SatTypeValueMap stvData;
    if (epochFlag == 0 || epochFlag == 1 || epochFlag == 6) {

        std::vector<SatID> satIndex(numSats);
        for (int isv = 0; isv < numSats; ++isv) {
            getline(*pFileStream, line); // 修改了这里的变量名以匹配上下文

            if (debug) {
                cout << "parseRinexObs:" << line << endl;
            }

            if ((*pFileStream).eof()) {
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
                        n = safeStoi(obsTypeStr.substr(1, 1));
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
    obsData.antennaPosition = rinexHeader.antennaPosition;

    // choose observations you selected
    chooseObs(obsData);

    return obsData;
}


CommonTime RinexObsReader::parseTime(const string &line) {

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

    year = safeStoi(line.substr(2, 4));
    month = safeStoi(line.substr(7, 2));
    day = safeStoi(line.substr(10, 2));
    hour = safeStoi(line.substr(13, 2));
    min = safeStoi(line.substr(16, 2));
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

// 这里需要修改，保证每个卫星需要的字段必须存在观测值，
// 否则删除，因为涉及后续观测值基准选择问题；
void RinexObsReader::chooseObs(ObsData &obsData) {
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
