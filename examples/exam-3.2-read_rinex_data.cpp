/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: shoujian zhang，shjzhang@sgg.whu.edu.cn， 2024-10-10
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernández-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */


#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <set>
#include "GnssStruct.h"
#include "TimeConvert.h"
#include "GnssFunc.h"
#include "RinexNavStore.hpp"
#include "RinexObsReader.h"
#include "SPPIFCode.h"

#define debug 1

using namespace std;

int main() {

    //--------------------
    // 打开文件流
    //--------------------

    // Replace with your actual RINEX file path
    string dirPath = "D:\\documents\\Source\\gnssLab-2.2\\data\\Zero-baseline\\";
//    string dirPath = "D:\\documents\\Source\\gnssLab-2.1\\data\\";

    // rover obs file name
    std::string roverFile = dirPath + "oem719-202203031500-1.obs";
    //std::string roverFile = dirPath + "ABMF00GLP_R_20210010000_01D_30S_MO.rnx";

    cout << roverFile << endl;

    // nav file name, download from IGS ftp site:ftp://gssc.esa.int/gnss/data/daily/YYYY/brdc
    std::string navFile = dirPath + "BRDC00IGS_R_20220620000_01D_MN.rnx";
    // std::string navFile = dirPath + "ABMF00GLP_R_20210010000_01D_MN.rnx";

    // 2022 03 03 06 48 37.0000000  0 45
    CivilTime stopCivilTime = CivilTime(2022, 03, 03, 06, 48 , 37);
    //   CivilTime stopCivilTime = CivilTime(2021, 01, 01, 10, 00 , 00);

    CommonTime stopEpoch = CivilTime2CommonTime(stopCivilTime);

    std::fstream roverObsStream(roverFile);
    if (!roverObsStream) {
        cerr << "rover file open error!" << strerror(errno) << endl;
        exit(-1);
    }

    // read nav file data before rtk
    RinexNavStore navStore;
    navStore.loadFile(navFile);

    cout << "after NavStore" << endl;

    std::map<string, std::set<string>> selectedTypes;
    selectedTypes["G"].insert("C1C");
    selectedTypes["G"].insert("C2W");
    selectedTypes["G"].insert("L1C");
    selectedTypes["G"].insert("L2W");


    std::map<string, std::pair<string, string>> ifCodeTypes;
    ifCodeTypes["G"].first = "C1";
    ifCodeTypes["G"].second = "C2";

    //-------------------
    // 定义数据处理的对象
    //-------------------
    //>>> classes for rover
    RinexObsReader readObsRover;
    readObsRover.setFileStream(&roverObsStream);
    readObsRover.setSelectedTypes(selectedTypes);


    std::string solFile = roverFile + ".spp.out";
    if(debug)
        cout << solFile << endl;

    while (true) {

        // solve spp for rover
        ObsData roverData;

        try {
            roverData = readObsRover.parseRinexObs();
            cout << "roverData:" << roverData << endl;
        }
        catch (EndOfFile &e) { break; }

        CommonTime epoch = roverData.epoch;


        // 调试代码时，设置一个stopEpoch，有助于快速得到结果
        if (roverData.epoch > stopEpoch)
            break;

    }

    roverObsStream.close();

}

