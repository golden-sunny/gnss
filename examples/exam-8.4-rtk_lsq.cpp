/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author:
 *  Shoujian Zhang，shjzhang@sgg.whu.edu.cn， 2024-10-10
 *
 * References:
 *  1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernández-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 *  2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
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
#include "SPPUCCodePhase.h"

using namespace std;

#define debug 1

int main() {

    //--------------------
    // 打开文件流
    //--------------------

    // Replace with your actual RINEX file path
    string dirPath = "D:\\documents\\Source\\gnssLab-2.2\\data\\Zero-baseline\\";

    // rover obs file name
    std::string roverFile = dirPath + "oem719-202203031500-1.obs";
    std::string baseFile = dirPath + "oem719-202203031500-2.obs";
    cout << roverFile << endl;

    // nav file name, download from IGS ftp site:ftp://gssc.esa.int/gnss/data/daily/YYYY/brdc
    std::string navFile = dirPath + "BRDC00IGS_R_20220620000_01D_MN.rnx";

    std::fstream roverObsStream(roverFile);
    if (!roverObsStream) {
        cerr << "rover file open error!" << strerror(errno) << endl;
        exit(-1);
    }

    std::fstream baseObsStream(baseFile);
    if (!baseObsStream) {
        cerr << "base file open error!" << strerror(errno) << endl;
        exit(-1);
    }

    // read nav file data before rtk
    RinexNavStore navStore;
    navStore.loadFile(navFile);

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

    std::map<string, std::pair<string, string>> dualCodeTypes;
    dualCodeTypes["G"].first = "C1";
    dualCodeTypes["G"].second = "C2";

    SPPUCCodePhase sppUCCodePhaseRover;
    sppUCCodePhaseRover.setRinexNavStore(&navStore);
    sppUCCodePhaseRover.setDualCodeTypes(dualCodeTypes);

    //>>> classes for base
    RinexObsReader readObsBase;
    readObsBase.setFileStream(&baseObsStream);
    readObsBase.setSelectedTypes(selectedTypes);
    
    SPPUCCodePhase sppUCCodePhaseBase;
    sppUCCodePhaseBase.setRinexNavStore(&navStore);
    sppUCCodePhaseBase.setStationAsBase();
    sppUCCodePhaseBase.setDualCodeTypes(dualCodeTypes);

    //>>> classes for rtk;    
    SolverLSQ solverRTK;

    CivilTime stopCivilTime = CivilTime(2022, 03, 03, 06, 49, 00);
    CommonTime stopEpoch = CivilTime2CommonTime(stopCivilTime);

    std::string solFile = roverFile + ".rtk.lsq.out";
    std::fstream solStream(solFile, ios::out);
    if (!solStream) {
        cerr << "solution file open error!" << strerror(errno) << endl;
        exit(-1);
    }

    while (true) {

        // solve spp for rover
        ObsData roverData;

        try {
            roverData = readObsRover.parseRinexObs();
        }
        catch (EndOfFile &e) { break; }

        if(debug)
        {
            cout << "rover>>>>:" << endl;
            cout << roverData << endl;
        }

        CommonTime epoch = roverData.epoch;

        sppUCCodePhaseRover.solve(roverData);
        EquSys equSysRover = sppUCCodePhaseRover.getEquSys();
        SatID datumSat = sppUCCodePhaseRover.getDatumSat();
        Vector3d xyzRover = sppUCCodePhaseRover.getXYZ();

        // solve spp for base
        ObsData baseObsData;
        try {
            baseObsData = readObsBase.parseRinexObs(epoch);
        }
        // 同步出现错误，跳过基准站处理，读取下一个流动站历元
        catch (SyncException &e) {
            continue;
        };

        if(debug)
        {
            cout << "Base data at epoch:"  << CommonTime2CivilTime(baseObsData.epoch) << endl;
            cout << baseObsData << endl;
        }

        sppUCCodePhaseBase.solve(baseObsData);
        EquSys equSysBase = sppUCCodePhaseBase.getEquSys();

        // compute between-station single-difference equation
        EquSys equSysSD;
        differenceStation(equSysRover, equSysBase, equSysSD);

        // compute between-sat single-difference equation
        EquSys equSysDD;
        differenceSat(datumSat, equSysSD, equSysDD);

        // solve solution
        solverRTK.solve(equSysDD);
        Vector3d dxyzRTK = solverRTK.getdxyz();

        if(debug)
            cout<< dxyzRTK << endl;

        // 获得rtk定位后的接收机位置
        // 其数值应该为接收机单点定位位置+rtk定位后的dxyz；
        Vector3d xyzRTKFloat;
        xyzRTKFloat = xyzRover + dxyzRTK;
        


        // print solution to files
        printSolution(solStream, epoch, xyzRover, xyzRTKFloat);

        // 调试代码时，设置一个stopEpoch，有助于快速得到结果
        if (roverData.epoch > stopEpoch)
            break;

    }

    roverObsStream.close();
    baseObsStream.close();
    solStream.close();

}

