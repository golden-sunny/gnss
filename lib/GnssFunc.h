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

#ifndef GNSSLAB_GNSSFUNC_H
#define GNSSLAB_GNSSFUNC_H

#include <string>
#include <set>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "CoordStruct.h"
#include "GnssStruct.h"
#include "StringUtils.h"
#include "RinexNavStore.hpp"

using namespace Eigen;

// Changing Sign
inline double sign(double x)
{
    return (x <= 0.0) ? -1.0 : 1.0;
};

// Rounding Values
inline double round(double x)
{
    return double(std::floor(x + 0.5));
};

// Swapping values
inline void swap(double& a, double& b)
{
    double t(a); a = b; b = t;
};

void parseRinexHeader(std::fstream &rinexFileStream,
                      RinexHeader &rinexHeader);

ObsData parseRinexObs(std::fstream &rinexFileStream);

CommonTime parseTime(const string &line);

void chooseObs(ObsData &obsData,
               std::map<string, std::set<string>> &sysTypes);

void convertObsType(ObsData &obsData);

struct ObsDelayCorrectionRecord {
    SatID sat;
    std::string obsType;
    double tgdSeconds = 0.0;
    double correctionSeconds = 0.0;
    double correctionMeters = 0.0;
};

struct SignalDelayTerm {
    double referenceDelaySeconds = 0.0;
    double iscSeconds = 0.0;
    bool hasReferenceDelay = false;
    bool hasIsc = false;
};

using ObsDelayBiasTable = std::map<SatID, std::map<std::string, SignalDelayTerm>>;

// Apply the textbook GPS code-delay correction layer (5.4.1) to code-like
// observations. The correction uses the broadcast TGD stored in the satellite
// ephemeris/Xvt and can be used as a standalone preprocessing step.
void correctTGD(ObsData &obsData,
                const std::map<SatID, Xvt> &satXvtTransTime,
                const ObsDelayBiasTable *biasTable = nullptr,
                std::vector<ObsDelayCorrectionRecord> *report = nullptr);

double wavelengthOfMW(string sys,
                      string L1Type,
                      string L2Type);

double varOfMW(string,
               string L1Type,
               string L2Type);

//================
// 系统误差相关函数
//================
std::map<SatID,Xvt> computeSatPos(ObsData &obsData, RinexNavStore& navStore);
Xvt computeAtTransmitTime(const CommonTime& tr,
                          const double& pr,
                          const SatID& sat,
                          RinexNavStore& navStore);

void computeElevAzim(Eigen::Vector3d& xyz,
                     std::map<SatID,Xvt> & satXvtTransTime,
                     SatValueMap& tempElevData,
                     SatValueMap& tempAzimData);

std::map<SatID,Xvt> earthRotation(Eigen::Vector3d& xyz,
                                  std::map<SatID,Xvt> & satXvtTransTime);

// todo
// ===========
// computeIonoDelay();
// computeTropDelay();
//================

void detectCSMW(ObsData &obsData,
                std::map<Variable, int> &csFlagData,
                SatEpochValueMap &satEpochMWData,
                SatEpochValueMap &satEpochMeanMWData,
                SatEpochValueMap &satEpochCSFlagData);

void differenceStation(EquSys& equSysRover, VariableDataMap& csFlagRover,
                       EquSys& equSysBase, VariableDataMap& csFlagBase,
                       EquSys& equSysSD, VariableDataMap& csFlagSD);

void differenceStation(EquSys& equSysRover,
                       EquSys& equSysBase,
                       EquSys& equSysSD);

SatID findDatumSat(bool& firstEpoch,
                   SatValueMap& satElevData);

void differenceSat( SatID& datumSat,
                    EquSys& equSysSD, VariableDataMap& csFlagSD,
                    EquSys& equSysDD, VariableDataMap& csFlagDD);

void differenceSat( SatID& datumSat,
                    EquSys& equSysSD,
                    EquSys& equSysDD);

void ambiguityDatum(bool& firstEpoch,
                    SatID& datumSat,
                    VariableDataMap& fixedAmbData,
                    EquSys& equSysDD);

void fixSolution(VectorXd& stateVec,
                 MatrixXd& covMatrix,
                 VariableSet& varSet,
                 double& ratio,
                 Vector3d& dxyzFixed,
                 VariableDataMap& fixedAmbData);

// print solution to files
void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover,
                   Eigen::Vector3d& xyzRTKFloat,
                   double& ratio,
                   Eigen::Vector3d& xyzRTKFixed);

// print solution to files
void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover,
                   Eigen::Vector3d& xyzRTKFloat);

// print solution to files
void printSolution(std::fstream & solStream,
                   CommonTime& ctTime,
                   Eigen::Vector3d& xyzRover);

#endif //GNSSLAB_GNSSFUNC_H
