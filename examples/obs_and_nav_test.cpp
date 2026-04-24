//
// Created by 1 on 2026/3/20.
//

#include "obs_and_nav_test.h"
#include "obsreadheader.h"
#include "obsreader.h"
#include "navreadheader.h"
#include "navreader.h"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <iomanip>

int main() {
    std::fstream fs("D:/softwareaddress/RNXCMP_4.2.0_Windows_mingw_64bit/RNXCMP_4.2.0_Windows_mingw_64bit/WUHN00CHN_R_20250010000_01D_30S_MO.25o");
    ObsReadHeader reader;
    reader.setFileStream(&fs);
    reader.parseHeader();
    reader.printHeader();

    ObsReader obsReader;
    obsReader.setFileStream(&fs);
    obsReader.setHeader(&reader.getHeader());

    long long epochCount = 0;
    long long totalSatObs = 0;
    std::set<std::string> uniqueSats;
    std::map<std::string, std::map<std::string, long long>> sysTypeCounts;

    try {
        while (true) {
            ObsEpochData epoch = obsReader.readNextEpoch();
            epochCount++;
            totalSatObs += epoch.numSat;
            if (epochCount <= 1) {
                std::cout << "  --- Detailed Obs (first epochs) ---\n";
                for (const auto& satEntry : epoch.satTypeValue) {
                    const std::string& sat = satEntry.first;
                    std::cout << "  " << sat << "\n";
                    std::string sys = sat.substr(0, 1);
                    auto itTypes = reader.getHeader().sysObsTypes.find(sys);
                    if (itTypes != reader.getHeader().sysObsTypes.end()) {
                        for (const auto& type : itTypes->second) {
                            auto itVal = satEntry.second.find(type);
                            std::cout << "    " << type << " : ";
                            if (itVal == satEntry.second.end() || std::isnan(itVal->second)) {
                                std::cout << "NA\n";
                            } else {
                                std::cout << std::fixed << std::setprecision(3) << itVal->second << "\n";
                            }
                        }
                    } else {
                        for (const auto& tv : satEntry.second) {
                            std::cout << "    " << tv.first << " : ";
                            if (std::isnan(tv.second)) {
                                std::cout << "NA\n";
                            } else {
                                std::cout << std::fixed << std::setprecision(3) << tv.second << "\n";
                            }
                        }
                    }
                }
            }

            for (const auto& satEntry : epoch.satTypeValue) {
                const std::string& sat = satEntry.first;
                uniqueSats.insert(sat);
                std::string sys = sat.substr(0, 1);
                for (const auto& tv : satEntry.second) {
                    if (!std::isnan(tv.second)) {
                        sysTypeCounts[sys][tv.first] += 1;
                    }
                }
            }
        }
    } catch (const EndOfObsFile&) {
        fs.close();
    }

    std::cout << "\n=== Statistics ===\n";
    std::cout << "Epochs: " << epochCount << "\n";
    std::cout << "Observation counts by system/type:\n";
    for (const auto& sysEntry : sysTypeCounts) {
        std::cout << "  System " << sysEntry.first << ":\n";
        for (const auto& typeEntry : sysEntry.second) {
            std::cout << "    " << typeEntry.first << " : " << typeEntry.second << "\n";
        }
    }

    std::cout << "\n=== NAV HEADER TEST ===\n";
    std::fstream navFs("D:/GNSSLAB/gnssLab-2.4/data/BRDC00IGS_R_20250010000_01D_MN.rnx");
    if (!navFs) {
        std::cerr << "Failed to open nav file.\n";
        return 0;
    }

    NavReadHeader navReader;
    navReader.setFileStream(&navFs);
    navReader.parseHeader();
    navReader.printHeader();

    NavReader navDataReader;
    navDataReader.setFileStream(&navFs);
    navDataReader.setHeader(&navReader.getHeader());
    try {
        NavGpsRecord firstRec = navDataReader.readNextGpsRecord();
        std::cout << "\n=== FIRST GPS NAV RECORD ===\n";
        std::cout << "Satellite: " << firstRec.sat << "\n";
        firstRec.eph.printData();
    } catch (const EndOfNavFile& e) {
        std::cerr << e.what() << "\n";
    }

    std::set<SatID> uniqueNavSats = navDataReader.collectUniqueSats();
    std::cout << "\n=== UNIQUE NAV SATELLITES ===\n";
    std::cout << "Total unique satellites in nav file: " << uniqueNavSats.size() << "\n";
    for (const auto& sat : uniqueNavSats) {
        std::cout << "  " << sat << "\n";
    }
    navFs.close();
    return 0;
}
