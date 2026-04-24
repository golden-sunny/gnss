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

#include "ConfigReader.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    try {

        string usages =
                " Parse configuration data from file! \n"
                " Usages: \n"
                "    parse_config [your_config_file] \n"
                " Examples: \n"
                "    parse_config ../examples/spp.ini \n";

        if (argc <= 1) {
            cout << "ERROR: you must input config file !" << endl;
            cout << usages << endl;
            exit(-1);
        }

        ConfigReader configReader(argv[1]); // Replace with your config file path

        try {
            int GPS = configReader.getValueAsInt("GPS");
            std::string navFile = configReader.getValueAsString("navFile");
            bool BD2 = configReader.getValueAsBool("BD2");
            double noiseGPS = configReader.getValueAsDouble("noiseGPSCode");

            std::cout << "GPS: " << GPS << std::endl;
            std::cout << "Nav File: " << navFile << std::endl;
            std::cout << "BD2: " << BD2 << std::endl;
            std::cout << "Noise of GPS Code: " << noiseGPS << std::endl;

            double noiseGLO = configReader.getValueAsDouble("noiseGLO");
            cout << noiseGLO << endl;
        }
        catch (std::runtime_error &e) {
            cerr << e.what() << endl;
            exit(-1);
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}