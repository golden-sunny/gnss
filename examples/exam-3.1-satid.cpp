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

#include <iostream>
#include <string>
#include "GnssStruct.h"

using namespace std;

int main() {
    SatID sat1("G15");
    SatID sat2("C03");

    cout << "Satellite ID: " << sat1 << endl;
    cout << "Satellite ID: " << sat2 << endl;

    if (sat1 < sat2) {
        cout << "sat1 is less than sat2" << endl;
    } else {
        cout << "sat1 is not less than sat2" << endl;
    }

    if (sat1 == sat2) {
        cout << "sat1 and sat2 are equal" << endl;
    } else {
        cout << "sat1 and sat2 are not equal" << endl;
    }

    return 0;
}//
// Created by shjzh on 2025/2/25.
//
