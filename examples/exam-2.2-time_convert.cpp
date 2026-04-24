
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

#include "TimeStruct.h"
#include "TimeConvert.h"

int main() {

    CivilTime civilTime(2025, 1, 4, 9, 0, 0.0);

    // convert yy/mm/dd to jd
    double jd;

    int yy = 2025;
    int month = 1;
    int day = 1;
    jd = convertYMD2JD(2025, 1, 1);

    cout << "jd for yy/mm/dd: " << yy << "/" << month << "/" << day << " is:" << fixed << jd << endl;

    JulianDate julianDate;
    CommonTime commonTime = CivilTime2CommonTime(civilTime);
    cout << "CommonTime is:" << commonTime << endl;

    julianDate = CommonTime2JulianDate(commonTime);

    cout << "julianData is:" << fixed << julianDate << endl;

    return 0;
}