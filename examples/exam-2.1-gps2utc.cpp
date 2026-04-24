/**
 * Exercise 2: GPS -> UTC conversion
 */

#include <iostream>
#include <iomanip>
#include "TimeStruct.h"
#include "TimeConvert.h"
#include <vector>
#include <string>
// 包含CivilTime、CommonTime等类定义

// 标准答案对照数组
// 标准答案对照数组
const std::vector<std::pair<std::string, std::string>> expectedResults = {
    {"1993.07.01 00:00:10", "1993.07.01 00:00:01"},
    {"1993.07.01 00:00:09", "1993.07.01 00:00:00"},
    {"1993.07.01 00:00:08", "1993.07.01 00:00:00"},
    {"1993.07.01 00:00:07", "1993.06.30 23:59:59"},
    {"1993.07.01 00:00:06", "1993.06.30 23:59:58"},
    {"1993.07.01 00:00:05", "1993.06.30 23:59:57"},
    {"1993.07.01 00:00:04", "1993.06.30 23:59:56"},
    {"1993.07.01 00:00:03", "1993.06.30 23:59:55"},
    {"1993.07.01 00:00:02", "1993.06.30 23:59:54"},
    {"1993.07.01 00:00:01", "1993.06.30 23:59:53"},
    {"1993.07.01 00:00:00", "1993.06.30 23:59:52"},
    {"1993.06.30 23:59:59", "1993.06.30 23:59:51"}
};

std::string civilToString(const CivilTime &ct) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%04d.%02d.%02d %02d:%02d:%02d",
             ct.year, ct.month, ct.day, ct.hour, ct.minute, (int)ct.second);
    return std::string(buf);
}

int main() {
    bool allCorrect = true;
    int idx = 0;

    CivilTime civil(1993, 7, 1, 0, 0, 10, TimeSystem::GPS);

    // 前11次循环
    for (int i = 0; i < 11; ++i) {
        CommonTime gpsCT = CivilTime2CommonTime(civil);
        CommonTime utcCT = convertTimeSystem(gpsCT, TimeSystem::UTC);
        CivilTime utcCivil = CommonTime2CivilTime(utcCT);

        std::string gpsStr = civilToString(civil);
        std::string utcStr = civilToString(utcCivil);
        std::string expStr = expectedResults[idx].second;

        bool ok = (utcStr == expStr);
        if (!ok) allCorrect = false;

        // 干净输出
        printf("=====================================\n");
        printf("GPS:  %s\n", gpsStr.c_str());
        printf("UTC:  %s\n", utcStr.c_str());
        printf("EXP:  %s\n", expStr.c_str());
        printf("RESULT:  %s\n", ok ? "YES" : "NO");

        idx++;
        civil.second -= 1.0;
    }

    // 最后一个测试用例
    CivilTime gpsCivil(1993, 6, 30, 23, 59, 59, TimeSystem::GPS);
    CommonTime gpsCT = CivilTime2CommonTime(gpsCivil);
    CommonTime utcCT = convertTimeSystem(gpsCT, TimeSystem::UTC);
    CivilTime utcCivil = CommonTime2CivilTime(utcCT);

    std::string gpsStr = civilToString(gpsCivil);
    std::string utcStr = civilToString(utcCivil);
    std::string expStr = expectedResults[idx].second;

    bool ok = (utcStr == expStr);
    if (!ok) allCorrect = false;

    printf("=====================================\n");
    printf("GPS:  %s\n", gpsStr.c_str());
    printf("UTC:  %s\n", utcStr.c_str());
    printf("EXP:  %s\n", expStr.c_str());
    printf("RESULT:  %s\n", ok ? "YES" : "NO");

    // 全部正确才输出
    if (allCorrect) {
        printf("\n=========================\n");
        printf("  All Correct!\n");
        printf("=========================\n");
    }

    return 0;
}