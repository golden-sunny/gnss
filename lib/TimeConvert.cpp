/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: Shoujian Zhang锛宻hjzhang@sgg.whu.edu.cn锛?2024-10-10
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hern谩ndez-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */
#include <cmath>
#include <map>
#include "TimeConvert.h"
#include "TimeStruct.h"
#include "Const.h"

#define debug 0
#define debugYDS 0

using namespace std;

double getLeapSeconds(const CommonTime &ct) {

    std::map<double, double> leapData;

    // mjd;
    //
    leapData[41317.0] = 10;
    leapData[41499.0] = 11;
    leapData[41683.0] = 12;
    leapData[42048.0] = 13;
    leapData[42413.0] = 14;
    leapData[42778.0] = 15;
    leapData[43144.0] = 16;
    leapData[43509.0] = 17;
    leapData[43874.0] = 18;
    leapData[44239.0] = 19;
    leapData[44786.0] = 20;
    leapData[45151.0] = 21;
    leapData[45516.0] = 22;
    leapData[46247.0] = 23;
    leapData[47161.0] = 24;
    leapData[47892.0] = 25;
    leapData[48257.0] = 26;
    leapData[48804.0] = 27;
    leapData[49169.0] = 28;
    leapData[49534.0] = 29;
    leapData[50083.0] = 30;
    leapData[50630.0] = 31;
    leapData[51179.0] = 32;
    leapData[53736.0] = 33;
    leapData[54832.0] = 34;
    leapData[56109.0] = 35;
    leapData[57204.0] = 36;
    leapData[57754.0] = 37;

    double mjd_ct = ct.m_day;

    // 1972.1.1
    if (mjd_ct < 41317.0) {
        InvalidRequest e("Time MUST be greater than 1972 for leapSec!");
        throw (e);
    }

    double leapSec;

    //  寰幆鎵€鏈夌殑璺崇璁板綍
    for (auto ld: leapData) {
        if (ld.first <= mjd_ct)
            leapSec = ld.second;
    }

    return leapSec;
}

CommonTime convertTimeSystem(
        const CommonTime &ct,
        const TimeSystem &outTS) {

    TimeSystem inTS = ct.m_timeSystem;

    double dt(0.0);

    // 1. 时间系统相同，直接返回
    if (inTS == outTS)
        return ct;

    // 2. 第一步：输入时间 → 转成 TAI 原子时（不变）
    if (inTS == TimeSystem::GPS)
        dt = 19.;
    else if (inTS == TimeSystem::UTC)
        dt = getLeapSeconds(ct);
    else if (inTS == TimeSystem::BDT)
        dt = 33.;
    else if (inTS == TimeSystem::TAI)
        dt = 0.;
    else {
        InvalidRequest e("Invalid input TimeSystem " + inTS.toString());
        throw (e);
    }

    CommonTime taiTime = ct + dt;  // 先算出准确的 TAI 时间

    // 3. 第二步：TAI → 目标时间系统（增加 UTC 两次跳秒修正）
    if (outTS == TimeSystem::GPS) {
        dt -= 19.;
    }
    else if (outTS == TimeSystem::UTC) {
        // 关键：UTC 转换必须做两次跳秒计算
        double leap1 = getLeapSeconds(taiTime);
        CommonTime utcTemp = taiTime - leap1;
        double leap2 = getLeapSeconds(utcTemp);
        dt -= leap2;
    }
    else if (outTS == TimeSystem::BDT) {
        dt -= 33.;
    }
    else if (outTS == TimeSystem::TAI) {
        dt -= 0.;
    }
    else {
        InvalidRequest e("Invalid output TimeSystem " + outTS.toString());
        throw (e);
    }

    return ct + dt;
}
//
void convertJD2YMD(double jd,
                   int &iyear,
                   int &imonth,
                   int &iday) {
    double a = std::floor(jd + 0.5);
    double b = a + 1537;
    double c = std::floor((b - 122.1) / 365.25);
    double d = std::floor(365.25 * c);
    double e = std::floor((b - d) / 30.6001);
    iday = b - d - std::floor(30.6001 * e) + (jd + 0.5) - std::floor(jd + 0.5);
    imonth = e - 1 - 12. * std::floor(e / 14);
    iyear = c - 4715 - std::floor((7 + imonth) / 10);
}

double convertYMD2JD(int yy, int mm, int dd) {

    if (mm <= 2) {
        mm += 12;
        yy -= 1;
    }

    double B = 2 - floor(yy / 100) + floor(yy / 400);
    double jd_double = floor(365.25 * (yy + 4716)) + floor(30.6001 * (mm + 1)) + B + dd - 1524.5;
    return jd_double;

}

void convertSOD2HMS(double sod,
                    int &hh,
                    int &mm,
                    double &sec) {
    // Get us to within one day.
    if (sod < 0) {
        sod += (1 +
                static_cast<unsigned long>(sod / SEC_PER_DAY)) * SEC_PER_DAY;
    } else if (sod >= SEC_PER_DAY) {
        sod -= static_cast<unsigned long>(sod / SEC_PER_DAY) * SEC_PER_DAY;
    }

    double temp;               // variable to hold the integer part of sod
    sod = modf(sod, &temp);    // sod holds the fraction, temp the integer
    long seconds = static_cast<long>(temp); // get temp into a real integer

    hh = seconds / 3600;
    mm = (seconds % 3600) / 60;
    sec = double(seconds % 60) + sod;

}

double convertHMS2SOD(int hh,
                      int mm,
                      double sec) {
    return (sec + 60. * (mm + 60. * hh));
}

CommonTime CivilTime2CommonTime(const CivilTime &civilt) {
    CommonTime ct;
    // get the julian day
    double jday = convertYMD2JD(civilt.year, civilt.month, civilt.day);

    // convert jday to mjd day.
    int mjd_day = jday - MJD_TO_JD;

    // get the second of day
    double sod = convertHMS2SOD(civilt.hour, civilt.minute, civilt.second);

    // mjd_day + sod
    ct.set(mjd_day, sod, civilt.timeSys);

    return ct;
}

CivilTime CommonTime2CivilTime(const CommonTime &ct) {
    CivilTime civilt;
    long mjd_day;

    double sod;
    TimeSystem sys;
    // get the julian day, second of day, and fractional second of day
    ct.get(mjd_day, sod, sys);
    civilt.timeSys = sys;

    double jday;
    jday = mjd_day + MJD_TO_JD;

    // convert the julian day to calendar "year/month/day of month"
    convertJD2YMD(jday, civilt.year, civilt.month, civilt.day);

    // convert the (whole) second of day to "hour/minute/second"
    convertSOD2HMS(sod, civilt.hour, civilt.minute, civilt.second);

    return civilt;
}

CommonTime JulianDate2CommonTime(JulianDate &jd) {
    CommonTime ct;
    long double temp_jd(jd.jd);

    long mjd_day = static_cast<long>(temp_jd - MJD_TO_JD);

    long double sod = (temp_jd - std::floor(temp_jd)) * SEC_PER_DAY;

    ct.set(mjd_day,
           static_cast<double>(sod),
           jd.timeSystem);
    return ct;
}

JulianDate CommonTime2JulianDate(CommonTime &ct) {
    JulianDate jd;
    long mjd_day;
    double sod;
    ct.get(mjd_day, sod, jd.timeSystem);

    double jday = mjd_day + MJD_TO_JD;
    jd.jd = static_cast<long double>(jday) +
            (static_cast<long double>(sod)) * DAY_PER_SEC;

    return jd;
};


CommonTime YDSTime2CommonTime(YDSTime &ydst) {
    CommonTime ct;
    long jday = convertYMD2JD(ydst.year, 1, 1) + ydst.doy - 1;
    ct.set(jday, ydst.sod, ydst.timeSystem);
    return ct;
}

YDSTime CommonTime2YDSTime(const CommonTime &ct) {
    YDSTime ydst;
    long mjday;
    double secDay;
    ct.get(mjday, secDay, ydst.timeSystem);
    ydst.sod = static_cast<double>(secDay);

    double jday;
    jday = mjday + MJD_TO_JD;

    int month, day;
    convertJD2YMD(jday, ydst.year, month, day);
    if (debugYDS)
        cout << "year:" << ydst.year << "month:" << month << "day:" << day;

    ydst.doy = jday - convertYMD2JD(ydst.year, 1, 1) + 1 ;
    return ydst;
}

void MJD2CommonTime(MJD& mjd, CommonTime& ct)
{
    try
    {

        long double mday = (mjd.mjd);
        // tmp now holds the partial days
        double sod =  mday - static_cast<long>(mday);
        // convert tmp to seconds of day
        sod *= SEC_PER_DAY;

        ct.set(mday,
               sod,
                mjd.timeSystem);
    }
    catch (InvalidRequest& ip)
    {
        InvalidRequest ir(ip);
        throw(ip);
    }
}

void CommonTime2MJD(const CommonTime& ct, MJD& mjd )
{
    long mday;
    double sod;
    double fsod;
    ct.get(mday, sod, mjd.timeSystem);
    mjd.mjd = static_cast<long double>(mday ) + static_cast<long double>(sod) * DAY_PER_SEC;
}

void CommonTime2WeekSecond(const CommonTime& ct, WeekSecond& wk )
{
    MJD mjd;
    CommonTime2MJD(ct,mjd );
    if (mjd.mjd < wk.MJDEpoch())
    {
        InvalidRequest ir("Unable to convert to Week/Second - before Epoch.");
        throw(ir);
    }

    long mday;
    double sod;
    ct.get(mday, sod, wk.timeSystem);
    // find the number of days since the beginning of the Epoch
    mday -= wk.MJDEpoch();
    // find out how many weeks that is
    wk.week = static_cast<int>(mday / 7);
    // find out what the day of week is
    mday %= 7;

    wk.sow = static_cast<double>(mday * SEC_PER_DAY + sod) ;
}

void WeekSecond2CommonTime(WeekSecond& wk, CommonTime& ct)
{
    try
    {
        //int dow = static_cast<int>( sow * DAY_PER_SEC );
        // Appears to have rounding issues on 32-bit platforms

        int dow = static_cast<int>(wk.sow / SEC_PER_DAY);
        // NB this assumes MJDEpoch is an integer - what if epoch H:M:S != 0:0:0 ?
        long mday = wk.MJDEpoch() + (7 * wk.week) + dow;
        double sod(wk.sow - SEC_PER_DAY * dow);
        ct.set(mday,sod,wk.timeSystem);
    }
    catch (InvalidRequest& ip)
    {
        throw(ip);
    }
}



