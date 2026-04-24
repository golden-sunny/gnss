/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: Shoujian Zhang
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hern谩ndez-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. 北斗广播星历算法整理自教材第 4 章中关于 BDS/MEO/IGSO/GEO 的位置与速度计算公式。
 */

#ifndef NavEphBDS_HPP
#define NavEphBDS_HPP

#include <string>
#include <cmath>

#include "TimeConvert.h"
#include "GnssStruct.h"

class NavEphBDS {
public:
    NavEphBDS(void)
            : beginValid(END_OF_TIME),
              endValid(BEGINNING_OF_TIME) {
        beginValid.m_timeSystem = TimeSystem::BDT;
        endValid.m_timeSystem = TimeSystem::BDT;
    }

    virtual ~NavEphBDS(void) {}

    void printData() const;

    double svClockBias(const CommonTime &t) const;
    double svClockDrift(const CommonTime &t) const;
    double svRelativity(const CommonTime &t) const;
    double svURA(const CommonTime &t) const;

    /// Compute satellite position/velocity at the given time.
    /// GEO satellites need a special rotation, so SatID is required.
    Xvt svXvt(const SatID &sat, const CommonTime &t) const;

    /// Return true for BDS GEO satellites that require the special 5-degree rotation.
    static bool isGeoSatellite(const SatID &sat);

    bool isValid(const CommonTime &ct) const;

    /// Ephemeris data
    ///   SV/EPOCH/SV CLK
    CivilTime CivilToc;
    double Toc;                ///< Time of clock (BDT seconds of week)
    double af0;                ///< SV clock bias(seconds)
    double af1;                ///< SV clock drift(sec/sec)
    double af2;                ///< SV clock drift rate (sec/sec2)

    ///   BROADCAST ORBIT-1
    double IODE;               ///< Issue of Data, Ephemeris / AODE
    double Crs;                ///< (meters)
    double Delta_n;            ///< Mean Motion Difference From Computed Value(semi-circles/sec)
    double M0;                 ///< Mean Anomaly at Reference Time(semi-circles)

    ///   BROADCAST ORBIT-2
    double Cuc;                ///< (radians)
    double ecc;                ///< Eccentricity
    double Cus;                ///< (radians)
    double sqrt_A;             ///< Square Root of the Semi-Major Axis(sqrt(m))

    ///   BROADCAST ORBIT-3
    double Toe;                ///< Time of Ephemeris(sec of BDS week)
    double Cic;                ///< (radians)
    double OMEGA_0;            ///< Longitude of Ascending Node of Orbit Plane(semi-circles)
    double Cis;                ///< (radians)

    ///   BROADCAST ORBIT-4
    double i0;                 ///< Inclination Angle at Reference Time(semi-circles)
    double Crc;                ///< (meters)
    double omega;              ///< Argument of Perigee(semi-circles)
    double OMEGA_DOT;          ///< Rate of Right Ascension(semi-circles/sec)

    ///   BROADCAST ORBIT-5
    double IDOT;               ///< Rate of Inclination Angle(semi-circles/sec)
    double L2Codes;
    double BDTWeek;            ///< BDS week, continuous number
    double L2Pflag;

    ///   BROADCAST ORBIT-6
    double URA;                ///< SV accuracy(meters)
    double SV_health;          ///< health / usable flag
    double TGD1;               ///< B1/B3 group delay (seconds)
    double TGD2;               ///< B2/B3 group delay (seconds)

    ///   BROADCAST ORBIT-7
    long HOWtime;              ///< Transmission time of message (sec of BDT week)
    double fitInterval;        ///< Fit Interval in hours
    double spare1;
    double spare2;

    /// member data
    CommonTime ctToc;          ///< Toc in CommonTime form
    CommonTime ctToe;          ///< Toe in CommonTime form
    CommonTime transmitTime;   ///< Transmission time in CommonTime form
    CommonTime beginValid;     ///< Time at beginning of validity
    CommonTime endValid;       ///< Time at end of fit validity
};

#endif // NavEphBDS_HPP
