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
 * 2. 北斗广播星历位置/速度推导，参考教材第 4 章 BDS 16 参数算法。
 */

#include <iostream>
#include <iomanip>

#include "NavEphBDS.hpp"
#include "Const.h"

using namespace std;

#define debug 0

namespace {

class BdsEllipsoid : public WGS84 {
public:
    /// BDS broadcast uses CGCS2000-like constants; use the book's values.
    double angVelocity() const throw() {
        return 7.2921150e-5;
    }

    double gm() const throw() {
        return 3.986004418e14;
    }
};

inline double wrapBdsLikeTime(double tk) {
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    return tk;
}

inline double solveKepler(double Mk, double ecc) {
    double twoPI = 2.0 * PI;
    Mk = fmod(Mk, twoPI);
    double Ek = Mk + ecc * ::sin(Mk);
    int loop_cnt = 1;
    double F, G, delea;
    do {
        F = Mk - (Ek - ecc * ::sin(Ek));
        G = 1.0 - ecc * ::cos(Ek);
        delea = F / G;
        Ek = Ek + delea;
        loop_cnt++;
    } while ((fabs(delea) > 1.0e-11) && (loop_cnt <= 20));
    return Ek;
}

inline Eigen::Matrix3d rotX(double angle) {
    Eigen::Matrix3d R;
    double c = ::cos(angle);
    double s = ::sin(angle);
    R << 1.0, 0.0, 0.0,
         0.0, c,   s,
         0.0, -s,  c;
    return R;
}

inline Eigen::Matrix3d rotZ(double angle) {
    Eigen::Matrix3d R;
    double c = ::cos(angle);
    double s = ::sin(angle);
    R << c,   s,   0.0,
         -s,  c,   0.0,
         0.0, 0.0, 1.0;
    return R;
}

inline Eigen::Matrix3d rotZDot(double angle, double omega) {
    Eigen::Matrix3d R;
    double c = ::cos(angle);
    double s = ::sin(angle);
    R << -s,  c,   0.0,
         -c, -s,   0.0,
          0.0, 0.0, 0.0;
    return omega * R;
}

inline void printVec(const string &label, const Eigen::Vector3d &v) {
    cout << label << ":" << setprecision(12) << v.transpose() << endl;
}

} // namespace

bool NavEphBDS::isGeoSatellite(const SatID &sat) {
    return (sat.system == "C" || sat.system == "B")
           && ((sat.id >= 1 && sat.id <= 5) || (sat.id >= 59 && sat.id <= 62));
}

void NavEphBDS::printData() const {
    cout << "****************************************************************"
         << "************" << endl
         << "BDS Broadcast Ephemeris Data: " << endl;
    cout << "Toc: " << this->CivilToc.year << " " << this->CivilToc.month << " "
         << this->CivilToc.day << " " << this->CivilToc.hour << " "
         << this->CivilToc.minute << " " << this->CivilToc.second << endl;
    cout << fixed << setprecision(12)
         << "af0: " << af0 << endl
         << "af1: " << af1 << endl
         << "af2: " << af2 << endl;

    cout << "IODE/AODE: " << IODE << endl
         << "Crs:  " << Crs << endl
         << "Delta_n: " << Delta_n << endl
         << "M0: " << M0 << endl;

    cout << "Cuc: " << Cuc << endl
         << "ecc: " << ecc << endl
         << "Cus: " << Cus << endl
         << "sqrt_A: " << sqrt_A << endl;

    cout << "Toe: " << Toe << endl;
    cout << "Cic: " << Cic << endl
         << "OMEGA_0: " << OMEGA_0 << endl
         << "Cis: " << Cis << endl;

    cout << "i0: " << i0 << endl
         << "Crc: " << Crc << endl
         << "omega: " << omega << endl
         << "OMEGA_DOT: " << OMEGA_DOT << endl;

    cout << "IDOT: " << IDOT << endl;
    cout << "Codes_On_L2_Channel: " << L2Codes << endl;
    cout << "BDTWeek: " << BDTWeek << endl;
    cout << "L2P_data_flag: " << L2Pflag << endl;

    cout << "URA: " << URA << endl
         << "SV_health: " << SV_health << endl
         << "TGD1: " << TGD1 << endl
         << "TGD2: " << TGD2 << endl;

    cout << "HOWtime: " << HOWtime << endl
         << "fitInterval: " << fitInterval << endl
         << "spare1: " << spare1 << endl
         << "spare2: " << spare2 << endl;

    cout << "ctToc: " << ctToc.toString() << endl
         << "ctToe: " << ctToe.toString() << endl;
}

double NavEphBDS::svClockBias(const CommonTime &t) const {
    CommonTime bdtTime = convertTimeSystem(t, TimeSystem::BDT);
    bdtTime.setTimeSystem(TimeSystem::BDT);
    double elaptc = bdtTime - ctToc;
    if (debug) {
        cout << "elaptc:" << elaptc << endl;
        cout << "af0:" << af0 << "af1:" << af1 << "af2:" << af2 << endl;
    }
    double dtc = af0 + elaptc * (af1 + elaptc * af2);
    if (debug) {
        cout << "dtc:" << dtc << endl;
    }
    return dtc;
}

double NavEphBDS::svClockDrift(const CommonTime &t) const {
    CommonTime bdtTime = convertTimeSystem(t, TimeSystem::BDT);
    bdtTime.setTimeSystem(TimeSystem::BDT);
    double elaptc = bdtTime - ctToc;
    return af1 + elaptc * af2;
}

double NavEphBDS::svRelativity(const CommonTime &t) const {
    BdsEllipsoid ell;
    double A = sqrt_A * sqrt_A;
    CommonTime bdtTime = convertTimeSystem(t, TimeSystem::BDT);
    bdtTime.setTimeSystem(TimeSystem::BDT);
    double tk = wrapBdsLikeTime(bdtTime - ctToe);
    double n0 = std::sqrt(ell.gm() / (A * A * A));
    double n = n0 + Delta_n;
    double Mk = M0 + n * tk;
    double Ek = solveKepler(Mk, ecc);
    return (REL_CONST_BDS * ecc * std::sqrt(A) * ::sin(Ek));
}

double NavEphBDS::svURA(const CommonTime &t) const {
    (void)t;
    return URA;
}

Xvt NavEphBDS::svXvt(const SatID &sat, const CommonTime &t) const {
    Xvt sv;
    BdsEllipsoid ell;
    const bool geo = isGeoSatellite(sat);
    CommonTime bdtTime = convertTimeSystem(t, TimeSystem::BDT);
    bdtTime.setTimeSystem(TimeSystem::BDT);

    double A = sqrt_A * sqrt_A;
    const double omegaE = ell.angVelocity();
    double n0 = std::sqrt(ell.gm() / (A * A * A));
    double tk = wrapBdsLikeTime(bdtTime - ctToe);
    double n = n0 + Delta_n;
    double Mk = M0 + n * tk;
    double Ek = solveKepler(Mk, ecc);

    sv.relcorr = svRelativity(t);
    sv.clkbias = svClockBias(t);
    sv.clkdrift = svClockDrift(t);
    sv.typeTGDData["TGD1"] = TGD1;
    sv.typeTGDData["TGD2"] = TGD2;

    double beta = std::sqrt(1.0 - ecc * ecc);
    double sinEk = ::sin(Ek);
    double cosEk = ::cos(Ek);
    double vk = atan2(beta * sinEk, cosEk - ecc);
    double Phi_k = vk + omega;
    double cos2Phi_k = ::cos(2.0 * Phi_k);
    double sin2Phi_k = ::sin(2.0 * Phi_k);

    double duk = Cuc * cos2Phi_k + Cus * sin2Phi_k;
    double drk = Crc * cos2Phi_k + Crs * sin2Phi_k;
    double dik = Cic * cos2Phi_k + Cis * sin2Phi_k;

    double uk = Phi_k + duk;
    double rk = A * (1.0 - ecc * cosEk) + drk;
    double ik = i0 + dik + IDOT * tk;

    double xk = rk * ::cos(uk);
    double yk = rk * ::sin(uk);

    double Edot = n / (1.0 - ecc * cosEk);
    double vdot = beta * Edot / (1.0 - ecc * cosEk);
    double Phidot = vdot;
    double rdot = A * ecc * sinEk * Edot + 2.0 * (Crs * cos2Phi_k - Crc * sin2Phi_k) * Phidot;
    double udot = vdot + 2.0 * Phidot * (Cus * cos2Phi_k - Cuc * sin2Phi_k);
    double didot = IDOT + 2.0 * Phidot * (Cis * cos2Phi_k - Cic * sin2Phi_k);
    double xdot = rdot * ::cos(uk) - rk * ::sin(uk) * udot;
    double ydot = rdot * ::sin(uk) + rk * ::cos(uk) * udot;

    if (debug) {
        cout << setprecision(12);
        cout << "A:" << A << endl;
        cout << "n0:" << n0 << endl;
        cout << "tk:" << tk << endl;
        cout << "n:" << n << endl;
        cout << "Mk:" << Mk << endl;
        cout << "Ek:" << Ek << endl;
        cout << "vk:" << vk << endl;
        cout << "Phi_k:" << Phi_k << endl;
        cout << "duk:" << duk << endl;
        cout << "drk:" << drk << endl;
        cout << "dik:" << dik << endl;
        cout << "uk:" << uk << endl;
        cout << "rk:" << rk << endl;
        cout << "ik:" << ik << endl;
        cout << "xk:" << xk << endl;
        cout << "yk:" << yk << endl;
        cout << "Edot:" << Edot << endl;
        cout << "vdot:" << vdot << endl;
        cout << "Phidot:" << Phidot << endl;
        cout << "rdot:" << rdot << endl;
        cout << "udot:" << udot << endl;
        cout << "xdot:" << xdot << endl;
        cout << "ydot:" << ydot << endl;
        cout << "didot:" << didot << endl;
    }

    double Omegadot = geo ? OMEGA_DOT : (OMEGA_DOT - omegaE);
    double Omega_k = geo
        ? (OMEGA_0 + Omegadot * tk - omegaE * Toe)
        : (OMEGA_0 + Omegadot* tk - omegaE * Toe);

    double sinOMG_k = ::sin(Omega_k);
    double cosOMG_k = ::cos(Omega_k);
    double cosik = ::cos(ik);
    double sinik = ::sin(ik);

    Eigen::Vector3d posOrb;
    Eigen::Vector3d velOrb;

    posOrb[0] = xk * cosOMG_k - yk * cosik * sinOMG_k;
    posOrb[1] = xk * sinOMG_k + yk * cosik * cosOMG_k;
    posOrb[2] = yk * sinik;

    velOrb[0] = xdot * cosOMG_k - xk * sinOMG_k * Omegadot
                - ydot * cosik * sinOMG_k
                + yk * (sinik * sinOMG_k * didot - cosik * cosOMG_k * Omegadot);
    velOrb[1] = xdot * sinOMG_k + xk * cosOMG_k * Omegadot
                + ydot * cosik * cosOMG_k
                - yk * (sinik * cosOMG_k * didot + cosik * sinOMG_k * Omegadot);
    velOrb[2] = ydot * sinik + yk * cosik * didot;

    if (debug) {
        cout << "OMEGA_k:" << Omega_k << endl;
        cout << "Xk:" << posOrb[0] << endl;
        cout << "Yk:" << posOrb[1] << endl;
        cout << "Zk:" << posOrb[2] << endl;
    }

    if (geo) {
        // The textbook writes Rx(-5°); with the column-vector convention used
        // by this code, the equivalent matrix is rotX(+5°).
        Eigen::Matrix3d R1 = rotX(-5.0 * PI / 180.0);
        Eigen::Matrix3d R2 = rotZ(omegaE * tk);
        Eigen::Matrix3d R2dot = rotZDot(omegaE * tk, omegaE);

        sv.x = R2 * R1 * posOrb;
        sv.v = R2dot * R1 * posOrb + R2 * R1 * velOrb;

        if (debug) {
            cout << "Xgk:" << sv.x[0] << endl;
            cout << "Ygk:" << sv.x[1] << endl;
            cout << "Zgk:" << sv.x[2] << endl;
            cout << "惯性系下速度为:" << endl;
            cout << "(" << velOrb[0] << ", " << velOrb[1] << ", " << velOrb[2] << ")" << endl;
            cout << "地固系下速度为:" << endl;
            cout << "(" << sv.v[0] << ", " << sv.v[1] << ", " << sv.v[2] << ")" << endl;
        }
    } else {
        sv.x = posOrb;
        sv.v = velOrb;
    }

    return sv;
}

bool NavEphBDS::isValid(const CommonTime &ct) const {
    if (ct < beginValid || ct > endValid) return false;
    return true;
}
