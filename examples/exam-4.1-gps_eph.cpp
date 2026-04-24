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
#include "GnssStruct.h"
#include "NavEphGPS.hpp"
#include "StringUtils.h"
#include "SP3Store.hpp"
#include "Const.h"

#include <cmath>
#include <iomanip>

struct GpsBookComputed {
    double A{};
    double n0{};
    double tk{};
    double n{};
    double Mk{};
    double Ek{};
    double vk{};
    double Phi_k{};
    double duk{};
    double drk{};
    double dik{};
    double uk{};
    double rk{};
    double ik{};
    double xk{};
    double yk{};
    double Omega_k{};
    double Mdot{};
    double Edot{};
    double vdot{};
    double Phidot{};
    double rdot{};
    double udot{};
    double xdot{};
    double ydot{};
    double didot{};
    double Omegadot{};
    Vector3d ecefPos;
    Vector3d ecefVel;
    double clockBias{};
    double clockDrift{};
    double relcorr{};
    Xvt xvt;
};

struct GpsBookReference {
    double A{};
    double n0{};
    double tk{};
    double n{};
    double Mk{};
    double Ek{};
    double vk{};
    double Phi_k{};
    double uk{};
    double rk{};
    double ik{};
    double xk{};
    double yk{};
    double Omega_k{};
    double Mdot{};
    double Edot{};
    double vdot{};
    double Phidot{};
    double rdot{};
    double udot{};
    double xdot{};
    double ydot{};
    double didot{};
    double Omegadot{};
    Vector3d ecefPos;
    Vector3d ecefVel;
    double clockBias{};
    double clockDrift{};
    double relcorr{};
    Vector3d sp3Pos;
    Vector3d sp3Vel;
};

static GpsBookComputed computeGpsBookValues(const NavEphGPS &gpsEph, const CommonTime &predictedTime) {
    GpsBookComputed out;
    GPSEllipsoid ell;

    out.A = gpsEph.sqrt_A * gpsEph.sqrt_A;
    out.n0 = std::sqrt(ell.gm() / (out.A * out.A * out.A));

    out.tk = predictedTime - gpsEph.ctToe;
    if (out.tk > 302400.0) out.tk -= 604800.0;
    if (out.tk < -302400.0) out.tk += 604800.0;

    out.n = out.n0 + gpsEph.Delta_n;
    out.Mk = gpsEph.M0 + out.n * out.tk;

    double twoPI = 2.0 * PI;
    out.Mk = std::fmod(out.Mk, twoPI);
    out.Ek = out.Mk + gpsEph.ecc * std::sin(out.Mk);
    int loop_cnt = 1;
    double F, G, delea;
    do {
        F = out.Mk - (out.Ek - gpsEph.ecc * std::sin(out.Ek));
        G = 1.0 - gpsEph.ecc * std::cos(out.Ek);
        delea = F / G;
        out.Ek += delea;
        loop_cnt++;
    } while ((std::fabs(delea) > 1.0e-11) && (loop_cnt <= 20));

    double q = std::sqrt(1.0 - gpsEph.ecc * gpsEph.ecc);
    double sinEk = std::sin(out.Ek);
    double cosEk = std::cos(out.Ek);
    out.vk = std::atan2(q * sinEk, cosEk - gpsEph.ecc);
    out.Phi_k = out.vk + gpsEph.omega;

    double cos2phi_k = std::cos(2.0 * out.Phi_k);
    double sin2phi_k = std::sin(2.0 * out.Phi_k);
    out.duk = gpsEph.Cuc * cos2phi_k + gpsEph.Cus * sin2phi_k;
    out.drk = gpsEph.Crc * cos2phi_k + gpsEph.Crs * sin2phi_k;
    out.dik = gpsEph.Cic * cos2phi_k + gpsEph.Cis * sin2phi_k;

    out.uk = out.Phi_k + out.duk;
    out.rk = out.A * (1.0 - gpsEph.ecc * cosEk) + out.drk;
    out.ik = gpsEph.i0 + gpsEph.IDOT * out.tk + out.dik;

    out.xk = out.rk * std::cos(out.uk);
    out.yk = out.rk * std::sin(out.uk);

    out.Omegadot = gpsEph.OMEGA_DOT - ell.angVelocity();
    out.Omega_k = gpsEph.OMEGA_0 + out.Omegadot * out.tk - ell.angVelocity() * gpsEph.Toe;

    double sinOMG_k = std::sin(out.Omega_k);
    double cosOMG_k = std::cos(out.Omega_k);
    double cosik = std::cos(out.ik);
    double sinik = std::sin(out.ik);

    out.ecefPos[0] = out.xk * cosOMG_k - out.yk * cosik * sinOMG_k;
    out.ecefPos[1] = out.xk * sinOMG_k + out.yk * cosik * cosOMG_k;
    out.ecefPos[2] = out.yk * sinik;

    out.Mdot = out.n;
    out.Edot = out.Mdot / (1.0 - gpsEph.ecc * cosEk);
    out.vdot = q * out.Edot / (1.0 - gpsEph.ecc * cosEk);
    out.Phidot = out.vdot;
    out.rdot = gpsEph.sqrt_A * gpsEph.sqrt_A * gpsEph.ecc * out.Edot * sinEk +
               2.0 * (gpsEph.Crs * cos2phi_k - gpsEph.Crc * sin2phi_k) * out.Phidot;
    out.udot = out.vdot + 2.0 * (gpsEph.Cus * cos2phi_k - gpsEph.Cuc * sin2phi_k) * out.Phidot;
    out.xdot = out.rdot * std::cos(out.uk) - out.rk * out.udot * std::sin(out.uk);
    out.ydot = out.rdot * std::sin(out.uk) + out.rk * out.udot * std::cos(out.uk);
    out.didot = gpsEph.IDOT + 2.0 * out.Phidot * (gpsEph.Cis * cos2phi_k - gpsEph.Cic * sin2phi_k);

    out.ecefVel[0] = -out.xk * out.Omegadot * sinOMG_k + out.xdot * cosOMG_k -
                     out.ydot * sinOMG_k * cosik -
                     out.yk * (out.Omegadot * cosOMG_k * cosik -
                               out.didot * sinOMG_k * sinik);
    out.ecefVel[1] = out.xk * out.Omegadot * cosOMG_k + out.xdot * sinOMG_k +
                     out.ydot * cosOMG_k * cosik -
                     out.yk * (out.Omegadot * sinOMG_k * cosik +
                               out.didot * cosOMG_k * sinik);
    out.ecefVel[2] = out.yk * out.didot * cosik + out.ydot * sinik;

    double dtc = predictedTime - gpsEph.ctToc;
    out.clockBias = gpsEph.af0 + dtc * (gpsEph.af1 + dtc * gpsEph.af2);
    out.clockDrift = gpsEph.af1 + dtc * gpsEph.af2;
    out.relcorr = REL_CONST * gpsEph.ecc * std::sqrt(out.A) * std::sin(out.Ek);

    out.xvt.x = out.ecefPos;
    out.xvt.v = out.ecefVel;
    out.xvt.clkbias = out.clockBias;
    out.xvt.clkdrift = out.clockDrift;
    out.xvt.relcorr = out.relcorr;

    return out;
}

static GpsBookReference textbookGpsReference() {
    GpsBookReference ref;
    ref.A = 26560484.172412;
    ref.n0 = 0.000146;
    ref.tk = 300.000000;
    ref.n = 0.000146;
    ref.Mk = 3.069198;
    ref.Ek = 3.070375;
    ref.vk = 3.071543;
    ref.Phi_k = 2.007411;
    ref.uk = 2.007405;
    ref.rk = 26998503.114566;
    ref.ik = 0.966387;
    ref.xk = -11416823.663205;
    ref.yk = 24465798.737636;
    ref.Omega_k = -20.870738;
    ref.Mdot = 0.000146;
    ref.Edot = 0.000143;
    ref.vdot = 0.000141;
    ref.Phidot = 0.000141;
    ref.rdot = 4.535263951190;
    ref.udot = 0.000141140401;
    ref.xdot = -3455.030455452242;
    ref.ydot = -1607.265249938396;
    ref.didot = 0.000000000177;
    ref.Omegadot = -0.000072929173;
    ref.ecefPos = Vector3d(17486772.338748, 4226022.560461, 20131385.866812);
    ref.ecefVel = Vector3d(989.929, 2232.799, -1322.516);
    ref.clockBias = -2.78710952997946e-04;
    ref.clockDrift = 8.981260180010e-12;
    ref.relcorr = 0.0;
    ref.sp3Pos = Vector3d(17486772.311, 4226022.115, 20131386.712);
    ref.sp3Vel = Vector3d(989.934, 2232.825, -1322.526);
    return ref;
}

static void printScalarCompare(const string &name, double calc, double theory, int precision = 12) {
    cout << left << setw(12) << name
         << " calc=" << right << setw(18) << fixed << setprecision(precision) << calc
         << " theory=" << setw(18) << theory
         << " diff=" << setw(18) << (calc - theory)
         << endl;
}

static void printVectorCompare(const string &name, const Vector3d &calc, const Vector3d &theory, int precision = 12) {
    cout << left << setw(12) << name
         << " calc=(" << fixed << setprecision(precision)
         << calc[0] << ", " << calc[1] << ", " << calc[2] << ")"
         << " theory=(" << theory[0] << ", " << theory[1] << ", " << theory[2] << ")"
         << " diff=(" << (calc[0] - theory[0]) << ", "
         << (calc[1] - theory[1]) << ", "
         << (calc[2] - theory[2]) << ")"
         << endl;
}

int main(int argc,char* argv[]) {

    CivilTime civilTimePrediced;
    civilTimePrediced = CivilTime(2025, 1, 1, 0, 5, 0.0);

    CommonTime predictedTime;
    predictedTime = CivilTime2CommonTime(civilTimePrediced);

    YDSTime ydsPrediced;
    ydsPrediced = CommonTime2YDSTime(predictedTime);

    cout << "epoch:" << ydsPrediced << endl;


    string line0=
"G01 2025 01 01 00 00 00 8.645467460160E-06 3.649347490860E-11 0.000000000000E+00";
    string line1=
"     3.900000000000E+01 9.565625000000E+01 4.723053877010E-09 3.125812576130E+00";
    string line2=
"     5.071982741360E-06 2.085076412190E-04 1.234933733940E-06 5.153755249020E+03";
    string line3=
"     2.592000000000E+05 7.636845111850E-08-1.782896012340E+00 8.195638656620E-08";
    string line4=
"     9.596454288100E-01 3.528437500000E+02-1.317975728420E+00-8.442851678660E-09";
    string line5=
"     1.582208762490E-10 1.000000000000E+00 2.347000000000E+03 0.000000000000E+00";
    string line6=
"     2.000000000000E+00 6.300000000000E+01-1.396983861920E-09 3.900000000000E+01";
    string line7=
"     2.520060000000E+05 4.000000000000E+00 0.000000000000E+00 0.000000000000E+00";

    // 将上述rinex 3.05格式导航电文存储到navGPS中
    NavEphGPS gpsEph;

    SatID sat(line0.substr(0,3));

    int yr = safeStoi(line0.substr(4, 4));
    int mo = safeStoi(line0.substr(9, 2));
    int day = safeStoi(line0.substr(12, 2));
    int hr = safeStoi(line0.substr(15, 2));
    int min = safeStoi(line0.substr(18, 2));
    double sec = safeStod(line0.substr(21, 2));

    CivilTime cvt(yr, mo, day, hr, min, sec);
    gpsEph.CivilToc = cvt;
    gpsEph.ctToe = CivilTime2CommonTime(cvt);;
    gpsEph.ctToe.setTimeSystem(TimeSystem::GPS);

    GPSWeekSecond gws;
    CommonTime2WeekSecond(gpsEph.ctToe, gws);     // sow is system-independent
    gpsEph.Toc = gws.sow;
    gpsEph.af0 = safeStod(line0.substr(23, 19));
    gpsEph.af1 = safeStod(line0.substr(42, 19));
    gpsEph.af2 = safeStod(line0.substr(61, 19));

    ///orbit-1
    int n = 4;

    replace(line1.begin(), line1.end(), 'D', 'e');
    gpsEph.IODE = safeStod(line1.substr(n, 19));
    n += 19;
    gpsEph.Crs = safeStod(line1.substr(n, 19));
    n += 19;
    gpsEph.Delta_n = safeStod(line1.substr(n, 19));
    n += 19;
    gpsEph.M0 = safeStod(line1.substr(n, 19));
    ///orbit-2
    n = 4;

    replace(line2.begin(), line2.end(), 'D', 'e');
    gpsEph.Cuc = safeStod(line2.substr(n, 19));
    n += 19;
    gpsEph.ecc = safeStod(line2.substr(n, 19));
    n += 19;
    gpsEph.Cus = safeStod(line2.substr(n, 19));
    n += 19;
    gpsEph.sqrt_A = safeStod(line2.substr(n, 19));

    ///orbit-3
    n = 4;
    replace(line3.begin(), line3.end(), 'D', 'e');
    gpsEph.Toe = safeStod(line3.substr(n, 19));
    n += 19;
    gpsEph.Cic = safeStod(line3.substr(n, 19));
    n += 19;
    gpsEph.OMEGA_0 = safeStod(line3.substr(n, 19));
    n += 19;
    gpsEph.Cis = safeStod(line3.substr(n, 19));

    ///orbit-4
    n = 4;
    replace(line4.begin(), line4.end(), 'D', 'e');
    gpsEph.i0 = safeStod(line4.substr(n, 19));
    n += 19;
    gpsEph.Crc = safeStod(line4.substr(n, 19));
    n += 19;
    gpsEph.omega = safeStod(line4.substr(n, 19));
    n += 19;
    gpsEph.OMEGA_DOT = safeStod(line4.substr(n, 19));

    ///orbit-5
    n = 4;
    replace(line5.begin(), line5.end(), 'D', 'e');
    gpsEph.IDOT = safeStod(line5.substr(n, 19));
    n += 19;
    gpsEph.L2Codes = safeStod(line5.substr(n, 19));
    n += 19;
    gpsEph.GPSWeek = safeStod(line5.substr(n, 19));
    n += 19;
    gpsEph.L2Pflag = safeStod(line5.substr(n, 19));

    ///orbit-6
    n = 4;
    replace(line6.begin(), line6.end(), 'D', 'e');
    gpsEph.URA = safeStod(line6.substr(n, 19));
    n += 19;
    gpsEph.SV_health = safeStod(line6.substr(n, 19));
    n += 19;
    gpsEph.TGD = safeStod(line6.substr(n, 19));
    n += 19;
    gpsEph.IODC = safeStod(line6.substr(n, 19));

    ///orbit-7
    n = 4;
    replace(line7.begin(), line7.end(), 'D', 'e');
    gpsEph.HOWtime = safeStod(line7.substr(n, 19));
    n += 19;
    gpsEph.fitInterval = safeStod(line7.substr(n, 19));

    GPSWeekSecond gws2 = GPSWeekSecond(gpsEph.GPSWeek, gpsEph.Toc, TimeSystem::GPS);
    WeekSecond2CommonTime(gws2, gpsEph.ctToc);
    gpsEph.ctToc.setTimeSystem(TimeSystem::GPS);

    Xvt xvtNav = gpsEph.svXvt(predictedTime);
    cout << "nav:" << xvtNav << endl;

    // 读取精密星历，内插出位置和速度以及钟差等数值

    string sp3File = "D:\\GNSSLAB\\gnssLab-2.4\\data\\COD0MGXFIN_20250010000_01D_05M_ORB.SP3";
    SP3Store sp3Store;
    sp3Store.loadSP3File(sp3File);

    Xvt xvtSP3 = sp3Store.getXvt(sat, predictedTime);
    cout << "sp3:" << xvtSP3 << endl;
    Vector3d xSP3 = xvtSP3.getPos();

    Vector3d diffXYZ = xvtNav.getPos() - xSP3;
    Vector3d diffVel = xvtNav.getVel() - xvtSP3.getVel();
    double diffClockBias = xvtNav.getClockBias() - xvtSP3.getClockBias();
    double diffRelCorr = xvtNav.getRelativityCorr() - xvtSP3.getRelativityCorr();

    cout << ydsPrediced << " \n"
         << " sat:" << sat << " \n"
         << " nav:\n" << xvtNav << " \n"
         << " sp3:\n" << xvtSP3 << " \n"
         << " diffXYZ:\n" << diffXYZ << " \n"
         << " diffVel:\n" << diffVel << " \n"
         << " diffClockBias:\n" << diffClockBias << " \n"
         << " diffRelCorr:\n" << diffRelCorr << " \n"
         << endl;

    cout << "\n============================================================" << endl;
    cout << "Example 4-1 textbook GPS G02" << endl;
    cout << "============================================================" << endl;

    string book0 =
            "G02 2025 01 01 00 00 00-2.787136472760E-04 8.981260180010E-12 0.000000000000E+00";
    string book1 =
            "     6.000000000000E+01 1.311562500000E+02 4.450542525820E-09 3.025440726570E+00";
    string book2 =
            "     6.888061761860E-06 1.654516137210E-02 2.384185791020E-06 5.153686464310E+03";
    string book3 =
            "     2.592000000000E+05 2.030283212660E-07-1.947697258140E+00 3.352761268620E-07";
    string book4 =
            "     9.663872750750E-01 3.332812500000E+02-1.064131669630E+00-8.021762710040E-09";
    string book5 =
            "     1.942938074030E-10 1.000000000000E+00 2.347000000000E+03 0.000000000000E+00";
    string book6 =
            "     2.000000000000E+00 0.000000000000E+00-1.769512891770E-08 6.000000000000E+01";
    string book7 =
            "     2.520060000000E+05 4.000000000000E+00 0.000000000000E+00 0.000000000000E+00";

    NavEphGPS gpsBookEph;
    SatID bookSat(book0.substr(0, 3));

    int bookYr = safeStoi(book0.substr(4, 4));
    int bookMo = safeStoi(book0.substr(9, 2));
    int bookDay = safeStoi(book0.substr(12, 2));
    int bookHr = safeStoi(book0.substr(15, 2));
    int bookMin = safeStoi(book0.substr(18, 2));
    double bookSec = safeStod(book0.substr(21, 2));

    CivilTime bookCvt(bookYr, bookMo, bookDay, bookHr, bookMin, bookSec);
    gpsBookEph.CivilToc = bookCvt;
    gpsBookEph.ctToe = CivilTime2CommonTime(bookCvt);
    gpsBookEph.ctToe.setTimeSystem(TimeSystem::GPS);

    GPSWeekSecond bookGws;
    CommonTime2WeekSecond(gpsBookEph.ctToe, bookGws);
    gpsBookEph.Toc = bookGws.sow;
    gpsBookEph.af0 = safeStod(book0.substr(23, 19));
    gpsBookEph.af1 = safeStod(book0.substr(42, 19));
    gpsBookEph.af2 = safeStod(book0.substr(61, 19));

    int nbook = 4;
    replace(book1.begin(), book1.end(), 'D', 'e');
    gpsBookEph.IODE = safeStod(book1.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Crs = safeStod(book1.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Delta_n = safeStod(book1.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.M0 = safeStod(book1.substr(nbook, 19));

    nbook = 4;
    replace(book2.begin(), book2.end(), 'D', 'e');
    gpsBookEph.Cuc = safeStod(book2.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.ecc = safeStod(book2.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Cus = safeStod(book2.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.sqrt_A = safeStod(book2.substr(nbook, 19));

    nbook = 4;
    replace(book3.begin(), book3.end(), 'D', 'e');
    gpsBookEph.Toe = safeStod(book3.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Cic = safeStod(book3.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.OMEGA_0 = safeStod(book3.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Cis = safeStod(book3.substr(nbook, 19));

    nbook = 4;
    replace(book4.begin(), book4.end(), 'D', 'e');
    gpsBookEph.i0 = safeStod(book4.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.Crc = safeStod(book4.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.omega = safeStod(book4.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.OMEGA_DOT = safeStod(book4.substr(nbook, 19));

    nbook = 4;
    replace(book5.begin(), book5.end(), 'D', 'e');
    gpsBookEph.IDOT = safeStod(book5.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.L2Codes = safeStod(book5.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.GPSWeek = safeStod(book5.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.L2Pflag = safeStod(book5.substr(nbook, 19));

    nbook = 4;
    replace(book6.begin(), book6.end(), 'D', 'e');
    gpsBookEph.URA = safeStod(book6.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.SV_health = safeStod(book6.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.TGD = safeStod(book6.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.IODC = safeStod(book6.substr(nbook, 19));

    nbook = 4;
    replace(book7.begin(), book7.end(), 'D', 'e');
    gpsBookEph.HOWtime = safeStod(book7.substr(nbook, 19));
    nbook += 19;
    gpsBookEph.fitInterval = safeStod(book7.substr(nbook, 19));

    GPSWeekSecond bookGws2 = GPSWeekSecond(gpsBookEph.GPSWeek, gpsBookEph.Toc, TimeSystem::GPS);
    WeekSecond2CommonTime(bookGws2, gpsBookEph.ctToc);
    gpsBookEph.ctToc.setTimeSystem(TimeSystem::GPS);

    cout << "\n--- textbook broadcast ephemeris ---" << endl;
    gpsBookEph.printData();

    GpsBookComputed bookCalc = computeGpsBookValues(gpsBookEph, predictedTime);
    GpsBookReference bookRef = textbookGpsReference();

    cout << "\n--- intermediate values vs textbook theory ---" << endl;
    printScalarCompare("A", bookCalc.A, bookRef.A, 6);
    printScalarCompare("n0", bookCalc.n0, bookRef.n0, 6);
    printScalarCompare("tk", bookCalc.tk, bookRef.tk, 6);
    printScalarCompare("n", bookCalc.n, bookRef.n, 6);
    printScalarCompare("Mk", bookCalc.Mk, bookRef.Mk, 6);
    printScalarCompare("Ek", bookCalc.Ek, bookRef.Ek, 6);
    printScalarCompare("vk", bookCalc.vk, bookRef.vk, 6);
    printScalarCompare("Phi_k", bookCalc.Phi_k, bookRef.Phi_k, 6);
    printScalarCompare("uk", bookCalc.uk, bookRef.uk, 6);
    printScalarCompare("rk", bookCalc.rk, bookRef.rk, 6);
    printScalarCompare("ik", bookCalc.ik, bookRef.ik, 6);
    printScalarCompare("Omega_k", bookCalc.Omega_k, bookRef.Omega_k, 6);
    printVectorCompare("XYZk", bookCalc.ecefPos, bookRef.ecefPos, 6);

    cout << "\n--- velocity intermediates vs textbook theory ---" << endl;
    printScalarCompare("Mdot", bookCalc.Mdot, bookRef.Mdot, 9);
    printScalarCompare("Edot", bookCalc.Edot, bookRef.Edot, 9);
    printScalarCompare("vdot", bookCalc.vdot, bookRef.vdot, 9);
    printScalarCompare("Phidot", bookCalc.Phidot, bookRef.Phidot, 9);
    printScalarCompare("rdot", bookCalc.rdot, bookRef.rdot, 9);
    printScalarCompare("udot", bookCalc.udot, bookRef.udot, 12);
    printScalarCompare("xdot", bookCalc.xdot, bookRef.xdot, 6);
    printScalarCompare("ydot", bookCalc.ydot, bookRef.ydot, 6);
    printScalarCompare("didot", bookCalc.didot, bookRef.didot, 12);
    printScalarCompare("Omegadot", bookCalc.Omegadot, bookRef.Omegadot, 12);
    printVectorCompare("XYZdot", bookCalc.ecefVel, bookRef.ecefVel, 6);

    Xvt bookNav = bookCalc.xvt;
    string bookSp3File = "D:\\GNSSLAB\\gnssLab-2.4\\data\\WUM0MGXFIN_20250010000_01D_05M_ORB.SP3";
    SP3Store bookSp3Store;
    bookSp3Store.loadSP3File(bookSp3File);

    Xvt bookSp3 = bookSp3Store.getXvt(bookSat, predictedTime);
    Vector3d diffBookXYZ = bookNav.getPos() - bookSp3.getPos();
    Vector3d diffBookVel = bookNav.getVel() - bookSp3.getVel();
    double diffBookClockBias = bookNav.getClockBias() - bookSp3.getClockBias();
    double diffBookRelCorr = bookNav.getRelativityCorr() - bookSp3.getRelativityCorr();

    cout << "\n--- final result vs precise ephemeris ---" << endl;
    printVectorCompare("nav-pos", bookNav.getPos(), bookRef.ecefPos, 6);
    printVectorCompare("sp3-pos", bookSp3.getPos(), bookRef.sp3Pos, 6);
    printVectorCompare("nav-vel", bookNav.getVel(), bookRef.ecefVel, 6);
    printVectorCompare("sp3-vel", bookSp3.getVel(), bookRef.sp3Vel, 6);
    printScalarCompare("clockBias", bookNav.getClockBias(), bookRef.clockBias, 12);
    printScalarCompare("clockDrift", bookNav.getClockDrift(), bookRef.clockDrift, 15);
    printScalarCompare("relCorr", bookNav.getRelativityCorr(), bookSp3.getRelativityCorr(), 12);

    cout << ydsPrediced << " \n"
         << " sat:" << bookSat << " \n"
         << " nav:\n" << bookNav << " \n"
         << " sp3:\n" << bookSp3 << " \n"
         << " diffXYZ:\n" << diffBookXYZ << " \n"
         << " diffVel:\n" << diffBookVel << " \n"
         << " diffClockBias:\n" << diffBookClockBias << " \n"
         << " diffRelCorr:\n" << diffBookRelCorr << " \n"
         << endl;

}
