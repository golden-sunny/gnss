//
// BDS broadcast ephemeris walkthrough for C01 and C06.
// The program reads the navigation message from a RINEX BRDC file, computes
// the satellite position/velocity at the textbook epoch, and prints the
// intermediate values together with the book's theoretical numbers.
//

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "NavEphBDS.hpp"
#include "SP3Store.hpp"
#include "navreadheader.h"
#include "navreader.h"
#include "TimeConvert.h"

using namespace std;
using namespace Eigen;

namespace {

struct BdsComputation {
    bool geo{false};
    CommonTime targetGps;
    CommonTime targetBdt;
    double omegaE{0.0};
    double A{0.0};
    double n0{0.0};
    double tk{0.0};
    double n{0.0};
    double Mk{0.0};
    double Ek{0.0};
    double vk{0.0};
    double Phi_k{0.0};
    double duk{0.0};
    double drk{0.0};
    double dik{0.0};
    double uk{0.0};
    double rk{0.0};
    double ik{0.0};
    double xk{0.0};
    double yk{0.0};
    double Omegadot{0.0};
    double Omega_k{0.0};
    double Edot{0.0};
    double vdot{0.0};
    double Phidot{0.0};
    double rdot{0.0};
    double udot{0.0};
    double didot{0.0};
    double xdot{0.0};
    double ydot{0.0};
    Vector3d posOrb{};
    Vector3d velOrb{};
    Vector3d finalPos{};
    Vector3d finalVel{};
    Xvt navXvt{};
};

struct ScalarRow {
    string formula;
    double actual{0.0};
    double theory{0.0};
    int precision{6};
};

struct VectorRow {
    string formula;
    Vector3d actual{};
    Vector3d theory{};
    int precision{6};
};

static string fmt(double v, int precision) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << v;
    return oss.str();
}

static string fmt(const Vector3d &v, int precision) {
    ostringstream oss;
    oss << fixed << setprecision(precision)
        << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return oss.str();
}

static void printScalarRow(const ScalarRow &row) {
    cout << "  " << row.formula << "\n";
    cout << "    actual : " << fmt(row.actual, row.precision) << "\n";
    cout << "    theory : " << fmt(row.theory, row.precision) << "\n";
}

static void printVectorRow(const VectorRow &row) {
    cout << "  " << row.formula << "\n";
    cout << "    actual : " << fmt(row.actual, row.precision) << "\n";
    cout << "    theory : " << fmt(row.theory, row.precision) << "\n";
}

static vector<ScalarRow> buildC01PositionRows(const BdsComputation& c) {
    return {
        {"A = (sqrt_A)^2", c.A, 42166175.565771, 6},
        {"n0 = sqrt(mu / A^3)", c.n0, 0.000073, 6},
        {"tk = t - toe", c.tk, 286.000000, 6},
        {"n = n0 + delta_n", c.n, 0.000073, 6},
        {"Mk = M0 + n*tk", c.Mk, -0.682851, 6},
        {"Ek = Mk + e*sin(Ek)", c.Ek, -0.683150, 6},
        {"vk = atan2(sqrt(1-e^2)sinEk, cosEk-e)", c.vk, -0.683449, 6},
        {"Phi_k = vk + omega", c.Phi_k, -0.869285, 6},
        {"delta_u_k = Cuc*cos(2Phi_k) + Cus*sin(2Phi_k)", c.duk, -0.000028, 6},
        {"delta_r_k = Crc*cos(2Phi_k) + Crs*sin(2Phi_k)", c.drk, 157.336631, 6},
        {"delta_i_k = Cic*cos(2Phi_k) + Cis*sin(2Phi_k)", c.dik, 0.000000, 6},
        {"u_k = Phi_k + delta_u_k", c.uk, -0.869314, 6},
        {"r_k = A(1 - e*cosEk) + delta_r_k", c.rk, 42150835.964427, 6},
        {"i_k = i0 + IDOT*tk + delta_i_k", c.ik, 0.058372, 6},
        {"x_k = r_k*cos(u_k)", c.xk, 27202087.124328, 6},
        {"y_k = r_k*sin(u_k)", c.yk, -32198438.294124, 6},
        {"Omega_k = OMEGA_0 + OMEGA_DOT*tk - omega_e*Toe", c.Omega_k, -21.721757, 6},
    };
}

static vector<ScalarRow> buildC01VelocityRows(const BdsComputation& c) {
    return {
        {"Mdot = n", c.n, 0.000072916034, 12},
        {"Edot = Mdot / (1 - e*cosEk)", c.Edot, 0.000072942842, 12},
        {"vdot = sqrt(1-e^2)*Edot / (1 - e*cosEk)", c.vdot, 0.000072969028, 12},
        {"Phidot = vdot", c.Phidot, 0.000072969028, 12},
        {"rdot = A*e*sinEk*Edot + 2(...) * Phidot", c.rdot, -1.045269503484, 12},
        {"udot = vdot + 2*Phidot*(Cus*cos2Phi - Cuc*sin2Phi)", c.udot, 0.000072968309, 12},
        {"xdot = rdot*cos(u_k) - r_k*udot*sin(u_k)", c.xdot, -1736.518640042023, 12},
        {"ydot = rdot*sin(u_k) + r_k*udot*cos(u_k)", c.ydot, -2535.962486858447, 12},
        {"didot = IDOT + 2*Phidot*(Cis*cos2Phi - Cic*sin2Phi)", c.didot, 0.000000000000, 12},
        {"Omega_dot_k = OMEGA_DOT", c.Omegadot, 0.000000000000, 12},
    };
}

static vector<ScalarRow> buildC06PositionRows(const BdsComputation& c) {
    return {
        {"A = (sqrt_A)^2", c.A, 42163633.240530, 6},
        {"n0 = sqrt(mu / A^3)", c.n0, 0.000073, 6},
        {"tk = t - toe", c.tk, 286.000000, 6},
        {"n = n0 + delta_n", c.n, 0.000073, 6},
        {"Mk = M0 + n*tk", c.Mk, -2.833020, 6},
        {"Ek = Mk + e*sin(Ek)", c.Ek, -2.834316, 6},
        {"vk = atan2(sqrt(1-e^2)sinEk, cosEk-e)", c.vk, -2.835610, 6},
        {"Phi_k = vk + omega", c.Phi_k, -5.554115, 6},
        {"delta_u_k = Cuc*cos(2Phi_k) + Cus*sin(2Phi_k)", c.duk, 0.000011, 6},
        {"delta_r_k = Crc*cos(2Phi_k) + Crs*sin(2Phi_k)", c.drk, 16.729332, 6},
        {"delta_i_k = Cic*cos(2Phi_k) + Cis*sin(2Phi_k)", c.dik, 0.000000, 6},
        {"u_k = Phi_k + delta_u_k", c.uk, -5.554104, 6},
        {"r_k = A(1 - e*cosEk) + delta_r_k", c.rk, 42335915.623071, 6},
        {"i_k = i0 + IDOT*tk + delta_i_k", c.ik, 0.946266, 6},
        {"x_k = r_k*cos(u_k)", c.xk, 31573575.262918, 6},
        {"y_k = r_k*sin(u_k)", c.yk, 28203529.863487, 6},
        {"Omega_k = OMEGA_0 + (OMEGA_DOT - omega_e)*tk - omega_e*Toe", c.Omega_k, -17.703253, 6},
    };
}

static vector<ScalarRow> buildC06VelocityRows(const BdsComputation& c) {
    return {
        {"Mdot = n", c.n, 0.000072924013, 12},
        {"Edot = Mdot / (1 - e*cosEk)", c.Edot, 0.000072627255, 12},
        {"vdot = sqrt(1-e^2)*Edot / (1 - e*cosEk)", c.vdot, 0.000072329589, 12},
        {"Phidot = vdot", c.Phidot, 0.000072329589, 12},
        {"rdot = A*e*sinEk*Edot + 2(...) * Phidot", c.rdot, -3.955580183010, 12},
        {"udot = vdot + 2*Phidot*(Cus*cos2Phi - Cuc*sin2Phi)", c.udot, 0.000072329616, 12},
        {"xdot = rdot*cos(u_k) - r_k*udot*sin(u_k)", c.xdot, -2042.900493825388, 12},
        {"ydot = rdot*sin(u_k) + r_k*udot*cos(u_k)", c.ydot, 2281.069415742512, 12},
        {"didot = IDOT + 2*Phidot*(Cis*cos2Phi - Cic*sin2Phi)", c.didot, 0.000000011450, 12},
        {"Omega_dot_k = OMEGA_DOT - omega_e", c.Omegadot, -0.000072922848, 12},
    };
}

static string trim(const string &s) {
    size_t b = s.find_first_not_of(' ');
    if (b == string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

static double wrapBdsLikeTime(double tk) {
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    return tk;
}

static double solveKepler(double Mk, double ecc) {
    const double twoPI = 2.0 * PI;
    Mk = fmod(Mk, twoPI);
    double Ek = Mk + ecc * sin(Mk);
    int loop_cnt = 1;
    double F = 0.0, G = 0.0, delea = 0.0;
    do {
        F = Mk - (Ek - ecc * sin(Ek));
        G = 1.0 - ecc * cos(Ek);
        delea = F / G;
        Ek += delea;
        ++loop_cnt;
    } while (fabs(delea) > 1.0e-11 && loop_cnt <= 20);
    return Ek;
}

static Matrix3d rotX(double angle) {
    Matrix3d R;
    const double c = cos(angle);
    const double s = sin(angle);
    R << 1.0, 0.0, 0.0,
         0.0, c,   s,
         0.0, -s,  c;
    return R;
}

static Matrix3d rotZ(double angle) {
    Matrix3d R;
    const double c = cos(angle);
    const double s = sin(angle);
    R << c,   s,   0.0,
         -s,  c,   0.0,
          0.0, 0.0, 1.0;
    return R;
}

static Matrix3d rotZDot(double angle, double omega) {
    Matrix3d R;
    const double c = cos(angle);
    const double s = sin(angle);
    R << -s,  c,   0.0,
         -c, -s,   0.0,
          0.0, 0.0, 0.0;
    return omega * R;
}

template <typename EphT>
static bool selectNearestEph(const map<CommonTime, EphT> &records,
                             const CommonTime &target,
                             EphT &selected,
                             double &absDt,
                             double maxAgeSec = 7201.0) {
    if (records.empty()) return false;
    absDt = numeric_limits<double>::max();
    bool found = false;
    for (const auto &it : records) {
        const double d = fabs(target - it.first);
        if (d < absDt) {
            absDt = d;
            selected = it.second;
            found = true;
        }
    }
    return found && absDt <= maxAgeSec;
}

static map<SatID, map<CommonTime, NavEphBDS>> loadBdsNavData(const string &navFile) {
    fstream navFs(navFile.c_str(), ios::in);
    if (!navFs.is_open()) {
        throw runtime_error("Cannot open navigation file: " + navFile);
    }

    NavReadHeader headerReader;
    headerReader.setFileStream(&navFs);
    headerReader.parseHeader();

    NavReader reader;
    reader.setFileStream(&navFs);
    reader.setHeader(&headerReader.getHeader());

    map<SatID, map<CommonTime, NavEphBDS>> data;
    try {
        while (true) {
            NavBdsRecord rec = reader.readNextBdsRecord();
            data[rec.sat][rec.eph.ctToe] = rec.eph;
        }
    } catch (const EndOfNavFile &) {
    }

    navFs.close();
    return data;
}

static BdsComputation computeBdsComputation(const NavEphBDS &eph,
                                            const SatID &sat,
                                            const CommonTime &targetGps) {
    BdsComputation out;
    out.geo = NavEphBDS::isGeoSatellite(sat);
    out.targetGps = targetGps;
    out.targetBdt = convertTimeSystem(targetGps, TimeSystem::BDT);
    out.targetBdt.setTimeSystem(TimeSystem::BDT);
    out.omegaE = 7.2921150e-5;
    out.A = eph.sqrt_A * eph.sqrt_A;
    out.n0 = sqrt(3.986004418e14 / (out.A * out.A * out.A));
    out.tk = wrapBdsLikeTime(out.targetBdt - eph.ctToe);
    out.n = out.n0 + eph.Delta_n;
    out.Mk = eph.M0 + out.n * out.tk;
    out.Ek = solveKepler(out.Mk, eph.ecc);

    out.navXvt = eph.svXvt(sat, targetGps);

    const double beta = sqrt(1.0 - eph.ecc * eph.ecc);
    const double sinEk = sin(out.Ek);
    const double cosEk = cos(out.Ek);
    out.vk = atan2(beta * sinEk, cosEk - eph.ecc);
    out.Phi_k = out.vk + eph.omega;
    const double cos2Phi_k = cos(2.0 * out.Phi_k);
    const double sin2Phi_k = sin(2.0 * out.Phi_k);

    out.duk = eph.Cuc * cos2Phi_k + eph.Cus * sin2Phi_k;
    out.drk = eph.Crc * cos2Phi_k + eph.Crs * sin2Phi_k;
    out.dik = eph.Cic * cos2Phi_k + eph.Cis * sin2Phi_k;
    out.uk = out.Phi_k + out.duk;
    out.rk = out.A * (1.0 - eph.ecc * cosEk) + out.drk;
    out.ik = eph.i0 + out.dik + eph.IDOT * out.tk;
    out.xk = out.rk * cos(out.uk);
    out.yk = out.rk * sin(out.uk);

    out.Omegadot = out.geo ? eph.OMEGA_DOT : (eph.OMEGA_DOT - out.omegaE);
    out.Omega_k = out.geo
        ? (eph.OMEGA_0 + out.Omegadot * out.tk - out.omegaE * eph.Toe)
        : (eph.OMEGA_0 + out.Omegadot * out.tk - out.omegaE * eph.Toe);

    out.Edot = out.n / (1.0 - eph.ecc * cosEk);
    out.vdot = beta * out.Edot / (1.0 - eph.ecc * cosEk);
    out.Phidot = out.vdot;
    out.rdot = out.A * eph.ecc * sinEk * out.Edot
               + 2.0 * (eph.Crs * cos2Phi_k - eph.Crc * sin2Phi_k) * out.Phidot;
    out.udot = out.vdot + 2.0 * out.Phidot * (eph.Cus * cos2Phi_k - eph.Cuc * sin2Phi_k);
    out.didot = eph.IDOT + 2.0 * out.Phidot * (eph.Cis * cos2Phi_k - eph.Cic * sin2Phi_k);
    out.xdot = out.rdot * cos(out.uk) - out.rk * sin(out.uk) * out.udot;
    out.ydot = out.rdot * sin(out.uk) + out.rk * cos(out.uk) * out.udot;

    out.posOrb = Vector3d(out.xk * cos(out.Omega_k) - out.yk * cos(out.ik) * sin(out.Omega_k),
                          out.xk * sin(out.Omega_k) + out.yk * cos(out.ik) * cos(out.Omega_k),
                          out.yk * sin(out.ik));

    out.velOrb = Vector3d(out.xdot * cos(out.Omega_k)
                          - out.xk * sin(out.Omega_k) * out.Omegadot
                          - out.ydot * cos(out.ik) * sin(out.Omega_k)
                          + out.yk * (sin(out.ik) * sin(out.Omega_k) * out.didot
                                      - cos(out.ik) * cos(out.Omega_k) * out.Omegadot),
                          out.xdot * sin(out.Omega_k)
                          + out.xk * cos(out.Omega_k) * out.Omegadot
                          + out.ydot * cos(out.ik) * cos(out.Omega_k)
                          - out.yk * (sin(out.ik) * cos(out.Omega_k) * out.didot
                                      + cos(out.ik) * sin(out.Omega_k) * out.Omegadot),
                          out.ydot * sin(out.ik) + out.yk * cos(out.ik) * out.didot);

    if (out.geo) {
        const Matrix3d R1 = rotX(-5.0 * PI / 180.0);
        const Matrix3d R2 = rotZ(out.omegaE * out.tk);
        const Matrix3d R2dot = rotZDot(out.omegaE * out.tk, out.omegaE);
        out.finalPos = R2 * R1 * out.posOrb;
        out.finalVel = R2dot * R1 * out.posOrb + R2 * R1 * out.velOrb;
    } else {
        out.finalPos = out.posOrb;
        out.finalVel = out.velOrb;
    }

    return out;
}

static void printHeader(const string &title, const string &satName, bool geo) {
    cout << "\n============================================================\n";
    cout << title << " (" << satName << ")"
         << "  [" << (geo ? "GEO" : "IGSO/MEO") << "]\n";
    cout << "============================================================\n";
}

static void printCase(const string &tableTitlePos,
                      const string &tableTitleVel,
                      const BdsComputation &c,
                      const vector<ScalarRow> &posRows,
                      const vector<ScalarRow> &velRows,
                      const Vector3d &theoryFinalPos,
                      const Vector3d &theoryFinalVel,
                      SP3Store &sp3,
                      const SatID &sat,
                      const CommonTime &targetGps) {
    printHeader(tableTitlePos, sat.toString(), c.geo);

    cout << "Target GPS time : " << CommonTime2YDSTime(targetGps) << "\n";
    cout << "Target BDT time : " << CommonTime2CivilTime(c.targetBdt).toString() << "\n";
    cout << "Branch          : " << (c.geo ? "special GEO rotation" : "IGSO/MEO") << "\n";
    cout << "------------------------------------------------------------\n";

    Xvt navXvt = c.navXvt;
    Xvt precise = sp3.getXvt(sat, targetGps);

    for (const auto &row : posRows) {
        printScalarRow(row);
    }

    cout << "  Final ECEF position/velocity from NavEphBDS::svXvt\n";
    cout << "    x actual : " << fmt(navXvt.getPos(), 12) << "\n";
    cout << "    x theory : " << fmt(theoryFinalPos, 12) << "\n";
    cout << "    v actual : " << fmt(navXvt.getVel(), 12) << "\n";
    cout << "    v theory : " << fmt(theoryFinalVel, 12) << "\n";

    cout << "\n" << tableTitleVel << "\n";
    cout << "------------------------------------------------------------\n";
    for (const auto &row : velRows) {
        printScalarRow(row);
    }

    const Vector3d posDiff = navXvt.getPos() - precise.getPos();
    const Vector3d velDiff = navXvt.getVel() - precise.getVel();

    cout << "\nSP3 comparison\n";
    cout << "  precise pos : " << fmt(precise.getPos(), 3) << "\n";
    cout << "  precise vel : " << fmt(precise.getVel(), 3) << "\n";
    cout << "  diff pos    : " << fmt(posDiff, 3) << "\n";
    cout << "  diff vel    : " << fmt(velDiff, 3) << "\n";
}

} // namespace

int main(int argc, char *argv[]) {
    const string navFile = (argc > 1) ? argv[1] : R"(D:\GNSSLAB\gnssLab-2.4\data\BRDC00IGS_R_20250010000_01D_MN.rnx)";
    const string sp3File = (argc > 2) ? argv[2] : R"(D:\GNSSLAB\gnssLab-2.4\data\WUM0MGXFIN_20250010000_01D_05M_ORB.SP3)";

    const CivilTime targetCivil(2025, 1, 1, 0, 5, 0.0, TimeSystem::GPS);
    CommonTime targetGps = CivilTime2CommonTime(targetCivil);
    targetGps.setTimeSystem(TimeSystem::GPS);
    CommonTime targetBdt = convertTimeSystem(targetGps, TimeSystem::BDT);
    targetBdt.setTimeSystem(TimeSystem::BDT);

    cout << "nav file : " << navFile << "\n";
    cout << "sp3 file : " << sp3File << "\n";
    cout << "target   : " << CommonTime2YDSTime(targetGps) << "\n";

    const auto navData = loadBdsNavData(navFile);
    SP3Store sp3;
    sp3.loadSP3File(sp3File);

    const auto printOne = [&](const string &satName,
                              const string &posTitle,
                              const string &velTitle,
                              const Vector3d &theoryFinalPos,
                              const Vector3d &theoryFinalVel) {
        const SatID sat(satName);
        auto it = navData.find(sat);
        if (it == navData.end()) {
            cout << "\nSatellite " << satName << " not found in nav file.\n";
            return;
        }

        NavEphBDS eph;
        double age = 0.0;
        if (!selectNearestEph(it->second, targetBdt, eph, age)) {
            cout << "\nSatellite " << satName << " has no ephemeris close to the target time.\n";
            return;
        }

        BdsComputation c = computeBdsComputation(eph, sat, targetGps);
        cout << "\nSatellite " << satName << " selected ephemeris age: " << fmt(age, 3) << " s\n";
        cout << "Broadcast epoch   : " << CommonTime2CivilTime(eph.ctToe).toString() << "\n";
        cout << "Clock epoch       : " << CommonTime2CivilTime(eph.ctToc).toString() << "\n";
        cout << "Computed branch   : " << (c.geo ? "GEO" : "IGSO/MEO") << "\n";
        vector<ScalarRow> posRows = (satName == "C01")
                                    ? buildC01PositionRows(c)
                                    : buildC06PositionRows(c);
        vector<ScalarRow> velRows = (satName == "C01")
                                    ? buildC01VelocityRows(c)
                                    : buildC06VelocityRows(c);
        printCase(posTitle, velTitle, c, posRows, velRows, theoryFinalPos, theoryFinalVel, sp3, sat, targetGps);
    };

    printOne(
        "C01",
        "Table 4-8 C01 satellite position intermediate values",
        "Table 4-9 C01 satellite velocity intermediate values",
        Vector3d(34271596.501952, 24537957.549023, 198204.514906),
        Vector3d(0.300914383807, -0.522354726871, -105.591181578189)
    );

    printOne(
        "C06",
        "Table 4-10 C06 satellite position intermediate values",
        "Table 4-11 C06 satellite velocity intermediate values",
        Vector3d(-2023579.797096, 35563336.759909, 22879771.741647),
        Vector3d(536.614, -1164.705, 1850.508)
    );

    return 0;
}
