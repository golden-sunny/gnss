/**
 * Exercise 6: ECEF <-> Local ENU coordinate conversion
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include "CoordConvert.h"

Eigen::Vector3d ecefToEnu(const XYZ& refXYZ, const XYZ& targetXYZ) {
    WGS84 wgs84;
    BLH refBLH = xyz2blh(refXYZ, wgs84);

    double lat = refBLH.B();
    double lon = refBLH.L();

    Eigen::Matrix3d R;
    R << -sin(lon),            cos(lon),           0,
         -sin(lat)*cos(lon),  -sin(lat)*sin(lon),  cos(lat),
          cos(lat)*cos(lon),   cos(lat)*sin(lon),  sin(lat);

    Eigen::Vector3d d = targetXYZ - refXYZ;
    return R * d;
}

XYZ enuToEcef(const XYZ& refXYZ, const Eigen::Vector3d& enu) {
    WGS84 wgs84;
    BLH refBLH = xyz2blh(refXYZ, wgs84);

    double lat = refBLH.B();
    double lon = refBLH.L();

    Eigen::Matrix3d R;
    R << -sin(lon),            cos(lon),           0,
         -sin(lat)*cos(lon),  -sin(lat)*sin(lon),  cos(lat),
          cos(lat)*cos(lon),   cos(lat)*sin(lon),  sin(lat);

    Eigen::Vector3d d = R.transpose() * enu;
    return XYZ(refXYZ.x() + d.x(), refXYZ.y() + d.y(), refXYZ.z() + d.z());
}

int main() {
    // Reference station ECEF (example)
    XYZ refXYZ(4081945.67, 2187689.34, 4767321.89);
    // Target point ECEF (example)
    XYZ targetXYZ(4081955.67, 2187699.34, 4767331.89);

    Eigen::Vector3d enu = ecefToEnu(refXYZ, targetXYZ);
    XYZ backXYZ = enuToEcef(refXYZ, enu);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "ENU (E, N, U): " << enu.transpose() << std::endl;
    std::cout << "Back ECEF     : " << backXYZ.transpose() << std::endl;

    return 0;
}
