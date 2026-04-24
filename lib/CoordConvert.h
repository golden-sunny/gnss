/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: Shoujian Zhang，shjzhang@sgg.whu.edu.cn， 2024-10-10
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernández-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */

#ifndef GNSSLAB_COORDCONVERT_H
#define GNSSLAB_COORDCONVERT_H

#include "CoordStruct.h"
#include "Exception.h"
#include <Eigen/Eigen>

#define debug 1

// 坐标转换函数
inline BLH xyz2blh(const XYZ &xyz, const ReferenceFrame &frame) {
    // 获取椭球参数
    double a = frame.getA();
    double e2 = frame.getE2();
    // 计算水平距离 rho（即 sqrt(x^2 + y^2)）
    double rho = sqrt(xyz.X() * xyz.X() + xyz.Y() * xyz.Y());

    // 定义阈值，用于判断是否在极点
    const double eps = 1.0e-13;

    // 判断是否在极点
    if (rho < eps) {
        // 在极点，根据 z 的符号判断是南极还是北极
        double B = (xyz.Z() > 0) ? M_PI / 2 : -M_PI / 2;  // 北极为 +90°，南极为 -90°
        double L = 0.0;  // 经度在极点无定义，通常设为 0
        double H = fabs(xyz.Z()) - a * sqrt(1 - e2);  // 高度计算

        return BLH(B, L, H);
    }

    // 不在极点，正常计算
    double B0 = atan2(xyz.Z(), rho);

    // 迭代计算大地纬度 B
    const int maxIterations = 100;
    int iterationCount = 0;
    double B1, N;
    do {
        N = a / sqrt(1 - e2 * sin(B0) * sin(B0));
        B1 = atan2(xyz.Z() + e2 * N * sin(B0), rho);

        if (fabs(B1 - B0) < eps) break;

        B0 = B1;
        iterationCount++;

        if (iterationCount > maxIterations) {
            throw std::runtime_error("Iteration did not converge.");
        }
    } while (true);

    // 计算大地经度 L
    double L = atan2(xyz.Y(), xyz.X());

    // 计算高度 H
    double H = rho / cos(B1) - N;

    // 返回大地坐标
    return BLH(B1, L, H);
}

// computes the elevation of the input (Target) position as seen from ref Position, using a Geodetic
// (i.e. ellipsoidal) system.
// @return the elevation in degrees
inline double elevation(const XYZ& refXYZ, const XYZ& targetXYZ)
noexcept(false)
{
    BLH refBLH;
    WGS84 wgs84;
    refBLH = xyz2blh(refXYZ, wgs84);

    if(debug)
        cout << "refBLH:" << refBLH.transpose() << endl;

    double lat = refBLH.B();
    double lon = refBLH.L();

    double localUp;
    double cosUp;

    Eigen::Vector3d z;

    // Let's get the slant vector, 这里需要修改接口
    z = targetXYZ - refXYZ;

    if (z.norm()<=1e-4) // if the positions are within .1 millimeter
    {
        InvalidRequest e("Positions are within .1 millimeter");
        throw(e);
    }

    // Compute k vector in local North-East-Up (NEU) system
    Eigen::Vector3d kVector(::cos(lat)*::cos(lon), ::cos(lat)*::sin(lon), ::sin(lat));

    // Take advantage of dot method to get Up coordinate in local NEU system
    localUp = z.dot(kVector);

    // Let's get cos(z), being z the angle with respect to local vertical (Up);
    cosUp = localUp/z.norm();

    if(debug)
        cout << "cosUp:" << cosUp << endl;

    double elev = 90.0 - ((::acos(cosUp))*RAD_TO_DEG);
    if(debug)
        cout << "elev:"<<  elev << endl;

    return elev;
}

// A member function that computes the azimuth of the input
// (Target) position as seen from this Position, using a Geodetic
// (i.e. ellipsoidal) system.
// @param Target the Position which is observed to have the
//        computed azimuth, as seen from this Position.
// @return the azimuth in degrees
inline double azimuth(const XYZ& refXYZ, const XYZ& targetXYZ)
noexcept(false)
{
    WGS84 wgs84;
    BLH refBLH = xyz2blh(refXYZ,wgs84);;

    double latRad = refBLH.B();
    double lonRad = refBLH.L();

    double localN, localE;

    Eigen::Vector3d z;
    // Let's get the slant vector
    z = targetXYZ - refXYZ;

    if (z.norm()<=1e-4) // if the positions are within .1 millimeter
    {
        GeometryException ge("azimuthGeodetic::Positions are within .1 millimeter");
        throw(ge);
    }

    // Compute i vector in local North-East-Up (NEU) system
    Eigen::Vector3d iVector(-::sin(latRad)*::cos(lonRad),
                            -::sin(latRad)*::sin(lonRad),
                            ::cos(latRad));

    // Compute j vector in local North-East-Up (NEU) system
    Eigen::Vector3d jVector(-::sin(lonRad),
                            ::cos(lonRad),
                            0);

    // Now, let's use dot product to get localN and localE unitary vectors
    localN = (z.dot(iVector))/z.norm();
    localE = (z.dot(jVector))/z.norm();

    // Let's test if computing azimuth has any sense
    double test = fabs(localN) + fabs(localE);

    // Warning: If elevation is very close to 90 degrees, we will return azimuth = 0.0
    if (test < 1.0e-16) return 0.0;

    double alpha = ((::atan2(localE, localN)) * RAD_TO_DEG);
    if (alpha < 0.0)
    {
        return alpha + 360.0;
    }
    else
    {
        return alpha;
    }
}

#endif //GNSSLAB_COORDCONVERT_H
