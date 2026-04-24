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
#include <cmath>
#include <stdexcept>
#include <memory>
#include <iomanip>

using namespace std;

// 参考框架基类
class ReferenceFrame {
public:
    virtual ~ReferenceFrame() = default;

    // 获取长半轴
    virtual double getA() const = 0;

    // 获取扁率
    virtual double getF() const = 0;

    // 获取地球自转角速度
    virtual double getOmega() const = 0;

    // 获取地球引力常数
    virtual double getGM() const = 0;

    // 获取二阶田谐系数
    virtual double getJ2() const = 0;

    // 计算第一偏心率平方 e^2
    double getE2() const {
        double f = getF();
        return 2 * f - f * f;
    }
};

// WGS84 参考框架
class WGS84 : public ReferenceFrame {
public:
    double getA() const override {
        return 6378137.0; // 长半轴 (米)
    }

    double getF() const override {
        return 1 / 298.257223563; // 扁率
    }

    double getOmega() const override {
        return 7.292115e-5; // 地球自转角速度 (弧度/秒)
    }

    double getGM() const override {
        return 3.986004418e14; // 地球引力常数 (米^3/秒^2)
    }

    double getJ2() const override {
        return 1.08263e-3; // 二阶田谐系数
    }
};

// PZ90 参考框架
class PZ90 : public ReferenceFrame {
public:
    double getA() const override {
        return 6378136.0; // 长半轴 (米)
    }

    double getF() const override {
        return 1 / 298.257839303; // 扁率
    }

    double getOmega() const override {
        return 7.2921150e-5; // 地球自转角速度 (弧度/秒)
    }

    double getGM() const override {
        return 3.9860044e14; // 地球引力常数 (米^3/秒^2)
    }

    double getJ2() const override {
        return 1.08262575e-3; // 二阶田谐系数
    }
};

// XYZ 坐标类
class XYZ {
public:
    double x, y, z;

    XYZ(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

// BLH 坐标类
class BLH {
public:
    double B, L, H;

    BLH(double B_, double L_, double H_) : B(B_), L(L_), H(H_) {}
};

// 坐标转换函数
BLH xyz2blh(const XYZ &xyz, const ReferenceFrame &frame) {
    // 获取椭球参数
    double a = frame.getA();
    double e2 = frame.getE2();

    // 计算水平距离 rho（即 sqrt(x^2 + y^2)）
    double rho = sqrt(xyz.x * xyz.x + xyz.y * xyz.y);

    // 定义阈值，用于判断是否在极点
    const double eps = 1.0e-13;

    // 判断是否在极点
    if (rho < eps) {
        // 在极点，根据 z 的符号判断是南极还是北极
        double B = (xyz.z > 0) ? M_PI / 2 : -M_PI / 2;  // 北极为 +90°，南极为 -90°
        double L = 0.0;  // 经度在极点无定义，通常设为 0
        double H = fabs(xyz.z) - a * sqrt(1 - e2);  // 高度计算

        return BLH(B, L, H);
    }

    // 不在极点，正常计算
    double B0 = atan2(xyz.z, rho);

    // 迭代计算大地纬度 B
    const int maxIterations = 100;
    int iterationCount = 0;
    double B1, N;
    do {
        N = a / sqrt(1 - e2 * sin(B0) * sin(B0));
        B1 = atan2(xyz.z + e2 * N * sin(B0), rho);

        if (fabs(B1 - B0) < eps) break;

        B0 = B1;
        iterationCount++;

        if (iterationCount > maxIterations) {
            throw std::runtime_error("Iteration did not converge.");
        }
    } while (true);

    // 计算大地经度 L
    double L = atan2(xyz.y, xyz.x);

    // 计算高度 H
    double H = rho / cos(B1) - N;

    // 返回大地坐标
    return BLH(B1, L, H);
}

int main() {
    try {
        // 示例 XYZ 坐标
        XYZ xyz_north_pole(0.0, 0.0, 6356752.314);  // 北极点
        XYZ xyz_south_pole(0.0, 0.0, -6356752.314); // 南极点
        XYZ xyz_normal(4081945.67, 2187689.34, 4767321.89); // 正常点

        // 使用动态绑定的参考框架
        WGS84 wgs84;
        PZ90 pz90;

        // 使用 PZ-90 参考框架进行坐标转换
        BLH blh_north_pole = xyz2blh(xyz_north_pole, pz90);
        std::cout << fixed << std::setprecision(14) << endl;
        std::cout << "PZ-90  north_pole (B, L, H): "
                  << blh_north_pole.B << ", " << blh_north_pole.L << ", " << blh_north_pole.H << std::endl;

        BLH blh_south_pole = xyz2blh(xyz_south_pole, pz90);
        std::cout << "PZ-90  south_pole (B, L, H): "
                  << blh_south_pole.B << ", " << blh_south_pole.L << ", " << blh_south_pole.H << std::endl;

        BLH blh_normal = xyz2blh(xyz_normal, pz90);
        std::cout << "PZ-90  normal_point  (B, L, H): "
                  << blh_normal.B << ", " << blh_normal.L << ", " << blh_normal.H << std::endl;

        // 使用 WGS84 参考框架进行坐标转换
        BLH blh_normal_wgs84 = xyz2blh(xyz_normal, wgs84);
        std::cout << "WGS84  normal_point  (B, L, H): "
                  << blh_normal_wgs84.B << ", " << blh_normal_wgs84.L << ", " << blh_normal_wgs84.H << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}


