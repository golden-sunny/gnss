//
// Created by shjzh on 2025/3/31.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

struct Satellite {
    double azimuth;    // 方位角（度）
    double elevation;  // 仰角（度）
    std::string id;    // 卫星编号
    std::string color; // 显示颜色
};

void polarToCartesian(double centerX, double centerY, double radius,
                      double azimuthDeg, double elevationDeg,
                      double& x, double& y) {
    const double elevMax = 90.0;
    double r = radius * (elevMax - elevationDeg) / elevMax;
    double azimuthRad = azimuthDeg * M_PI / 180.0;
    x = centerX + r * sin(azimuthRad);
    y = centerY - r * cos(azimuthRad);  // SVG坐标系Y轴向下
}

void generateSkyplot(const std::vector<Satellite>& satellites,
                     const std::string& filename = "skyplot.svg") {
    const double canvasSize = 500.0;
    const double center = canvasSize / 2;
    const double plotRadius = 200.0;

    std::ofstream svg(filename);
    if (!svg) {
        std::cerr << "无法创建文件: " << filename << std::endl;
        return;
    }

    // SVG头
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg width=\"" << canvasSize << "\" height=\"" << canvasSize
        << "\" viewBox=\"0 0 " << canvasSize << " " << canvasSize << "\" "
        << "xmlns=\"http://www.w3.org/2000/svg\">\n";

    // 绘制背景
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // 绘制仰角同心圆（15度间隔）
    for (int elev = 0; elev <= 90; elev += 15) {
        double r = plotRadius * (90.0 - elev) / 90.0;
        svg << "<circle cx=\"" << center << "\" cy=\"" << center
            << "\" r=\"" << r << "\" fill=\"none\" stroke=\"#e0e0e0\" "
            << "stroke-width=\"1\"/>\n";
    }

    // 绘制方位角射线（30度间隔）
    for (int az = 0; az < 360; az += 30) {
        double rad = az * M_PI / 180.0;
        double x = center + plotRadius * sin(rad);
        double y = center - plotRadius * cos(rad);
        svg << "<line x1=\"" << center << "\" y1=\"" << center
            << "\" x2=\"" << x << "\" y2=\"" << y << "\" "
            << "stroke=\"#e0e0e0\" stroke-width=\"1\"/>\n";
    }

    // 绘制方位角标签
    svg << "<text x=\"" << center << "\" y=\"30\" text-anchor=\"middle\">N</text>\n"
        << "<text x=\"" << (canvasSize-20) << "\" y=\"" << center
        << "\" text-anchor=\"middle\">E</text>\n"
        << "<text x=\"" << center << "\" y=\"" << (canvasSize-20)
        << "\" text-anchor=\"middle\">S</text>\n"
        << "<text x=\"20\" y=\"" << center << "\" text-anchor=\"middle\">W</text>\n";

    // 绘制卫星
    for (const auto& sat : satellites) {
        double x, y;
        polarToCartesian(center, center, plotRadius,
                         sat.azimuth, sat.elevation, x, y);

        // 卫星符号
        svg << "<circle cx=\"" << x << "\" cy=\"" << y
            << "\" r=\"6\" fill=\"" << sat.color << "\" "
            << "stroke=\"black\" stroke-width=\"0.5\"/>\n";

        // 卫星编号
        svg << "<text x=\"" << x + 8 << "\" y=\"" << y + 4
            << "\" font-size=\"10\" fill=\"black\">"
            << sat.id << "</text>\n";
    }

    svg << "</svg>";
    svg.close();
}

int main() {
    // 示例卫星数据
    std::vector<Satellite> satellites = {
            {45.0, 30.0, "G12", "#FF4444"},   // 红色：GPS
            {135.0, 60.0, "R05", "#44FF44"},   // 绿色：GLONASS
            {270.0, 45.0, "E21", "#4444FF"},   // 蓝色：Galileo
            {315.0, 15.0, "C07", "#FFAA00"}    // 橙色：北斗
    };

    string filename = "D:\\documents\\Source\\gnssLab-2.2\\data\\skyplot.svg";

    generateSkyplot(satellites, filename);

    std::cout << "天空图已生成到 skyplot.svg" << std::endl;
    return 0;
}