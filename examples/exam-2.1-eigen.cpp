//
// Created by shjzh on 2024/12/24.
//

/**
 * Copyright:
 *  This software is licensed under the Mulan Permissive Software License, Version 2 (MulanPSL-2.0).
 *  You may obtain a copy of the License at:http://license.coscl.org.cn/MulanPSL2
 *  As stipulated by the MulanPSL-2.0, you are granted the following freedoms:
 *      To copy, use, and modify the software;
 *      To use the software for commercial purposes;
 *      To redistribute the software.
 *
 * Author: shoujian zhang
 *
 * References:
 * 1. Sanz Subirana, J., Juan Zornoza, J. M., & Hernandez-Pajares, M. (2013).
 *    GNSS data processing: Volume I: Fundamentals and algorithms. ESA Communications.
 * 2. Eckel, Bruce. Thinking in C++. 2nd ed., Prentice Hall, 2000.
 */

#include <iostream>
#include <Eigen/Dense>

int main() {
    using Eigen::MatrixXd;

    MatrixXd A(2, 2);
    A << 3, -1,
         2.5, 1.5;

    MatrixXd B(2, 2);
    B << 1, 4,
         -2, 0.5;

    MatrixXd add = A + B;
    MatrixXd sub = A - B;
    MatrixXd mul = A * B;
    MatrixXd invA = A.inverse();

    std::cout << "A:\n" << A << "\n\n";
    std::cout << "B:\n" << B << "\n\n";
    std::cout << "A + B:\n" << add << "\n\n";
    std::cout << "A - B:\n" << sub << "\n\n";
    std::cout << "A * B:\n" << mul << "\n\n";
    std::cout << "A^{-1}:\n" << invA << "\n\n";
    std::cout << "A * A^{-1} (check):\n" << (A * invA) << std::endl;

    return 0;
}
