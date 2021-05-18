//
// Created by vadim on 13.04.2021.
//

#ifndef SPLINETEST_CUBICSPLINE_H
#define SPLINETEST_CUBICSPLINE_H

#include <cmath>
#include <vector>
#include <array>
#include <valarray>
#include <utility>

class CubicSpline {
    std::vector<double> _x;
    std::vector<double> _f_x;

    int n;

    std::vector<double> _m;

public:
    CubicSpline() = default;

    void Init(std::vector<double> &x, std::vector<double> &f_x,
              double alpha0, double alpha1, double beta0, double beta1,
              double (*fd1)(double), double (*fd2)(double));

    double CalcInPoint(double p);

    std::valarray<double> CalcInPoints(const std::valarray<double> &points);
};


#endif //SPLINETEST_CUBICSPLINE_H
