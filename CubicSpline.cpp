//
// Created by vadim on 13.04.2021.
//

#include "CubicSpline.h"
#include <iostream>
#include <iomanip>

void print(std::vector<double> &v) {
    for (auto i : v)
        std::cout << std::fixed << std::setprecision(5) << i << ' ';
    std::cout << std::endl;
}

static std::vector<double> _tridiagonal_matrix_solution(const std::vector<double> &a, const std::vector<double> &b,
                                                        const std::vector<double> &c, const std::vector<double> &d,
                                                        int n) {
    std::vector<double> solution;
    std::vector<double> p = {-c.at(0) / b.at(0)},
            q = {d.at(0) / b.at(0)};
    for (int i = 1; i < n - 1; ++i) {
        p.push_back(-c.at(i) / (a.at(i) * p.at(i - 1) + b.at(i)));
        q.push_back((d.at(i) - a.at(i) * q.at(i - 1)) / (b.at(i) + p.at(i - 1) * a.at(i)));
    }
//    std::cout << "p: ";
//    print(p);
//    std::cout << "q: ";
//    print(q);
    solution.resize(n + 2);
    solution[n - 1] =
            (d.at(n - 1) - a.at(n - 1) * q.at(n - 2)) / (b.at(n - 1) + a.at(n - 1) * p.at(n - 2));
    for (int i = n - 2; i >= 0; --i)
        solution.at(i) = p.at(i) * solution.at(i + 1) + q.at(i);
    return solution;
}

void CubicSpline::Init(std::vector<double> &x, std::vector<double> &f_x,
                       double alpha0, double alpha1, double beta0, double beta1,
                       double (*fd1)(double), double (*fd2)(double)) {
    n = x.size();
    _x.resize(n + 2);
    std::copy(x.begin(), x.end(), _x.begin());
    _f_x.resize(n + 2);
    std::copy(f_x.begin(), f_x.end(), _f_x.begin());

//    std::cout << "_f_x: ";
//    print(_f_x);
//    std::cout << "_x: ";
//    print(_x);


    double gamma0 = alpha0 * fd1(x[0]) + beta0 * fd2(x[0]),
            gamma1 = alpha1 * fd1(x.at(n - 1)) + beta1 * fd2(x.at(n - 1));

    std::vector<double> h;
    h.resize(n + 2);
    for (int i = 1; i <= n - 1; ++i)
        h.at(i) = _x[i] - _x[i - 1];

    std::vector<double> a, b, c, d;
    a.resize(n + 2);
    b.resize(n + 2);
    c.resize(n + 2);
    d.resize(n + 2);

    for (int i = 1; i <= n - 2; ++i) {
        a.at(i) = h.at(i) / 6;
        b.at(i) = (h.at(i) + h.at(i + 1)) / 3;
        c.at(i) = h.at(i + 1) / 6;
        d.at(i) = (_f_x.at(i + 1) - _f_x.at(i)) / h.at(i + 1) - (_f_x.at(i) - _f_x.at(i - 1)) / h.at(i);
    }
    c[0] = -alpha0 * h.at(1) / 6;
    b[0] = -alpha0 * h.at(1) / 3 + beta0;
    b[n - 1] = alpha1 * h[n - 1] / 3 + beta1;
    a[n - 1] = alpha1 * h[n - 1] / 6;
    d[0] = gamma0 - alpha0 * (_f_x[1] - _f_x[0]) / h[1];
    d[n - 1] = gamma1 - alpha1 * (_f_x[n - 1] - _f_x[n - 2]) / h[n - 1];
    a[0] = 0;
    c[n - 1] = 0;

//    std::cout << "h: ";
//    print(h);
//
//    std::cout << "a: ";
//    print(a);
//    std::cout << "b: ";
//    print(b);
//    std::cout << "c: ";
//    print(c);
//    std::cout << "d: ";
//    print(d);

    _m = _tridiagonal_matrix_solution(a, b, c, d, n);

//    std::cout << "_m: ";
//    print(_m);
}

double CubicSpline::CalcInPoint(double p) {
    double res = 0;
    for (int i = 1; i <= n - 1; ++i) {
        if (_x[i - 1] <= p && p < _x[i]) {
            double h = _x[i] - _x[i - 1];
            res = _m[i - 1] * ((_x[i] - p) * (_x[i] - p) * (_x[i] - p)) / (6 * h) +
                  _m[i] * ((p - _x[i - 1]) * (p - _x[i - 1]) * (p - _x[i - 1])) / (6 * h) +
                  (_f_x[i - 1] - (_m[i - 1] * h * h) / 6) * ((_x[i] - p) / h) + ((_f_x[i] -
                                                                                  (_m[i] * h * h) / 6) *
                                                                                 ((p - _x[i - 1]) / h));
            break;
        }
    }
    return res;
}

std::valarray<double> CubicSpline::CalcInPoints(const std::valarray<double> &points) {
    std::valarray<double> retv(points.size());
    for (int i = 0; i < points.size(); ++i)
        retv[i] = CalcInPoint(points[i]);
    return retv;
}
