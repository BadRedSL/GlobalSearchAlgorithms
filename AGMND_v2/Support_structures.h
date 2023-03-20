#pragma once


#include <algorithm>


struct result {
    double x;
    double z;
    double accuracy;
    double total_runtime;
    size_t num_iter;
    result();
    result(double _x,
        double _z,
        double _accuracy,
        double _total_runtime,
        size_t _num_iter);
};

struct characteristics {
    double z;
    double der_z;
    double R;

    characteristics();
    characteristics(double _z,
        double _der_z,
        double _R);
};

struct interval {
    std::pair<double, characteristics> left_point;
    std::pair<double, characteristics> right_point;
    interval();
    interval(std::pair<double, characteristics> _left_point,
        std::pair<double, characteristics> _right_point);

    friend bool operator< (const interval& a, const interval& b);
    friend bool operator> (const interval& a, const interval& b);
};