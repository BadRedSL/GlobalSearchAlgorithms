#include "Support_structures.h"

result::result(): x(0.0), z(0.0), accuracy(0.0), total_runtime(0.0), num_iter(0) {}

result::result(double _x, double _z, double _accuracy, double _total_runtime,
    size_t _num_iter) : x(_x), z(_z), accuracy(_accuracy),
    total_runtime(_total_runtime), num_iter(_num_iter) {}

characteristics::characteristics(): z(0.0), der_z(0.0), R(0.0) {}

characteristics::characteristics(double _z, double _der_z, double _R) : z(_z), der_z(_der_z), R(_R) {}


interval::interval() {}

interval::interval(std::pair<double, characteristics> _left_point,
    std::pair<double, characteristics> _right_point) : left_point(_left_point),
    right_point(_right_point) {}

bool operator<(const interval& a, const interval& b)
{
    return a.left_point.second.R < b.left_point.second.R;
}

bool operator>(const interval& a, const interval& b)
{
    return a.left_point.second.R > b.left_point.second.R;
}
