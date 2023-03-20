#pragma once


#include "Support_structures.h"
#include <map>
#include <queue>
#include <chrono>

class AGMND
{
private:
    double a;
    double b;
    double r;
    double eps;
    size_t max_iter_count;
    double (*fnc)(double);
    double (*fnc_der)(double);

    double m;
    double new_m;
    double M_Max;
    double R;

    std::map<double, characteristics> point_set;
    std::pair<double, characteristics> left_point;
    std::pair<double, characteristics> right_point;
    std::pair<double, characteristics> new_point;
    std::priority_queue < interval, std::vector<interval>, std::greater<interval>>* r_queue;
    std::chrono::system_clock::time_point start_time;
    result res;

    void compute_first_iteration();
    void clear_solution();
    void recalc_M_Max();
    void recalc_R();
    void update_result();
    void create_new_point();

    void insert_to_map(double _x, double _z, double _R, double _der_z);
    double calculate_M();
    void calculate_m();
    bool stop_condition();
    double calculate_x_line_1();
    double calculate_x_line_2();
    double calculate_x_roof();
    double calculate_phi_roof(double _x);
    double calculate_A();
    double calculate_B();
    double calculate_d();
    double calculate_phi_roof_1(double _x);
    double calculate_phi_roof_2(double _x);
    double calculate_phi_roof_3(double _x);
    double calculate_R();
    double calculate_new_x();

public:
    AGMND(const double _a,
        const double _b,
        const double _r,
        const double _eps,
        const  size_t _max_iter_count,
        double (*_fnc)(double),
        double (*_fnc_der)(double));

    void reset(const double _a,
        const double _b,
        const double _r,
        const double _eps,
        const  size_t _maxIterCount,
        double (*_fnc)(double),
        double (*_fnc_der)(double));

    result find_min_value();

    result get_result();
};

