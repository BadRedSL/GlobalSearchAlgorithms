#pragma once


#include "Support_structures.h"
#include <map>
#include <queue>
#include <chrono>
#include <thread>

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
    const double h = 0.1;
    double numerical_fnc;
    size_t num_thread;

    double m;
    double new_m;
    double M_Max;
    double R;

    std::map<double, characteristics> point_set;
    std::pair<double, characteristics> left_point;
    std::pair<double, characteristics> right_point;
    std::pair<double, characteristics> new_point;
    std::priority_queue < interval, std::vector<interval>, std::greater<interval>>* r_queue;
    std::chrono::steady_clock::time_point start_time;
    result res;

    std::vector<std::thread> thread_vector;
    std::vector<std::pair<double, characteristics>> left_point_vector;
    std::vector<std::pair<double, characteristics>> right_point_vector;
    std::vector<std::pair<double, characteristics>> new_point_vector;

    void compute_first_iteration();
    void clear_solution();
    void recalc_M_Max();
    void recalc_R();
    void update_result();
    void create_new_point();

    void parallel_calculate_new_point(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p, std::pair<double, characteristics> &new_point_p);
    double calculate_new_x_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_x_line_1_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_x_line_2_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_x_roof_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_phi_roof_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_A_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_B_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_d_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_phi_roof_1_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_phi_roof_2_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_phi_roof_3_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p);
    double calculate_M_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p);
    void update_result_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, std::pair<double, characteristics> new_point_p);

    void insert_to_map(double _x, double _z, double _R, double _der_z);
    bool stop_condition();
    void calculate_m();
    double calculate_M();    
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
    double calculate_der_fnc(double _x);

public:
    AGMND(const double _a,
        const double _b,
        const double _r,
        const double _eps,
        const  size_t _max_iter_count,
        double (*_fnc)(double),
        double (*_fnc_der)(double),
        size_t _num_thread = 1,
        bool _numerical_fnc = false);

    void reset(const double _a,
        const double _b,
        const double _r,
        const double _eps,
        const  size_t _maxIterCount,
        double (*_fnc)(double),
        double (*_fnc_der)(double),
        size_t _num_thread = 1,
        bool _numerical_fnc = false);

    result find_min_value();
    result find_min_value_parallel();

    result get_result();
};

