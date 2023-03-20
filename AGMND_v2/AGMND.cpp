#include "AGMND.h"

void AGMND::compute_first_iteration()
{
    start_time = std::chrono::system_clock::now();

    res.x = fnc(b) < fnc(a) ? b : a;
    res.z = fnc(b) < fnc(a) ? fnc(b) : fnc(a);
    res.accuracy = b - a;
    res.num_iter = 1;
    res.total_runtime = (std::chrono::system_clock::now() - start_time).count();

    if (a != b)
    {
        left_point.first = a;
        left_point.second = characteristics(fnc(a), fnc_der(a), 0);
        right_point.first = b;
        right_point.second = characteristics(fnc(b), fnc_der(b), 0);

        M_Max = calculate_M();
        calculate_m();

        R = calculate_R();
        left_point.second.R = R;
        right_point.second.R = R;

        insert_to_map(a, fnc(a), R, fnc_der(a));
        insert_to_map(b, fnc(b), R, fnc_der(b));

        r_queue->push(interval({ left_point.first, left_point.second },
            { right_point.first, right_point.second }));
    }
}

void AGMND::clear_solution()
{
    delete r_queue;
    r_queue = nullptr;
    point_set.clear();

    m = 0;
    new_m = 0;
    M_Max = 0;
    R = 0;
}

void AGMND::recalc_M_Max()
{
    M_Max = 0;

    for (auto cur = (++point_set.begin()), prev = point_set.begin(); cur != point_set.end(); ++cur, ++prev)
    {
        left_point.first = prev->first;
        left_point.second = prev->second;
        right_point.first = cur->first;
        right_point.second = cur->second;

        M_Max = calculate_M() > M_Max ? calculate_M() : M_Max;
    }
}

void AGMND::recalc_R()
{
    delete r_queue;
    r_queue = new std::priority_queue < interval, std::vector<interval>, std::greater<interval>>;

    for (auto cur = (++point_set.begin()), prev = point_set.begin(); cur != point_set.end(); ++cur, ++prev)
    {
        left_point.first = prev->first;
        left_point.second = prev->second;
        right_point.first = cur->first;
        right_point.second = cur->second;

        R = calculate_R();
        left_point.second.R = R;
        right_point.second.R = R;

        r_queue->push(interval({ left_point.first, left_point.second },
            { right_point.first, right_point.second }));
    }
}

void AGMND::update_result()
{
    if (new_point.second.z < res.z)
    {
        res.x = new_point.first;
        res.z = new_point.second.z;
    }
    res.num_iter++;
    res.accuracy = right_point.first - left_point.first;
}

void AGMND::create_new_point()
{
    new_point.first = calculate_new_x();
    new_point.second = characteristics(fnc(new_point.first), fnc_der(new_point.first), 0);
}

void AGMND::insert_to_map(double _x, double _z, double _R, double _der_z)
{
    characteristics ch(_z, _R, _der_z);
    point_set.insert({ _x, ch });
    if (_z < res.z) {
        res.x = _x;
        res.z = _z;
    }
}

double AGMND::calculate_M()
{
    double M1, M2, M3, dx, dy;

    dx = abs(right_point.first - left_point.first);
    dy = abs(right_point.second.z - left_point.second.z);

    M1 = abs(right_point.second.der_z - left_point.second.der_z) / dx;
    M2 = 2 * abs((-dy + left_point.second.der_z * dx)) / pow(dx, 2);
    M3 = 2 * abs((dy - right_point.second.der_z * dx)) / pow(dx, 2);

    return std::max({ M1, M2, M3 });
}

void AGMND::calculate_m()
{
    m = (M_Max > 1.0e-8) ? r * M_Max : 1.0;
}

bool AGMND::stop_condition()
{

    return (res.accuracy <= eps) || (res.num_iter >= max_iter_count);
}

double AGMND::calculate_x_line_1()
{
    double numerator = ((left_point.second.z - left_point.second.der_z * left_point.first) -
        (right_point.second.z - right_point.second.der_z * right_point.first)) +
        (m / 2) * (pow(right_point.first, 2) - pow(left_point.first, 2))
        - m * pow(calculate_d(), 2);

    double denominator = m * (right_point.first - left_point.first) +
        (right_point.second.der_z - left_point.second.der_z);

    return numerator / denominator;
}

double AGMND::calculate_x_line_2()
{
    double numerator = ((left_point.second.z - left_point.second.der_z * left_point.first) -
        (right_point.second.z - right_point.second.der_z * right_point.first)) +
        (m / 2) * (pow(right_point.first, 2) - pow(left_point.first, 2))
        + m * pow(calculate_d(), 2);

    double denominator = m * (right_point.first - left_point.first) +
        (right_point.second.der_z - left_point.second.der_z);

    return numerator / denominator;
}

double AGMND::calculate_x_roof()
{
    return (-left_point.second.der_z + m * (calculate_x_line_1() - left_point.first) +
        m * right_point.first) / m;
}

double AGMND::calculate_phi_roof(double _x)
{
    if ((left_point.first <= _x) && (_x <= calculate_x_line_1()))
    {
        return calculate_phi_roof_1(_x);
    }
    else if ((calculate_x_line_1() <= _x) && (_x <= calculate_x_line_2()))
    {
        return calculate_phi_roof_2(_x);
    }
    else if ((calculate_x_line_2() <= _x) && (_x <= right_point.first))
    {
        return calculate_phi_roof_3(_x);
    }
}

double AGMND::calculate_A()
{
    return left_point.second.der_z - m * (calculate_x_line_1() - left_point.first);
}

double AGMND::calculate_B()
{
    return calculate_phi_roof_1(calculate_x_line_1());
}

double AGMND::calculate_d()
{
    return ((right_point.first - left_point.first) / 2) -
        ((right_point.second.der_z - left_point.second.der_z) / (2 * m));
}

double AGMND::calculate_phi_roof_1(double _x)
{
    return left_point.second.z + left_point.second.der_z *
        (_x - left_point.first) - 0.5 * m * pow(_x - left_point.first, 2);
}

double AGMND::calculate_phi_roof_2(double _x)
{
    return calculate_A() * (_x - calculate_x_line_1()) +
        0.5 * m * pow(_x - calculate_x_line_1(), 2) + calculate_B();
}

double AGMND::calculate_phi_roof_3(double _x)
{
    return right_point.second.z + right_point.second.der_z *
        (right_point.first - _x) - 0.5 * m * pow(right_point.first - _x, 2);
}

double AGMND::calculate_R()
{
    if ((calculate_x_line_1() <= calculate_x_roof()) && (calculate_x_roof() <= calculate_x_line_2()))
    {
        return calculate_phi_roof(calculate_x_roof());
    }
    else
    {
        return std::min(calculate_phi_roof(calculate_x_line_1()), calculate_phi_roof(calculate_x_line_2()));
    }
}

double AGMND::calculate_new_x()
{
    left_point = r_queue->top().left_point;
    right_point = r_queue->top().right_point;

    if ((calculate_x_line_1() <= calculate_x_roof()) && (calculate_x_roof() <= calculate_x_line_2()))
    {
        return calculate_x_roof();
    }
    else if (calculate_phi_roof(calculate_x_line_1()) <= calculate_phi_roof(calculate_x_line_2()))
    {
        return calculate_x_line_1();
    }
    else
    {
        return calculate_x_line_2();
    }

}

AGMND::AGMND(const double _a, const double _b, const double _r, const double _eps, const size_t _max_iter_count,
    double(*_fnc)(double), double(*_fnc_der)(double)) : a(_a), b(_b), r(_r), eps(_eps), max_iter_count(_max_iter_count), fnc(_fnc), fnc_der(_fnc_der)
{
    M_Max = 0;
    m = 0;
    new_m = 0;
    R = 0;
    r_queue = new std::priority_queue < interval, std::vector<interval>, std::greater<interval>>;
}

void AGMND::reset(const double _a, const double _b, const double _r, const double _eps,
    const size_t _maxIterCount, double(*_fnc)(double), double(*_fnc_der)(double)) {
    a = _a;
    b = _b;
    r = _r;
    eps = _eps;
    max_iter_count = _maxIterCount;
    fnc = _fnc;
    fnc_der = _fnc_der;
    if (r_queue != nullptr)
    {
        delete r_queue;
    }
    r_queue = new std::priority_queue < interval, std::vector<interval>, std::greater<interval>>;
}

result AGMND::find_min_value()
{
    compute_first_iteration();

    while (!stop_condition())
    {
        create_new_point();
        insert_to_map(new_point.first, fnc(new_point.first), 0, fnc_der(new_point.first));
        update_result();
        recalc_M_Max();
        calculate_m();
        recalc_R();
    }

    clear_solution();
    res.total_runtime = (std::chrono::system_clock::now() - start_time).count();
    return res;
}

result AGMND::get_result()
{
    return res;
}
