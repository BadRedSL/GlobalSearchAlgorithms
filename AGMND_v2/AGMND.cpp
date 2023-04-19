#include "AGMND.h"

void AGMND::compute_first_iteration()
{
    start_time = std::chrono::steady_clock::now();

    res.x = fnc(b) < fnc(a) ? b : a;
    res.z = fnc(b) < fnc(a) ? fnc(b) : fnc(a);
    res.accuracy = b - a;
    res.num_iter = 1;
    res.total_runtime = (std::chrono::steady_clock::now() - start_time).count();

    if (a != b)
    {
        left_point.first = a;
        left_point.second = (numerical_fnc)? characteristics(fnc(a), calculate_der_fnc(a), 0) : characteristics(fnc(a), fnc_der(a), 0);
        right_point.first = b;
        right_point.second = (numerical_fnc) ? characteristics(fnc(b), calculate_der_fnc(b), 0) : characteristics(fnc(b), fnc_der(b), 0);

        M_Max = calculate_M();
        calculate_m();

        R = calculate_R();
        left_point.second.R = R;
        right_point.second.R = R;

        insert_to_map(a, left_point.second.z, R, left_point.second.der_z);
        insert_to_map(b, right_point.second.z, R, right_point.second.der_z);

        r_queue->push(interval({ left_point.first, left_point.second },
            { right_point.first, right_point.second }));
    }
}


void AGMND::clear_solution()
{
    delete r_queue;
    r_queue = nullptr;
    point_set.clear();
    left_point_vector.clear();
    right_point_vector.clear();
    new_point_vector.clear();

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
    new_point.second = (numerical_fnc)? characteristics(fnc(new_point.first), calculate_der_fnc(new_point.first), 0) : characteristics(fnc(new_point.first), fnc_der(new_point.first), 0);
}

void AGMND::parallel_calculate_new_point(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p, std::pair<double, characteristics>& new_point_p)
{
    new_point_p.first = calculate_new_x_p(left_point_p, right_point_p, m_p);
    new_point_p.second = (numerical_fnc) ? characteristics(fnc(new_point_p.first), calculate_der_fnc(new_point_p.first), 0) : characteristics(fnc(new_point_p.first), fnc_der(new_point_p.first), 0);
}

double AGMND::calculate_new_x_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    if ((calculate_x_line_1_p(left_point_p, right_point_p, m_p) <= calculate_x_roof_p(left_point_p, right_point_p, m_p))
        && (calculate_x_roof_p(left_point_p, right_point_p, m_p) <= calculate_x_line_2_p(left_point_p, right_point_p, m_p)))
    {
        return calculate_x_roof_p(left_point_p, right_point_p, m_p);
    }
    else if (calculate_phi_roof_p(calculate_x_line_1_p(left_point_p, right_point_p, m_p), left_point_p, right_point_p, m_p)
        <= calculate_phi_roof_p(calculate_x_line_2_p(left_point_p, right_point_p, m_p), left_point_p, right_point_p, m_p))
    {
        return calculate_x_line_1_p(left_point_p, right_point_p, m_p);
    }
    else
    {
        return calculate_x_line_2_p(left_point_p, right_point_p, m_p);
    }
}

double AGMND::calculate_d_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return ((right_point_p.first - left_point_p.first) / 2) -
        ((right_point_p.second.der_z - left_point_p.second.der_z) / (2 * m_p));
}

double AGMND::calculate_x_line_1_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    double numerator = ((left_point_p.second.z - left_point_p.second.der_z * left_point_p.first) -
        (right_point_p.second.z - right_point_p.second.der_z * right_point_p.first)) +
        (m / 2) * (pow(right_point_p.first, 2) - pow(left_point_p.first, 2))
        - m * pow(calculate_d_p(left_point_p, right_point_p, m_p), 2);

    double denominator = m * (right_point_p.first - left_point_p.first) +
        (right_point_p.second.der_z - left_point_p.second.der_z);

    return numerator / denominator;
}

double AGMND::calculate_x_line_2_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    double numerator = ((left_point_p.second.z - left_point_p.second.der_z * left_point_p.first) -
        (right_point_p.second.z - right_point_p.second.der_z * right_point_p.first)) +
        (m / 2) * (pow(right_point_p.first, 2) - pow(left_point_p.first, 2))
        + m * pow(calculate_d_p(left_point_p, right_point_p, m_p), 2);

    double denominator = m * (right_point_p.first - left_point_p.first) +
        (right_point_p.second.der_z - left_point_p.second.der_z);

    return numerator / denominator;
}

double AGMND::calculate_x_roof_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return (-left_point_p.second.der_z + m * (calculate_x_line_1_p(left_point_p, right_point_p, m_p) - left_point_p.first) +
        m * right_point_p.first) / m;
}

double AGMND::calculate_phi_roof_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    if ((left_point_p.first <= _x) && (_x <= calculate_x_line_1_p(left_point_p, right_point_p, m_p)))
    {
        return calculate_phi_roof_1_p(_x, left_point_p, right_point_p, m_p);
    }
    else if ((calculate_x_line_1_p(left_point_p, right_point_p, m_p) <= _x) && (_x <= calculate_x_line_2_p(left_point_p, right_point_p, m_p)))
    {
        return calculate_phi_roof_2_p(_x, left_point_p, right_point_p, m_p);
    }
    else if ((calculate_x_line_2_p(left_point_p, right_point_p, m_p) <= _x) && (_x <= right_point_p.first))
    {
        return calculate_phi_roof_3_p(_x, left_point_p, right_point_p, m_p);
    }
}

double AGMND::calculate_phi_roof_1_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return left_point_p.second.z + left_point_p.second.der_z *
        (_x - left_point_p.first) - 0.5 * m * pow(_x - left_point_p.first, 2);
}

double AGMND::calculate_phi_roof_2_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return calculate_A_p(left_point_p, right_point_p, m_p) * (_x - calculate_x_line_1_p(left_point_p, right_point_p, m_p)) +
        0.5 * m * pow(_x - calculate_x_line_1_p(left_point_p, right_point_p, m_p), 2) + calculate_B_p(left_point_p, right_point_p, m_p);
}

double AGMND::calculate_phi_roof_3_p(double _x, std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return right_point_p.second.z + right_point_p.second.der_z *
        (right_point_p.first - _x) - 0.5 * m * pow(right_point_p.first - _x, 2);
}

double AGMND::calculate_A_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return left_point_p.second.der_z - m * (calculate_x_line_1_p(left_point_p, right_point_p, m_p) - left_point_p.first);
}

double AGMND::calculate_B_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, double m_p)
{
    return calculate_phi_roof_1_p(calculate_x_line_1_p(left_point_p, right_point_p, m_p), left_point_p, right_point_p, m_p);
}

double AGMND::calculate_M_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p)
{
    double M1, M2, M3, dx, dy;

    dx = abs(right_point_p.first - left_point_p.first);
    dy = abs(right_point_p.second.z - left_point_p.second.z);

    M1 = abs(right_point_p.second.der_z - left_point_p.second.der_z) / dx;
    M2 = 2 * abs((-dy + left_point_p.second.der_z * dx)) / pow(dx, 2);
    M3 = 2 * abs((dy - right_point_p.second.der_z * dx)) / pow(dx, 2);

    return std::max({ M1, M2, M3 });
}

void AGMND::update_result_p(std::pair<double, characteristics> left_point_p, std::pair<double, characteristics> right_point_p, std::pair<double, characteristics> new_point_p)
{
    if (new_point_p.second.z < res.z)
    {
        res.x = new_point_p.first;
        res.z = new_point_p.second.z;
    }
    res.num_iter++;
    res.accuracy = right_point_p.first - left_point_p.first;
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

double AGMND::calculate_der_fnc(double _x)
{
    for (size_t i = 0; i < 100000; ++i)
    {
        auto tmp = std::chrono::system_clock::now();
    }
    return (fnc(_x + h) - fnc(_x - h)) / (2 * h);
}

AGMND::AGMND(const double _a, const double _b, const double _r, const double _eps, const size_t _max_iter_count,
    double(*_fnc)(double), double(*_fnc_der)(double), size_t _num_thread, bool _numerical_fnc) : a(_a), b(_b), r(_r), eps(_eps),
    max_iter_count(_max_iter_count), fnc(_fnc), fnc_der(_fnc_der), num_thread(_num_thread), numerical_fnc(_numerical_fnc)
{
    M_Max = 0;
    m = 0;
    new_m = 0;
    R = 0;
    r_queue = new std::priority_queue < interval, std::vector<interval>, std::greater<interval>>;
    left_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
    right_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
    new_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
}

void AGMND::reset(const double _a, const double _b, const double _r, const double _eps,
    const size_t _maxIterCount, double(*_fnc)(double), double(*_fnc_der)(double), size_t _num_thread, bool _numerical_fnc) {
    numerical_fnc = numerical_fnc;
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
    num_thread = _num_thread;
    left_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
    right_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
    new_point_vector = std::vector<std::pair<double, characteristics>>(num_thread);
}

result AGMND::find_min_value()
{
    compute_first_iteration();

    while (!stop_condition())
    {
        create_new_point();
        insert_to_map(new_point.first, new_point.second.z, 0, new_point.second.der_z);
        update_result();
        recalc_M_Max();
        calculate_m();
        recalc_R();
    }

    clear_solution();
    res.total_runtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();
    return res;
}

result AGMND::find_min_value_parallel()
{
    compute_first_iteration();

    while (!stop_condition())
    {
        if (num_thread > r_queue->size())
        {
            create_new_point();
            insert_to_map(new_point.first, new_point.second.z, 0, new_point.second.der_z);
            update_result();
        }
        else
        {
            for (size_t t = 0; t < num_thread; t++)
            {
                left_point_vector[t] = r_queue->top().left_point;
                right_point_vector[t] = r_queue->top().right_point;
                r_queue->pop();
                thread_vector.push_back(std::thread(&AGMND::parallel_calculate_new_point, this, left_point_vector[t], right_point_vector[t], m, std::ref(new_point_vector[t])));
            }
            for (size_t t = 0; t < num_thread; t++)
            {
                thread_vector[t].join();
                insert_to_map(new_point_vector[t].first, new_point_vector[t].second.z, 0, new_point_vector[t].second.der_z);
                update_result_p(left_point_vector[t], right_point_vector[t], new_point_vector[t]);
            }
            thread_vector.clear();
        }
        recalc_M_Max();
        calculate_m();
        recalc_R();
    }
    clear_solution();
    res.total_runtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();
    return res;
}


result AGMND::get_result()
{
    return res;
}
