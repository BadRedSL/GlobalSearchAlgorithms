#include <iostream>
#include <map>
#include <cmath>
#include <climits>
#include <queue>
#include <utility>
#include <chrono>
#include <thread>
#include <fstream>

#include "Point.h"
#include "Functions.h"
#include "ForR.h"

// последовательный алгоритм (без функций)
// sequential algorithm (no functions)

double minValue(const double a, const double b, const double r, const double accuracy,
    size_t nMax, double (*fnc)(double), double& result)
{
    std::ofstream out;          // поток для записи
    out.open(".\\standart_algorithm.txt"); // окрываем файл для записи

    std::map<Point, double> set;
    auto* rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;
    std::priority_queue<double, std::vector<double>, std::less<double>> mQueue;

    double minX;
    (fnc(a) < fnc(b)) ? minX = a : minX = b;
    double min;
    (fnc(a) < fnc(b)) ? min = fnc(a) : min = fnc(b);

    Point rightPoint(b, fnc(b));
    Point leftPoint(a, fnc(a));
    Point newPoint(b, fnc(b));
    set.emplace(leftPoint, 0);
    auto curIter = set.emplace(rightPoint, 0);

    double M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));
    double prevM = M;

    mQueue.push(M);

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) + ((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z) /
        (m * (rightPoint.x - leftPoint.x))) - 2 * (rightPoint.z + leftPoint.z);

    ForR newR(R, rightPoint, leftPoint);
    (*rQueue).push(newR);

    size_t count = 0;

    while ((count < nMax) && (rightPoint.x - leftPoint.x > accuracy))
    {
        double newX = ((rightPoint.x + leftPoint.x) / 2) -
            ((rightPoint.z - leftPoint.z) / (2 * m));

        double newY = fnc(newX);

        newPoint.setPoint(newX, newY);

        if (newY < min)
        {
            minX = newX;
            min = newY;
        }

        curIter = set.emplace(newPoint, 0);

        auto curTmp = curIter.first;
        auto prevTmp = std::prev(curTmp);

        if (mQueue.top() == abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x)))
        {
            mQueue.pop();
        }

        for (size_t i = 0; i < 2; ++i, ++curTmp, ++prevTmp)
        {
            M = abs(((*curTmp).first.z - (*prevTmp).first.z) /
                ((*curTmp).first.x - (*prevTmp).first.x));

            mQueue.push(M);
        }

        if (mQueue.top() > 0)
        {
            m = r * mQueue.top();
        }
        else
        {
            m = 1;
        }

        if (mQueue.top() == prevM)
        {
            curTmp = curIter.first;
            prevTmp = std::prev(curTmp);

            for (size_t i = 0; i < 2; ++i, ++curTmp, ++prevTmp)
            {
                R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
                    ((((*curTmp).first.z - (*prevTmp).first.z) *
                        ((*curTmp).first.z - (*prevTmp).first.z)) /
                        (m * (*curTmp).first.x - (*prevTmp).first.x)) -
                    (2 * (*curTmp).first.z + (*prevTmp).first.z);

                newR.setForR(R, (*curTmp).first, (*prevTmp).first);
                if (i == 0)
                {
                    (*rQueue).pop();
                }
                (*rQueue).push(newR);
            }
        }
        else
        {
            delete rQueue;
            rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

            for (curTmp = (++set.begin()), prevTmp = set.begin(); curTmp != set.end(); ++curTmp, ++prevTmp)
            {
                R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
                    ((((*curTmp).first.z - (*prevTmp).first.z) *
                        ((*curTmp).first.z - (*prevTmp).first.z)) /
                        (m * (*curTmp).first.x - (*prevTmp).first.x)) -
                    (2 * (*curTmp).first.z + (*prevTmp).first.z);

                newR.setForR(R, (*curTmp).first, (*prevTmp).first);
                (*rQueue).push(newR);
            }
        }
        rightPoint = (*rQueue).top().rightPoint;
        leftPoint = (*rQueue).top().leftPoint;
        prevM = mQueue.top();
        //std::cout << count << std::endl;
        ++count;
        out << "итерация: " << count << " значение x: " << minX << " верхнее R : " << (*rQueue).top().R << std::endl;
    }

    //std::cout << "Х при котором достигается минимум: " << minX << std::endl;

    //std::cout << "Время работы алгоритма: " << sec.count() << " сек " << std::endl;
    std::cout << "Число итераций последовательного алгоритма: " << count << std::endl;
    result = minX;
    return minX;
}


// параллельный (по отразкам) алгоритм (использует последовательный алгоритм)
// parallel (by reflections) algorithm (uses a sequential algorithm)

double inefficientParallelAlgorithm(size_t numThread, double leftBorder, double rightBorder,
    double r, double accuracy, size_t maxIteration, double (*fnc)(double))
{
    double border = leftBorder;
    double step = (rightBorder - leftBorder) / numThread;
    std::vector<std::thread> vThread;
    std::vector<double> vResult(numThread, 0);
    double minX = INT_MAX;
    double minY = INT_MAX;

    auto start = std::chrono::system_clock::now();

    for (size_t i = 0; i < numThread; ++i)
    {
        vThread.push_back(std::thread(minValue, border, border + step, r, accuracy, maxIteration, fnc, std::ref(vResult[i])));
        border += step;

    }

    for (size_t i = 0; i < vThread.size(); ++i)
    {
        vThread[i].join();
    }

    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    for (const auto& i : vResult)
    {
        //std::cout << "Полученная точка: " << i << std::endl;
        if (fnc(i) < minY)
        {
            minY = fnc(i);
            minX = i;
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> sec = end - start;
    //std::cout << "Время работы параллельного алгоритма: " << sec.count() << " сек " << std::endl;
    //std::cout << "Х при котором достигается минимум: " << minX << std::endl;

    return minX;
}

/*=====================================================================================================================================*/

void newPointCalculator(const Point leftPoint, const Point rightPoint, double (*fnc)(double), const double m, Point& newPoint)
{
    double newX = ((rightPoint.x + leftPoint.x) / 2) -
        ((rightPoint.z - leftPoint.z) / (2 * m));
    double newY = fnc(newX);

    newPoint.setPoint(newX, newY);
}

void newMCalculator(const Point leftPoint, const Point rightPoint, double& M)
{
    M = abs((rightPoint.z - leftPoint.z) /
        (rightPoint.x - leftPoint.x));
}

void newRCalculator(const Point leftPoint, const Point rightPoint, const double m, double& R)
{
    R = (m * (rightPoint.x - leftPoint.x)) +
        (((rightPoint.z - leftPoint.z) *
            (rightPoint.z - leftPoint.z)) /
            (m * rightPoint.x - leftPoint.x)) -
        (2 * rightPoint.z + leftPoint.z);
}

void allParallelCalculator(const Point leftPoint, const Point rightPoint, const double r, double (*fnc)(double), bool recountM, double maxM, Point& newPoint,
    double& M1, double& M2, double& R1, double& R2, ForR& ForR1, ForR& ForR2, double& vm)
{
    newPointCalculator(leftPoint, rightPoint, fnc, vm, newPoint);
    newMCalculator(leftPoint, newPoint, M1);
    newMCalculator(newPoint, rightPoint, M2);
    double M;
    if (recountM)
    {
        M = (M1 > M2) ? M1 : M2;
    }
    else
    {
        M = (M1 > M2) ? M1 : M2;
        M = (maxM > M) ? maxM : M;
    }
    vm = (M > 0) ? r * M : 1;
    newRCalculator(leftPoint, newPoint, vm, R1);
    ForR1.setForR(R1, newPoint, leftPoint);
    newRCalculator(newPoint, rightPoint, vm, R2);
    ForR2.setForR(R2, rightPoint, newPoint);
}

double functionNewMinValue(const double a, const double b,
    const double r, const double accuracy, size_t nMax,
    double (*fnc)(double))
{
    std::ofstream out;          // поток для записи
    out.open(".\\functional_algorithm.txt"); // окрываем файл для записи

    bool changeM = false;
    std::map<Point, double> set;
    auto* rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

    double minX = (fnc(a) < fnc(b)) ? a : minX = b;
    double min = (fnc(a) < fnc(b)) ? fnc(a) : min = fnc(b);

    Point rightPoint(b, fnc(b));
    Point leftPoint(a, fnc(a));
    Point newPoint(b, fnc(b));

    set.emplace(leftPoint, 0);
    set.emplace(rightPoint, 0);

    double M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) +
        ((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z) /
            (m * (rightPoint.x - leftPoint.x))) -
        2 * (rightPoint.z + leftPoint.z);

    ForR tmpR(R, rightPoint, leftPoint);
    (*rQueue).push(tmpR);

    size_t count = 0;

    while ((count < nMax) && (rightPoint.x - leftPoint.x > accuracy))
    {
        rightPoint = (*rQueue).top().rightPoint;
        leftPoint = (*rQueue).top().leftPoint;

        newPointCalculator(leftPoint, rightPoint, fnc, m, newPoint);
        set.emplace(newPoint, 0);

        double tmpM = 0;
        double prevM = M;

        //============NEW===========================
        newMCalculator(leftPoint, rightPoint, tmpM);
        if (M == tmpM)
        {
            M = -1;
            for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
            {
                newMCalculator((*prev).first, (*cur).first, tmpM);
                M = (tmpM > M) ? tmpM : M;
            }
        }
        else
        {
            newMCalculator(leftPoint, newPoint, tmpM);
            M = (tmpM > M) ? tmpM : M;
            newMCalculator(newPoint, rightPoint, tmpM);
            M = (tmpM > M) ? tmpM : M;
            changeM = (M == prevM) ? false : true;
        }
        //============NEW===========================
        m = (M > 0) ? r * M : 1;

        if (!changeM)
        {
            (*rQueue).pop();
            newRCalculator(leftPoint, newPoint, m, R);
            tmpR.setForR(R, newPoint, leftPoint);
            (*rQueue).push(tmpR);
            newRCalculator(newPoint, rightPoint, m, R);
            tmpR.setForR(R, rightPoint, newPoint);
            (*rQueue).push(tmpR);
        }
        else
        {
            delete rQueue;
            rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

            for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
            {
                newRCalculator((*prev).first, (*cur).first, m, R);
                tmpR.setForR(R, (*cur).first, (*prev).first);
                (*rQueue).push(tmpR);
            }
        }
        //std::cout << count << std::endl;
        ++count;
        if (newPoint.z < min)
        {
            min = newPoint.z;
            minX = newPoint.x;
        }
        out << "итерация: " << count << " значение x: " << minX << " верхнее R : " << (*rQueue).top().R << std::endl;
    }
    std::cout << "Число итераций функционального последовательного алгоритма: " << count << std::endl;
    return minX;
}

double parallelNewMinValue(size_t numThread, const double a, const double b,
    const double r, const double accuracy, size_t nMax,
    double (*fnc)(double))
{
    std::ofstream out;          // поток для записи
    out.open(".\\parallel_algorithm.txt"); // окрываем файл для записи

    std::vector<std::thread> vThread;
    std::vector<Point> vNewPoint(numThread);
    std::vector<double> vM1(numThread);
    std::vector<double> vM2(numThread);
    std::vector<double> vR1(numThread);
    std::vector<double> vR2(numThread);
    std::vector<ForR> vForR1(numThread);
    std::vector<ForR> vForR2(numThread);
    std::vector<double> vm(numThread);

    bool changeM = false;
    bool exitFlag = false;
    std::map<Point, double> set;
    auto* rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

    double minX = (fnc(a) < fnc(b)) ? a : minX = b;
    double min = (fnc(a) < fnc(b)) ? fnc(a) : min = fnc(b);

    Point rightPoint(b, fnc(b));
    Point leftPoint(a, fnc(a));
    Point newPoint(b, fnc(b));

    set.emplace(leftPoint, 0);
    set.emplace(rightPoint, 0);

    double M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) +
        ((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z) /
            (m * (rightPoint.x - leftPoint.x))) -
        2 * (rightPoint.z + leftPoint.z);

    ForR tmpR(R, rightPoint, leftPoint);
    (*rQueue).push(tmpR);

    size_t count = 0;

    while (count < nMax) {

        if (numThread >= (*rQueue).size())
        {
            rightPoint = (*rQueue).top().rightPoint;
            leftPoint = (*rQueue).top().leftPoint;

            if (rightPoint.x - leftPoint.x <= accuracy)
            {
                break;
            }

            newPointCalculator(leftPoint, rightPoint, fnc, m, newPoint);
            set.emplace(newPoint, 0);

            double tmpM = 0;
            double prevM = M;

            //============NEW===========================
            newMCalculator(leftPoint, rightPoint, tmpM);
            if (M == tmpM)
            {
                M = -1;
                for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
                {
                    newMCalculator((*prev).first, (*cur).first, tmpM);
                    M = (tmpM > M) ? tmpM : M;
                }
            }
            else
            {
                newMCalculator(leftPoint, newPoint, tmpM);
                M = (tmpM > M) ? tmpM : M;
                newMCalculator(newPoint, rightPoint, tmpM);
                M = (tmpM > M) ? tmpM : M;
            }
            //============NEW===========================
            changeM = (M == prevM) ? false : true;
            m = (M > 0) ? r * M : 1;
            for (auto& iter : vm)
            {
                iter = m;
            }

            if (!changeM)
            {
                (*rQueue).pop();
                newRCalculator(leftPoint, newPoint, m, R);
                tmpR.setForR(R, newPoint, leftPoint);
                (*rQueue).push(tmpR);
                newRCalculator(newPoint, rightPoint, m, R);
                tmpR.setForR(R, rightPoint, newPoint);
                (*rQueue).push(tmpR);
            }
            else
            {
                delete rQueue;
                rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

                for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
                {
                    newRCalculator((*prev).first, (*cur).first, m, R);
                    tmpR.setForR(R, (*cur).first, (*prev).first);
                    (*rQueue).push(tmpR);
                }
            }
            ++count;
            if (newPoint.z < min)
            {
                min = newPoint.z;
                minX = newPoint.x;
            }
            out << "итерация: " << count << " значение x: " << minX << " верхнее R : " << (*rQueue).top().R << std::endl;
        }
        else
        {
            double prevM = M; //надо посмотреть внимательно
            double tmpM = 0;
            bool recountM = false;
            for (size_t i = 0; i < numThread; ++i)
            {
                rightPoint = (*rQueue).top().rightPoint;
                leftPoint = (*rQueue).top().leftPoint;

                //====================NEW==============================
                newMCalculator(leftPoint, rightPoint, tmpM);

                if (M == tmpM)
                {
                    recountM = true;
                }
                //====================NEW==============================

                if (rightPoint.x - leftPoint.x <= accuracy)
                {
                    exitFlag = true;
                }
                vThread.push_back(std::thread(allParallelCalculator, leftPoint, rightPoint, r, fnc, recountM, M, std::ref(vNewPoint[i]), std::ref(vM1[i]), 
                    std::ref(vM2[i]), std::ref(vR1[i]), std::ref(vR2[i]), std::ref(vForR1[i]), std::ref(vForR2[i]), std::ref(vm[i])));
                (*rQueue).pop();
            }
            for (size_t i = 0; i < numThread; ++i)
            {
                vThread[i].join();
            }
            vThread.clear();
            if (exitFlag)
            {
                break;
            }

            //====================NEW==============================
            if (recountM)
            {
                for (size_t i = 0; i < numThread; ++i)
                {
                    set.emplace(vNewPoint[i], 0);
                    (*rQueue).push(vForR1[i]);
                    (*rQueue).push(vForR2[i]);
                }
                M = -1;
                for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
                {
                    newMCalculator((*prev).first, (*cur).first, tmpM);
                    M = (tmpM > M) ? tmpM : M;
                }
                recountM = false;
            }
            else
            {
                for (size_t i = 0; i < numThread; ++i)
                {
                    set.emplace(vNewPoint[i], 0);
                    M = (vM1[i] > M) ? vM1[i] : M;
                    M = (vM2[i] > M) ? vM2[i] : M;
                    (*rQueue).push(vForR1[i]);
                    (*rQueue).push(vForR2[i]);
                }
            }
            //====================NEW==============================

            for (auto& iter : vm)
            {
                iter = (M>0) ? r*M : 1;
            }
            if (prevM != M)
            {
                delete rQueue;
                rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

                for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end(); ++cur, ++prev)
                {
                    newRCalculator((*prev).first, (*cur).first, vm[0], R);
                    tmpR.setForR(R, (*cur).first, (*prev).first);
                    (*rQueue).push(tmpR);
                }
            }
            count += numThread;
            for (const auto iter : vNewPoint)
            {
                if (iter.z < min)
                {
                    min = iter.z;
                    minX = iter.x;
                }
            } 
            out << "итерация: " << count << " значение x: " << minX << " верхнее R : " << (*rQueue).top().R << std::endl;
        }
    }
    std::cout << "Число итераций параллельного алгоритма: " << count << std::endl;
    return minX;
}
