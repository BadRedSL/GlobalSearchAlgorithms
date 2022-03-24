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

    double M = abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x));
    double prevM = M;

    mQueue.push(M);

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) + ((rightPoint.y - leftPoint.y) * (rightPoint.y - leftPoint.y) /
        (m * (rightPoint.x - leftPoint.x))) - 2 * (rightPoint.y + leftPoint.y);

    ForR newR(R, rightPoint, leftPoint);
    (*rQueue).push(newR);

    size_t count = 0;

    while ((count < nMax) && (rightPoint.x - leftPoint.x > accuracy))
    {
        double newX = ((rightPoint.x + leftPoint.x) / 2) -
            ((rightPoint.y - leftPoint.y) / (2 * m));

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

        if (mQueue.top() == abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x)))
        {
            mQueue.pop();
        }

        for (size_t i = 0; i < 2; ++i, ++curTmp, ++prevTmp)
        {
            M = abs(((*curTmp).first.y - (*prevTmp).first.y) /
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
                    ((((*curTmp).first.y - (*prevTmp).first.y) *
                        ((*curTmp).first.y - (*prevTmp).first.y)) /
                        (m * (*curTmp).first.x - (*prevTmp).first.x)) -
                    (2 * (*curTmp).first.y + (*prevTmp).first.y);

                newR.setForR(R, (*curTmp).first, (*prevTmp).first);
                if (i == 1)
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
                    ((((*curTmp).first.y - (*prevTmp).first.y) *
                        ((*curTmp).first.y - (*prevTmp).first.y)) /
                        (m * (*curTmp).first.x - (*prevTmp).first.x)) -
                    (2 * (*curTmp).first.y + (*prevTmp).first.y);

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
        ((rightPoint.y - leftPoint.y) / (2 * m));
    double newY = fnc(newX);

    newPoint.setPoint(newX, newY);
}

void newMCalculator(const Point leftPoint, const Point rightPoint, double& M)
{
    M = abs((rightPoint.y - leftPoint.y) /
        (rightPoint.x - leftPoint.x));
}

void newRCalculator(const Point leftPoint, const Point rightPoint, const double m, double& R)
{
    R = (m * (rightPoint.x - leftPoint.x)) +
        (((rightPoint.y - leftPoint.y) *
            (rightPoint.y - leftPoint.y)) /
            (m * rightPoint.x - leftPoint.x)) -
        (2 * rightPoint.y + leftPoint.y);
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

    double M = abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x));

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) +
        ((rightPoint.y - leftPoint.y) * (rightPoint.y - leftPoint.y) /
            (m * (rightPoint.x - leftPoint.x))) -
        2 * (rightPoint.y + leftPoint.y);

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
        if (newPoint.y < min)
        {
            min = newPoint.y;
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

    double M = abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x));

    double m = (M > 0) ? r * M : 1;

    double R = (m * (rightPoint.x - leftPoint.x)) +
        ((rightPoint.y - leftPoint.y) * (rightPoint.y - leftPoint.y) /
            (m * (rightPoint.x - leftPoint.x))) -
        2 * (rightPoint.y + leftPoint.y);

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
                changeM = (M == prevM) ? false : true;
            }
            //============NEW===========================
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
            if (newPoint.y < min)
            {
                min = newPoint.y;
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
                if (iter.y < min)
                {
                    min = iter.y;
                    minX = iter.x;
                }
            } 
            out << "итерация: " << count << " значение x: " << minX << " верхнее R : " << (*rQueue).top().R << std::endl;
        }
    }
    std::cout << "Число итераций параллельного алгоритма: " << count << std::endl;
    return minX;
}


int main()
{
    setlocale(LC_ALL, "rus");

    int numOfThreads = 1;
    int maxIteration = 10000;
    double accuracy = 0.0001;
    double r = 2.0;

    auto start11 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(2.7, 7.5, r, accuracy, maxIteration, function1) << std::endl;
    auto end11 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec11 = end11 - start11;
    std::cout << "Время работы функционального последовательного алгоритма (функция1): " << sec11.count() << " сек " << std::endl << std::endl;

    double result12 = 0;
    auto start12 = std::chrono::system_clock::now();
    std::cout << minValue(2.7, 7.5, r, accuracy, maxIteration, function1, result12) << std::endl;
    auto end12 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec12 = end12 - start12;
    std::cout << "Время работы обычного последовательного алгоритма (функция1): " << sec12.count() << " сек " << std::endl << std::endl;

    auto start13 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 2.7, 7.5, r, accuracy, maxIteration, function1) << std::endl;
    auto end13 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec13 = end13 - start13;
    std::cout << "Время работы параллельного алгоритма (функция1): " << sec13.count() << " сек " << std::endl << std::endl;


    auto start21 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(0, 10, r, accuracy, maxIteration, function2) << std::endl;
    auto end21 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec21 = end21 - start21;
    std::cout << "Время работы функционального последовательного алгоритма (функция2): " << sec21.count() << " сек " << std::endl << std::endl;

    double result22 = 0;
    auto start22 = std::chrono::system_clock::now();
    std::cout << minValue(0, 10, r, accuracy, maxIteration, function2, result22) << std::endl;
    auto end22 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec22 = end22 - start22;
    std::cout << "Время работы обычного последовательного алгоритма (функция2): " << sec22.count() << " сек " << std::endl << std::endl;

    auto start23 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 0, 10, r, accuracy, maxIteration, function2) << std::endl;
    auto end23 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec23 = end23 - start23;
    std::cout << "Время работы параллельного алгоритма (функция2): " << sec23.count() << " сек " << std::endl << std::endl;


    auto start31 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(0, 1.2, r, accuracy, maxIteration, function3) << std::endl;
    auto end31 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec31 = end31 - start31;
    std::cout << "Время работы функционального последовательного алгоритма (функция3): " << sec31.count() << " сек " << std::endl << std::endl;

    double result32 = 0;
    auto start32 = std::chrono::system_clock::now();
    std::cout << minValue(0, 1.2, r, accuracy, maxIteration, function3, result32) << std::endl;
    auto end32 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec32 = end22 - start22;
    std::cout << "Время работы обычного последовательного алгоритма (функция3): " << sec32.count() << " сек " << std::endl << std::endl;

    auto start33 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 0, 1.2, r, accuracy, maxIteration, function3) << std::endl;
    auto end33 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec33 = end33 - start33;
    std::cout << "Время работы параллельного алгоритма (функция3): " << sec33.count() << " сек " << std::endl << std::endl;


    auto start41 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(-10, 10, r, accuracy, maxIteration, function4) << std::endl;
    auto end41 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec41 = end41 - start41;
    std::cout << "Время работы функционального последовательного алгоритма (функция4): " << sec41.count() << " сек " << std::endl << std::endl;

    double result42 = 0;
    auto start42 = std::chrono::system_clock::now();
    std::cout << minValue(-10, 10, r, accuracy, maxIteration, function4, result42) << std::endl;
    auto end42 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec42 = end42 - start42;
    std::cout << "Время работы обычного последовательного алгоритма (функция4): " << sec42.count() << " сек " << std::endl << std::endl;

    auto start43 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, -10, 10, r, accuracy, maxIteration, function4) << std::endl;
    auto end43 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec43 = end43 - start43;
    std::cout << "Время работы параллельного алгоритма (функция4): " << sec43.count() << " сек " << std::endl << std::endl;


    auto start51 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(2.7, 7.5, r, accuracy, maxIteration, function5) << std::endl;
    auto end51 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec51 = end51 - start51;
    std::cout << "Время работы функционального последовательного алгоритма (функция5): " << sec51.count() << " сек " << std::endl << std::endl;

    double result52 = 0;
    auto start52 = std::chrono::system_clock::now();
    std::cout << minValue(2.7, 7.5, r, accuracy, maxIteration, function5, result52) << std::endl;
    auto end52 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec52 = end52 - start52;
    std::cout << "Время работы обычного последовательного алгоритма (функция5): " << sec52.count() << " сек " << std::endl << std::endl;

    auto start53 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 2.7, 7.5, r, accuracy, maxIteration, function5) << std::endl;
    auto end53 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec53 = end53 - start53;
    std::cout << "Время работы параллельного алгоритма (функция5): " << sec53.count() << " сек " << std::endl << std::endl;


    auto start71 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(0, 4, r, accuracy, maxIteration, function7) << std::endl;
    auto end71 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec71 = end71 - start71;
    std::cout << "Время работы функционального последовательного алгоритма (функция7): " << sec71.count() << " сек " << std::endl << std::endl;

    double result72 = 0;
    auto start72 = std::chrono::system_clock::now();
    std::cout << minValue(0, 4, r, accuracy, maxIteration, function7, result72) << std::endl;
    auto end72 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec72 = end72 - start72;
    std::cout << "Время работы обычного последовательного алгоритма (функция7): " << sec72.count() << " сек " << std::endl << std::endl;

    auto start73 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 0, 4, r, accuracy, maxIteration, function7) << std::endl;
    auto end73 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec73 = end73 - start73;
    std::cout << "Время работы параллельного алгоритма (функция7): " << sec73.count() << " сек " << std::endl << std::endl;


    auto start81 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(-5, 5, r, accuracy, maxIteration, function8) << std::endl;
    auto end81 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec81 = end81 - start81;
    std::cout << "Время работы функционального последовательного алгоритма (функция8): " << sec81.count() << " сек " << std::endl << std::endl;

    double result82 = 0;
    auto start82 = std::chrono::system_clock::now();
    std::cout << minValue(-5, 5, r, accuracy, maxIteration, function8, result82) << std::endl;
    auto end82 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec82 = end82 - start82;
    std::cout << "Время работы обычного последовательного алгоритма (функция8): " << sec82.count() << " сек " << std::endl << std::endl;

    auto start83 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, -5, 5, r, accuracy, maxIteration, function8) << std::endl;
    auto end83 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec83 = end83 - start83;
    std::cout << "Время работы параллельного алгоритма (функция8): " << sec83.count() << " сек " << std::endl << std::endl;


    auto start91 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(0, 6.5, r, accuracy, maxIteration, function9) << std::endl;
    auto end91 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec91 = end91 - start91;
    std::cout << "Время работы функционального последовательного алгоритма (функция9): " << sec91.count() << " сек " << std::endl << std::endl;

    double result92 = 0;
    auto start92 = std::chrono::system_clock::now();
    std::cout << minValue(0, 6.5, r, accuracy, maxIteration, function9, result92) << std::endl;
    auto end92 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec92 = end92 - start92;
    std::cout << "Время работы обычного последовательного алгоритма (функция9): " << sec92.count() << " сек " << std::endl << std::endl;

    auto start93 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, 0, 6.5, r, accuracy, maxIteration, function9) << std::endl;
    auto end93 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec93 = end93 - start93;
    std::cout << "Время работы параллельного алгоритма (функция9): " << sec93.count() << " сек " << std::endl << std::endl;


    auto start101 = std::chrono::system_clock::now();
    std::cout << functionNewMinValue(-3, 3, r, accuracy, maxIteration, function10) << std::endl;
    auto end101 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec101 = end101 - start101;
    std::cout << "Время работы функционального последовательного алгоритма (функция10): " << sec101.count() << " сек " << std::endl << std::endl;

    double result102 = 0;
    auto start102 = std::chrono::system_clock::now();
    std::cout << minValue(-3, 3, r, accuracy, maxIteration, function10, result102) << std::endl;
    auto end102 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec102 = end102 - start102;
    std::cout << "Время работы обычного последовательного алгоритма (функция10): " << sec102.count() << " сек " << std::endl << std::endl;

    auto start103 = std::chrono::system_clock::now();
    std::cout << parallelNewMinValue(numOfThreads, -3, 3, r, accuracy, maxIteration, function10) << std::endl;
    auto end103 = std::chrono::system_clock::now();
    std::chrono::duration<double> sec103 = end103 - start103;
    std::cout << "Время работы параллельного алгоритма (функция10): " << sec103.count() << " сек " << std::endl << std::endl;

    return 0;
}
