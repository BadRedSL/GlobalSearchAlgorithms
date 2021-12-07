#include <iostream>
#include <map>
#include <cmath>
#include <climits>
#include <queue>
#include <utility>
#include <chrono>
#include <thread>

#include "Point.h"
#include "Functions.h"
#include "ForR.h"

// последовательный алгоритм (без функций)
// sequential algorithm (no functions)

double minValue(const double a, const double b, const double r, const double accuracy,
    size_t nMax, double (*fnc)(double), double& result)
{
    auto start = std::chrono::system_clock::now();

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

        for (size_t i = 0; i < 2; ++i)
        {
            M = abs(((*curTmp).first.y - (*prevTmp).first.y) /
                ((*curTmp).first.x - (*prevTmp).first.x));

            mQueue.push(M);
            ++curTmp;
            ++prevTmp;
        }

        if (mQueue.top() > 0)
        {
            m = r * mQueue.top();
        }
        else
        {
            m = 1;
        }

        if (M == prevM)
        {
            curTmp = curIter.first;
            prevTmp = std::prev(curTmp);

            for (size_t i = 0; i < 2; ++i)
            {
                R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
                    ((((*curTmp).first.y - (*prevTmp).first.y) *
                        ((*curTmp).first.y - (*prevTmp).first.y)) /
                        (m * (*curTmp).first.x - (*prevTmp).first.x)) -
                    (2 * (*curTmp).first.y + (*prevTmp).first.y);

                ++curTmp;
                ++prevTmp;

                newR.setForR(R, (*curTmp).first, (*prevTmp).first);
                (*rQueue).pop();
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
        prevM = M;
        std::cout << count << std::endl;
        ++count;
    }

    //std::cout << "Х при котором достигается минимум: " << minX << std::endl;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> sec = end - start;

    //std::cout << "Время работы алгоритма: " << sec.count() << " сек " << std::endl;

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


// функции для последовательного алгоритма
// functions for sequential algorithm

void setCalculateNewPoint(std::pair<std::map<Point, double>::iterator, bool>& curIter, std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>* rQueue, std::map<Point, double>& set,
    double& min, double& minX, const double m, double (*fnc)(double))
{
    Point rightPoint = (*rQueue).top().rightPoint;
    Point leftPoint = (*rQueue).top().leftPoint;
    Point newPoint;

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
}

void calculateM(std::pair<std::map<Point, double>::iterator, bool>& curIter, std::priority_queue<double, std::vector<double>,
    std::less<double>> &mQueue, bool& changeM, const Point rightPoint, const Point leftPoint)
{
    auto curTmp = curIter.first;
    auto prevTmp = std::prev(curTmp);
    double M;

    if (mQueue.top() == abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x)))
    {
        mQueue.pop();
    }

    for (size_t i = 0; i < 2; ++i)
    {
        M = abs(((*curTmp).first.y - (*prevTmp).first.y) /
            ((*curTmp).first.x - (*prevTmp).first.x));
        if (!mQueue.empty())
        {
            changeM = (M == mQueue.top()) ? true : false;
        }
        mQueue.push(M);
        ++curTmp;
        ++prevTmp;
    }
}

void calculateR(std::pair<std::map<Point, double>::iterator, bool>& curIter, std::vector<ForR>& newR, const double m)
{
    auto curTmp = curIter.first;
    auto prevTmp = std::prev(curTmp);
    double R;

    for (size_t i = 0; i < 2; ++i)
    {
        R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
            ((((*curTmp).first.y - (*prevTmp).first.y) *
                ((*curTmp).first.y - (*prevTmp).first.y)) /
                (m * (*curTmp).first.x - (*prevTmp).first.x)) -
            (2 * (*curTmp).first.y + (*prevTmp).first.y);

        ++curTmp;
        ++prevTmp;

        newR[i].setForR(R, (*curTmp).first, (*prevTmp).first);

    }
}

void reNewR(std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>** rQueue, std::map<Point, double>& set, const double m)
{
    delete* rQueue;
    *rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;
    double R;
    ForR newR;

    for (auto curTmp = (++set.begin()), prevTmp = set.begin(); curTmp != set.end(); ++curTmp, ++prevTmp)
    {
        R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
            ((((*curTmp).first.y - (*prevTmp).first.y) *
                ((*curTmp).first.y - (*prevTmp).first.y)) /
                (m * (*curTmp).first.x - (*prevTmp).first.x)) -
            (2 * (*curTmp).first.y + (*prevTmp).first.y);

        newR.setForR(R, (*curTmp).first, (*prevTmp).first);
        (**rQueue).push(newR);
    }
}


// последовательный алгоритм на функциях
// sequential algorithm on functions

double functionMinValue(const double a, const double b, const double r, const double accuracy, size_t nMax, double (*fnc)(double))
{
    bool changeM = false;
    std::vector<ForR> newR(2);
    std::map<Point, double> set;
    std::priority_queue<double, std::vector<double>, std::less<double>> mQueue;
    auto* rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

    double minX = (fnc(a) < fnc(b)) ? a : minX = b;
    double min = (fnc(a) < fnc(b)) ? fnc(a) : min = fnc(b);

    Point rightPoint(b, fnc(b));
    Point leftPoint(a, fnc(a));
    Point newPoint(b, fnc(b));

    set.emplace(leftPoint, 0);
    auto curIter = set.emplace(rightPoint, 0);

    double M = abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x));
    mQueue.push(M);

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
        setCalculateNewPoint(curIter, rQueue, set, min, minX, m, fnc);

        calculateM(curIter, mQueue, changeM, rightPoint, leftPoint);

        if (mQueue.top() > 0)
        {
            m = r * mQueue.top();
        }
        else
        {
            m = 1;
        }

        if (changeM)
        {
            calculateR(curIter, newR, m);
            (*rQueue).pop();
            (*rQueue).push(newR[0]);
            (*rQueue).pop();
            (*rQueue).push(newR[1]);
        }
        else
        {
            reNewR(&rQueue, set, m);
        }

        rightPoint = (*rQueue).top().rightPoint;
        leftPoint = (*rQueue).top().leftPoint;
        std::cout << count << std::endl;
        ++count;
    }
    return minX;
}

// функции для параллельного алгоритма
// functions for the parallel algorithm

void parallelCalculateNewPoint(const ForR top, double (*fnc)(double), const double m, Point& newPoint)
{
    Point rightPoint = top.rightPoint;
    Point leftPoint = top.leftPoint;

    double newX = ((rightPoint.x + leftPoint.x) / 2) -
        ((rightPoint.y - leftPoint.y) / (2 * m));
    double newY = fnc(newX);

    newPoint.setPoint(newX, newY);
}

void parallelCalculateM(const std::pair<std::map<Point, double>::iterator, bool> iter, double& firstM, double& secondM)
{
    auto curTmp = iter.first;
    auto prevTmp = std::prev(curTmp);
    double M;

    for (size_t i = 0; i < 2; ++i)
    {
        M = abs(((*curTmp).first.y - (*prevTmp).first.y) /
            ((*curTmp).first.x - (*prevTmp).first.x));

        if (i == 0)
        {
            firstM = M;
        }
        if (i == 1)
        {
            secondM = M;
        }

        ++curTmp;
        ++prevTmp;
    }
}

void parallelCalculateR(const double r, const double M, const std::pair<std::map<Point, double>::iterator, bool> iter, ForR& firstR, ForR& secondR)
{
    double m;
    double R;
    ForR newR;

    if (M > 0)
    {
        m = r * M;
    }
    else
    {
        m = 1;
    }

    auto curTmp = iter.first;
    auto prevTmp = std::prev(curTmp);

    for (size_t i = 0; i < 2; ++i)
    {
        R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
            ((((*curTmp).first.y - (*prevTmp).first.y) *
                ((*curTmp).first.y - (*prevTmp).first.y)) /
                (m * (*curTmp).first.x - (*prevTmp).first.x)) -
            (2 * (*curTmp).first.y + (*prevTmp).first.y);

        ++curTmp;
        ++prevTmp;

        newR.setForR(R, (*curTmp).first, (*prevTmp).first);
        if (i == 0)
        {
            firstR = newR;
        }
        if (i == 1)
        {
            secondR = newR;
        }
        
    }
}


// параллельный алгоритм (в разработке)
// parallel algorithm (in development)

//double parallelMinValue(size_t numThread, const double a, const double b, const double r, const double accuracy, size_t nMax, double (*fnc)(double))
//{
//    std::vector<std::thread> vThread(numThread);
//    std::vector<Point> vPoint(numThread);
//    std::vector<std::pair<std::map<Point, double>::iterator, bool>> vIter(numThread);
//    std::vector<double> firstM(numThread);
//    std::vector<double> secondM(numThread);
//    std::vector<ForR> firstR(numThread);
//    std::vector<ForR> secondR(numThread);
//    
//
//    bool changeM = false;
//    std::vector<ForR> newR(2);
//    std::map<Point, double> set;
//    std::priority_queue<double, std::vector<double>, std::less<double>> mQueue;
//    auto* rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;
//
//    double minX = (fnc(a) < fnc(b)) ? a : minX = b;
//    double min = (fnc(a) < fnc(b)) ? fnc(a) : min = fnc(b);
//
//    Point rightPoint(b, fnc(b));
//    Point leftPoint(a, fnc(a));
//    Point newPoint(b, fnc(b));
//
//    set.emplace(leftPoint, 0);
//    auto curIter = set.emplace(rightPoint, 0);
//
//    double M = abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x));
//    mQueue.push(M);
//
//    double m = (M > 0) ? r * M : 1;
//
//    double R = (m * (rightPoint.x - leftPoint.x)) +
//        ((rightPoint.y - leftPoint.y) * (rightPoint.y - leftPoint.y) /
//            (m * (rightPoint.x - leftPoint.x))) -
//        2 * (rightPoint.y + leftPoint.y);
//    ForR tmpR(R, rightPoint, leftPoint);
//    (*rQueue).push(tmpR);
//
//    size_t count = 0;
//
//    while ((count < nMax) && (rightPoint.x - leftPoint.x > accuracy))
//    {
//        if ((*rQueue).size() >= numThread)
//        {
//            for (size_t i = 0; i < numThread; ++i)
//            {
//                vThread[i] = (std::thread(parallelCalculateNewPoint, (*rQueue).top(), fnc, m, std::ref(vPoint[i]))); //нужно добавить флаг выхода по accuracy
//                (*rQueue).pop();
//            }
//            for (size_t i = 0; i < numThread; ++i)
//            {
//                vThread[i].join();
//                min = (vPoint[i].y < min) ? vPoint[i].y : min;
//                minX = (vPoint[i].y < min) ? vPoint[i].x : minX; //совместить функции
//                vIter[i] = set.emplace(vPoint[i], 0);
//            }
//
//            for (size_t i = 0; i < numThread; ++i)
//            {
//                /*if (mQueue.top() == abs((rightPoint.y - leftPoint.y) / (rightPoint.x - leftPoint.x))) // убрать очередь с М
//                {
//                    mQueue.pop();
//                }*/
//                vThread[i] = (std::thread(parallelCalculateM, vIter[i], std::ref(firstM[i]), std::ref(secondM[i])));
//            }
//            for (size_t i = 0; i < numThread; ++i)
//            {
//                vThread[i].join();
//                mQueue.push(firstM[i]);
//                mQueue.push(secondM[i]);
//                if (!mQueue.empty())
//                {
//                    changeM = (firstM[i] == mQueue.top() || secondM[i] == mQueue.top()) ? true : false;
//                }
//            }
//            
//            if (changeM)
//            {
//                for (size_t i = 0; i < numThread; ++i)
//                {
//                    vThread[i] = (std::thread(parallelCalculateR, ));
//                }
//                for (size_t i = 0; i < numThread; ++i)
//                {
//                    vThread[i].join();
//                   
//                }
//            }
//            else
//            {
//                (*rQueue).pop();
//                (*rQueue).push(newR);
//            }
//
//        }
//        else
//        {
//            setCalculateNewPoint(curIter, rQueue, set, min, minX, m, fnc);
//
//            calculateM(curIter, mQueue, changeM, rightPoint, leftPoint);
//
//            if (mQueue.top() > 0)
//            {
//                m = r * mQueue.top();
//            }
//            else
//            {
//                m = 1;
//            }
//
//            if (changeM)
//            {
//                calculateR(curIter, newR, m);
//                (*rQueue).pop();
//                (*rQueue).push(newR[0]);
//                (*rQueue).pop();
//                (*rQueue).push(newR[1]);
//            }
//            else
//            {
//                reNewR(&rQueue, set, m);
//            }
//
//            rightPoint = (*rQueue).top().rightPoint;
//            leftPoint = (*rQueue).top().leftPoint;
//            std::cout << count << std::endl;
//            ++count;
//        }
//    }
//    return minX;
//}

int main()
{
    setlocale(LC_ALL, "rus");

    //std::cout << inefficientParallelAlgorithm(1, 2.7, 7.5, 4.29, 0.0001, 10000, function1);
    std::cout << parallelMinValue(2, 2.7, 7.5, 4.29, 0.0001, 10000, function1);
    return 0;
}
