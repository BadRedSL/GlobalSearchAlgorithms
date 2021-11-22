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

double minValue(const double a, const double b, const double r, const double accuracy,
	size_t nMax, double (*fnc)(double), double &result)
{
	auto start = std::chrono::system_clock::now();

	std::map<Point, double> set;
	auto *rQueue = new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;
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
				(*rQueue).push(newR);

				
				rightPoint = (*rQueue).top().rightPoint;
				leftPoint = (*rQueue).top().leftPoint;

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

				rightPoint = (*rQueue).top().rightPoint;
				leftPoint = (*rQueue).top().leftPoint;
				
			}

		}

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

int main()
{
	setlocale(LC_ALL, "rus");

	size_t numThread = 12;
	double leftBorder = 2.7;
	double rightBorder = 7.5;
	double r = 4.29;
	double accuracy = 0.0001;
	size_t maxIteration = 10000;
	double (*fnc)(double) = function1;
	
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

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	for (const auto &i: vResult)
	{
		std::cout << "Полученная точка: " << i << std::endl;
		if (fnc(i) < minY)
		{
			minY = fnc(i);
			minX = i;
		}
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> sec = end - start;
	std::cout << "Время работы параллельного алгоритма: " << sec.count() << " сек " << std::endl;

	std::cout <<"Х при котором достигается минимум: " << minX << std::endl;

	return 0;
}
