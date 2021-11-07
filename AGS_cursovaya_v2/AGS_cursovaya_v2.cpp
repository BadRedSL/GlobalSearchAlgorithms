#include <iostream>
#include <map>
#include <cmath>
#include <climits>
#include <queue>
#include <utility>
#include <chrono>

#include "Point.h"
#include "Functions.h"
#include "ForR.h"


double minValueOld(const double a, const double b, const double r, const double accuracy,
	size_t nMax, double (*fnc)(double))
{
	auto start = std::chrono::system_clock::now();

	std::map<Point, double> set;
	std::priority_queue<double, std::vector<double>, std::less<double>> rQueue;
	std::priority_queue<double, std::vector<double>, std::less<double>> mQueue;

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

	rQueue.push(R);

	size_t count = 0;
	
	while ((count < nMax) && (rightPoint.x - leftPoint.x > accuracy))
	{
		double newX = ((rightPoint.x + leftPoint.x) / 2) - 
			((rightPoint.y - leftPoint.y) / (2 * m));

		double newY = fnc(newX);

		newPoint.setPoint(newX, newY);

		if (newY < min)
		{
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
				rQueue.push(R);

				if (rQueue.top() == R)
				{
					rightPoint = (*curTmp).first;
					leftPoint = (*prevTmp).first;
				}

			}
		}
		else
		{
			while (!rQueue.empty())
			{
				rQueue.pop();
			}

			for (curTmp = (++set.begin()), prevTmp = set.begin(); curTmp != set.end(); ++curTmp, ++prevTmp)
			{
				R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
					((((*curTmp).first.y - (*prevTmp).first.y) *
						((*curTmp).first.y - (*prevTmp).first.y)) /
						(m * (*curTmp).first.x - (*prevTmp).first.x)) -
					(2 * (*curTmp).first.y + (*prevTmp).first.y);

				rQueue.push(R);

				if (rQueue.top() == R)
				{
					rightPoint = (*curTmp).first;
					leftPoint = (*prevTmp).first;
				}

			}

		}

		prevM = M;
		std::cout << count << std::endl;
		++count;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> sec = end - start;
	std::cout << "Время работы алгоритма: " << sec.count() << " сек " << std::endl;

	return min;
}


double minValue(const double a, const double b, const double r, const double accuracy,
	size_t nMax, double (*fnc)(double))
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
	std::cout << "Время работы алгоритма: " << sec.count() << " сек " << std::endl;

	return minX;
}

int main()
{
	setlocale(LC_ALL, "rus");

	//std::cout << minValue(2.7, 7.5, 4.29, 0.001, 1000, function1); 
	//std::cout << minValue(0, 10, 67, 0.001, 1000, function2); 
	//std::cout << minValue(0, 1.2, 36, 0.001, 1000, function3); 
	//std::cout << minValue(-10, 10, 2.5, 0.001, 1000, function4); 
	//std::cout << minValue(2.7, 7.5, 6, 0.001, 1000, function5); 
	//std::cout << minValue(0, 4, 6.5, 0.001, 1000, function7);
	//std::cout << minValue(-5, 5, 6.5, 0.001, 1000, function8);
	//std::cout << minValue(0, 6.5, 4, 0.00001, 1000, function9);
	//std::cout << minValue(-3, 3, 85, 0.001, 1000, function10);

	//std::cout << minValue(0, 6.5, 4, 0.00000000001, 12000, function9);//329 sec debug-mode
	//std::cout << minValueOld(0, 6.5, 4, 0.00000000001, 12000, function9);//408 sec debug-mode

	std::cout << minValue(0, 6.5, 4, 0.00000000001, 12000, function9);//9 sec release-mode
	//std::cout << minValueOld(0, 6.5, 4, 0.00000000001, 12000, function9);//10,5 sec release-mode

	return 0;
}
