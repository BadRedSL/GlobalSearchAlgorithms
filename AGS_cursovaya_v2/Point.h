#pragma once


class Point
{
public:

	double x;
	double y;

	Point();
	Point(double x, double y);
	Point(const Point& a);

	void setPoint(const double x, const double y);

	friend bool operator< (const Point& a, const Point& b);

	friend bool operator> (const Point& a, const Point& b);

};

