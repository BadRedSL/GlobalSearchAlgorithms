#pragma once


struct Point
{
public:

	double x;
	double z;

	Point();
	Point(double x, double z);
	Point(const Point& a);

	void setPoint(const double x, const double z);

	friend bool operator< (const Point& a, const Point& b);

	friend bool operator> (const Point& a, const Point& b);

};

