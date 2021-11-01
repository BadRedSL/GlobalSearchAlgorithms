#include "Point.h"

Point::Point()
{
	this->x = 0;
	this->y = 0;
}

Point::Point(double x, double y)
{
	this->x = x;
	this->y = y;
}

Point::Point(const Point& a)
{
	this->x = a.x;
	this->y = a.y;
}

void Point::setPoint(const double x, const double y)
{
	this->x = x;
	this->y = y;
}

bool operator<(const Point& a, const Point& b)
{
	return (a.x < b.x);
}

bool operator>(const Point& a, const Point& b)
{
	return (a.x > b.x);
}
