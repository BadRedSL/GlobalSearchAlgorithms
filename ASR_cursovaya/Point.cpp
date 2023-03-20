#include "Point.h"

Point::Point()
{
	this->x = 0;
	this->z = 0;
}

Point::Point(double x, double z)
{
	this->x = x;
	this->z = z;
}

Point::Point(const Point& a)
{
	this->x = a.x;
	this->z = a.z;
}

void Point::setPoint(const double x, const double z)
{
	this->x = x;
	this->z = z;
}

bool operator<(const Point& a, const Point& b)
{
	return (a.x < b.x);
}

bool operator>(const Point& a, const Point& b)
{
	return (a.x > b.x);
}
