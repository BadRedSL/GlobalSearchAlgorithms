#pragma once
#include "Point.h"

class ForR
{
public:

	double R;
	Point rightPoint;
	Point leftPoint;

	ForR(double R, Point rightPoint, Point leftPoint);

	void setForR(double R, Point rightPoint, Point leftPoint);

	friend bool operator< (const ForR& a, const ForR& b);

	friend bool operator> (const ForR& a, const ForR& b);

};

