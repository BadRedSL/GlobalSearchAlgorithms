#pragma once
#include "Point.h"

struct ForR
{
public:

	double R;
	Point rightPoint;
	Point leftPoint;

	ForR();

	ForR(double R, Point rightPoint, Point leftPoint);

	void setForR(double R, Point rightPoint, Point leftPoint);

	friend bool operator< (const ForR& a, const ForR& b);

	friend bool operator> (const ForR& a, const ForR& b);

};

