#include "ForR.h"

ForR::ForR()
{
	R = 0;
}

ForR::ForR(double R, Point rightPoint, Point leftPoint)
{
	this->R = R;
	this->rightPoint = rightPoint;
	this->leftPoint = leftPoint;
}

void ForR::setForR(double R, Point rightPoint, Point leftPoint)
{
	this->R = R;
	this->rightPoint = rightPoint;
	this->leftPoint = leftPoint;
}

bool operator<(const ForR& a, const ForR& b)
{
	return a.R < b.R;
}

bool operator>(const ForR& a, const ForR& b)
{
	return a.R > b.R;
}
