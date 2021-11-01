#define _USE_MATH_DEFINES
#include "Functions.h"

double function1(const double x)
{
	double y = std::sin(x) + std::sin(10 * x / 3);
	return y;
}

double function2(const double x)
{
	double y = 0;

	for (size_t k = 1; k <= 5; ++k)
	{
		y += (k * std::sin((k + 1) * x + k));
	}

	return (-y);
}

double function3(const double x)
{
	double y = (3 * x - 1.4) * std::sin(18 * x);
	return y;
}

double function4(const double x)
{
	double y = (-1) * (x + std::sin(x)) * std::exp(-(x * x));
	return y;
}

double function5(const double x)
{
	double y = std::sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3;
	return y;
}

double function7(const double x)
{
	double y = -std::sin(2 * M_PI * x) * std::exp(-x);
	return y;
}

double function8(const double x)
{
	double y = (x * x - 5 * x + 6) / (x * x + 1);
	return y;
}

double function9(const double x)
{
	double y = -x + std::sin(3 * x) - 1;
	return y;
}

double function10(const double x)
{
	double y = 2 * (x - 3) * (x - 3) + std::exp((x * x) / 2);
	return y;
}
