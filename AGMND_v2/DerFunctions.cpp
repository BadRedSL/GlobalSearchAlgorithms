#define _USE_MATH_DEFINES
#include "DerFunctions.h"

double derFunction1(const double x)
{
	double y = std::cos(x) + (10/3) * std::cos((10/3)*x);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction2(const double x)
{
	double y = 0;

	for (size_t k = 1; k <= 5; ++k)
	{
		y += (k * std::cos((k + 1) * x + k)*(k+1));
	}

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction3(const double x)
{
	double y = 18*(3*x-1.4)*std::cos(18*x)+3*std::sin(18*x);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction4(const double x)
{
	double y = 2*x*(x+std::sin(x))*std::exp(-(x*x))-(std::cos(x)+1)*std::exp(-(x*x));

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction5(const double x)
{
	double y = -21 / 25 + 1 / x + 10 * std::cos((10 * x) / 3) / 3 + std::cos(x);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction7(const double x)
{
	double y = exp(-x) * sin((2 * M_PI) * x) - 2 * M_PI * cos((2 * M_PI) * x) * exp(-x);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction8(const double x)
{
	double y = ((-5 + 2 * x) / (x * x + 1) - (2 * x * (x * x - 5 * x + 6))) / ((x * x + 1) * (x * x + 1));

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction9(const double x)
{
	double y = 3 * std::cos(3 * x);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction10(const double x)
{
	double y = -6 + 2 * x + 2 * (x - 3) + x * std::exp((x * x) / 2);

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

double derFunction11(const double x)
{
	double y = 2*x;

	for (size_t i = 0; i < 100000; ++i)
	{
		auto tmp = std::chrono::system_clock::now();
	}

	return y;
}

