#include <fstream>
#include "AGMND.h"
#include "Functions.h"
#include <iostream>
#include "DerFunctions.h"


int main()
{
    setlocale(LC_ALL, "rus");
    //Sequential version===========================================================
    {
        std::ofstream out;
        out.open(".\\sequential_agmnd.txt");

        //Constants===========================================================
        double r = 1.5;
        const double eps = 0.01;
        const  size_t maxIterCount = 1000;
        const bool numerical = true;
        //====================================================================

        //Function 1 =========================================================
        double a = 2.7;
        double b = 7.5;
        double (*fnc)(double) = function1;
        double (*fncDer)(double) = derFunction1;

        AGMND test = AGMND(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);

        auto res = test.find_min_value();

        double minX = res.x;
        double minZ = res.z;
        size_t iterCount = res.num_iter;
        double runTime = res.total_runtime;
        double accuracy = res.accuracy;

        out << "Function1:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;

        //====================================================================

        //Function 2 =========================================================
        a = 0.0;
        b = 10.0;
        fnc = function2;
        fncDer = derFunction2;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function2:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 3 =========================================================
        a = 0.0;
        b = 1.2;
        fnc = function3;
        fncDer = derFunction3;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function3:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 4 =========================================================
        a = -10.0;
        b = 10.0;
        fnc = function4;
        fncDer = derFunction4;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function4:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 5 =========================================================
        a = 2.7;
        b = 7.5;
        fnc = function5;
        fncDer = derFunction5;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function5:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 7 =========================================================
        a = 0.0;
        b = 4.0;
        fnc = function7;
        fncDer = derFunction7;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function7:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 8 =========================================================
        a = -5.0;
        b = 5.0;
        fnc = function8;
        fncDer = derFunction8;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function8:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 9 =========================================================
        a = 0.0;
        b = 6.5;
        fnc = function9;
        fncDer = derFunction9;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function9:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 10 ========================================================
        a = -3.0;
        b = 3.0;
        fnc = function10;
        fncDer = derFunction10;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, numerical);
        res = test.find_min_value();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function10:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================
    }

    //Parallel version===========================================================

    {
        std::ofstream out;
        out.open(".\\parallel_agmnd.txt");

        //Constants===========================================================
        double r = 1.5;
        const double eps = 0.01;
        const  size_t maxIterCount = 1000;
        size_t num_threads = 6;
        const bool numerical = true;
        //====================================================================

        //Function 1 =========================================================
        double a = 2.7;
        double b = 7.5;
        double (*fnc)(double) = function1;
        double (*fncDer)(double) = derFunction1;

        AGMND test = AGMND(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);

        auto res = test.find_min_value_parallel();

        double minX = res.x;
        double minZ = res.z;
        size_t iterCount = res.num_iter;
        double runTime = res.total_runtime;
        double accuracy = res.accuracy;

        out << "Function1:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;

        //====================================================================

        //Function 2 =========================================================
        a = 0.0;
        b = 10.0;
        fnc = function2;
        fncDer = derFunction2;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function2:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 3 =========================================================
        a = 0.0;
        b = 1.2;
        fnc = function3;
        fncDer = derFunction3;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function3:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 4 =========================================================
        a = -10.0;
        b = 10.0;
        fnc = function4;
        fncDer = derFunction4;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function4:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 5 =========================================================
        a = 2.7;
        b = 7.5;
        fnc = function5;
        fncDer = derFunction5;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function5:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 7 =========================================================
        a = 0.0;
        b = 4.0;
        fnc = function7;
        fncDer = derFunction7;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function7:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 8 =========================================================
        a = -5.0;
        b = 5.0;
        fnc = function8;
        fncDer = derFunction8;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function8:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 9 =========================================================
        a = 0.0;
        b = 6.5;
        fnc = function9;
        fncDer = derFunction9;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function9:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================

        //Function 10 ========================================================
        a = -3.0;
        b = 3.0;
        fnc = function10;
        fncDer = derFunction10;

        test.reset(a, b, r, eps, maxIterCount, fnc, fncDer, num_threads, numerical);
        res = test.find_min_value_parallel();

        minX = res.x;
        minZ = res.z;
        iterCount = res.num_iter;
        runTime = res.total_runtime;
        accuracy = res.accuracy;

        out << "Function10:" << std::endl;
        out << "MinX: " << minX << std::endl;
        out << "MinZ: " << minZ << std::endl;
        out << "IterCount: " << iterCount << std::endl;
        out << "Accuracy: " << accuracy << std::endl;
        out << "RunTime: " << runTime << std::endl;
        out << std::endl;
        //====================================================================
    }
}
