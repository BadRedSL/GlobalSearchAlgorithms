#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "ForR.h"
#include "Functions.h"
#include "Point.h"
#include "evolvent.h"
#include "grishagin_function.hpp"

// z - function value
// z - argument in multidimensional space
// x - argument in one-dimensional space

static const double rand_minimums[] = {
    0.603052, 0.408337, /*f(min1)=-13.51436*/
    0.652988, 0.320592, /*f(min2)=-11.28447*/
    1.000000, 0.000000, /*f(min3)=-13.20907*/
    0.066182, 0.582587, /*f(min4)=-11.54117*/
    0.904308, 0.872639, /*f(min5)=-9.969261*/
    0.344375, 0.524932, /*f(min6)=-9.180137*/
    0.000000, 1.000000, /*f(min7)=-9.359486*/
    0.948275, 0.887031, /*f(min8)=-11.46999*/
    0.226047, 0.520153, /*f(min9)=-11.41470*/
    0.341732, 0.197620, /*f(min10)=-12.35783*/
    0.069264, 0.430955, /*f(min11)=-8.298100*/
    0.000000, 1.000000, /*f(min12)=-10.77891*/
    0.45221,  0.07292,  /*f(min13)=-9.639918*/
    0.579769, 0.046396, /*f(min14)=-10.84688*/
    0.000000, 1.000000, /*f(min15)=-9.998392*/
    0.310179, 1.000000, /*f(min16)=-13.27447*/
    0.909758, 0.926195, /*f(min17)=-11.20579*/
    0.434562, 0.825608, /*f(min18)=-8.555728*/
    0.06686,  0.77051,  /*f(min19)=-11.08050*/
    0.641337, 0.135186, /*f(min20)=-10.84137*/
    0.885029, 0.390289, /*f(min21)=-10.09915*/
    0.649650, 0.414282, /*f(min22)=-8.821688*/
    0.142623, 0.157327, /*f(min23)=-11.20434*/
    0.862953, 1.000000, /*f(min24)=-9.774841*/
    0.46036,  0.99314,  /*f(min25)=-9.446216*/
    0.379189, 0.688051, /*f(min26)=-9.922234*/
    0.845292, 0.424546, /*f(min27)=-9.353030*/
    0.441160, 0.016803, /*f(min28)=-8.927842*/
    1.000000, 1.000000, /*f(min29)=-11.97038*/
    0.303295, 0.134722, /*f(min30)=-11.05922*/
    0.109520, 0.265486, /*f(min31)=-8.961329*/
    1.000000, 0.000000, /*f(min32)=-10.25347*/
    0.593726, 0.503014, /*f(min33)=-13.75610*/
    0.694905, 1.000000, /*f(min34)=-10.07243*/
    0.051975, 0.409344, /*f(min35)=-10.66758*/
    0.125664, 0.518969, /*f(min36)=-9.606340*/
    0.000000, 0.000000, /*f(min37)=-11.10867*/
    0.155081, 0.238663, /*f(min38)=-8.586483*/
    0.53707,  0.46181,  /*f(min39)=-8.823492*/
    0.110985, 0.917791, /*f(min40)=-10.22533*/
    1.000000, 0.000000, /*f(min41)=-13.84155*/
    0.776095, 0.764724, /*f(min42)=-10.76901*/
    0.087367, 0.677632, /*f(min43)=-8.574448*/
    0.308037, 0.536113, /*f(min44)=-10.40137*/
    0.042100, 0.563607, /*f(min45)=-8.889051*/
    0.287025, 0.159219, /*f(min46)=-10.44960*/
    0.451926, 0.169839, /*f(min47)=-10.45448*/
    0.884761, 0.245341, /*f(min48)=-9.749494*/
    0.047782, 0.171633, /*f(min49)=-10.80496*/
    0.00000,  0.41596,  /*f(min50)=-12.16739*/
    0.192108, 0.303789, /*f(min51)=-11.14192*/
    0.554153, 0.809821, /*f(min52)=-10.06221*/
    0.91475,  0.54149,  /*f(min53)=-9.518639*/
    0.661592, 0.925902, /*f(min54)=-9.404589*/
    0.962924, 0.436680, /*f(min55)=-10.19342*/
    0.000000, 0.000000, /*f(min56)=-9.009641*/
    0.616058, 0.560244, /*f(min57)=-10.02504*/
    0.439890, 0.343722, /*f(min58)=-10.95753*/
    0.218146, 0.677192, /*f(min59)=-9.850361*/
    1.000000, 1.000000, /*f(min60)=-11.74782*/
    0.198145, 0.317876, /*f(min61)=-10.73907*/
    0.875874, 0.653336, /*f(min62)=-9.060382*/
    0.22999,  0.33624,  /*f(min63)=-10.74852*/
    0.169351, 0.015656, /*f(min64)=-8.455091*/
    0.760073, 0.906035, /*f(min65)=-11.29555*/
    0.702941, 0.308403, /*f(min66)=-10.47617*/
    0.365371, 0.282325, /*f(min67)=-9.285640*/
    0.314012, 0.651377, /*f(min68)=-9.592745*/
    0.237687, 0.374368, /*f(min69)=-10.82136*/
    0.586334, 0.508672, /*f(min70)=-9.353659*/
    0.000000, 0.000000, /*f(min71)=-9.803835*/
    0.383319, 1.000000, /*f(min72)=-10.86034*/
    0.780103, 0.103783, /*f(min73)=-10.98085*/
    0.350265, 0.566946, /*f(min74)=-8.966131*/
    0.798535, 0.478706, /*f(min75)=-11.00382*/
    0.31759,  0.06967,  /*f(min76)=-9.959365*/
    0.715929, 0.704778, /*f(min77)=-11.38632*/
    0.563040, 0.442557, /*f(min78)=-11.44656*/
    0.565078, 0.322618, /*f(min79)=-14.17049*/
    0.146731, 0.510509, /*f(min80)=-10.61985*/
    0.000000, 0.543167, /*f(min81)=-11.31970*/
    0.208533, 0.454252, /*f(min82)=-10.95230*/
    0.155111, 0.972329, /*f(min83)=-11.41650*/
    0.000000, 1.000000, /*f(min84)=-12.05415*/
    0.336467, 0.909056, /*f(min85)=-10.45893*/
    0.57001,  0.90847,  /*f(min86)=-9.429290*/
    0.296290, 0.540579, /*f(min87)=-10.45261*/
    0.172262, 0.332732, /*f(min88)=-10.75306*/
    0.000000, 1.000000, /*f(min89)=-9.816527*/
    1.000000, 0.000000, /*f(min90)=-12.20688*/
    1.000000, 1.000000, /*f(min91)=-10.39147*/
    0.674061, 0.869954, /*f(min92)=-9.689030*/
    1.000000, 1.000000, /*f(min93)=-11.91809*/
    0.852506, 0.637278, /*f(min94)=-10.42941*/
    0.877491, 0.399780, /*f(min95)=-10.34945*/
    0.835605, 0.751888, /*f(min96)=-9.013732*/
    0.673378, 0.827427, /*f(min97)=-8.916823*/
    0.831754, 0.367117, /*f(min98)=-8.747803*/
    0.601971, 0.734465, /*f(min99)=-12.76981*/
    0.000000, 0.000000  /*f(min100)=-11.43793*/
};

static const double rand_minimums_values[] = {
    -13.51436,  /*f(min1)=-13.51436*/
    -11.28447,  /*f(min2)=-11.28447*/
    -13.20907,  /*f(min3)=-13.20907*/
    -11.54117,  /*f(min4)=-11.54117*/
    -9.969261,  /*f(min5)=-9.969261*/
    -9.180137,  /*f(min6)=-9.180137*/
    -9.359486,  /*f(min7)=-9.359486*/
    -11.46999,  /*f(min8)=-11.46999*/
    -11.41470,  /*f(min9)=-11.41470*/
    -12.35783,  /*f(min10)=-12.35783*/
    -8.298100,  /*f(min11)=-8.298100*/
    -10.77891,  /*f(min12)=-10.77891*/
    -9.639918,  /*f(min13)=-9.639918*/
    -10.84688,  /*f(min14)=-10.84688*/
    -9.998392,  /*f(min15)=-9.998392*/
    -13.27447,  /*f(min16)=-13.27447*/
    -11.20579,  /*f(min17)=-11.20579*/
    -8.555728,  /*f(min18)=-8.555728*/
    -11.08050,  /*f(min19)=-11.08050*/
    -10.84137,  /*f(min20)=-10.84137*/
    -10.09915,  /*f(min21)=-10.09915*/
    -8.821688,  /*f(min22)=-8.821688*/
    -11.20434,  /*f(min23)=-11.20434*/
    -9.774841,  /*f(min24)=-9.774841*/
    -9.446216,  /*f(min25)=-9.446216*/
    -9.922234,  /*f(min26)=-9.922234*/
    -9.353030,  /*f(min27)=-9.353030*/
    -8.927842,  /*f(min28)=-8.927842*/
    -11.97038,  /*f(min29)=-11.97038*/
    -11.05922,  /*f(min30)=-11.05922*/
    -8.961329,  /*f(min31)=-8.961329*/
    -10.25347,  /*f(min32)=-10.25347*/
    -13.75610,  /*f(min33)=-13.75610*/
    -10.07243,  /*f(min34)=-10.07243*/
    -10.66758,  /*f(min35)=-10.66758*/
    -9.606340,  /*f(min36)=-9.606340*/
    -11.10867,  /*f(min37)=-11.10867*/
    -8.586483,  /*f(min38)=-8.586483*/
    -8.823492,  /*f(min39)=-8.823492*/
    -10.22533,  /*f(min40)=-10.22533*/
    -13.84155,  /*f(min41)=-13.84155*/
    -10.76901,  /*f(min42)=-10.76901*/
    -8.574448,  /*f(min43)=-8.574448*/
    -10.40137,  /*f(min44)=-10.40137*/
    -8.889051,  /*f(min45)=-8.889051*/
    -10.44960,  /*f(min46)=-10.44960*/
    -10.45448,  /*f(min47)=-10.45448*/
    -9.749494,  /*f(min48)=-9.749494*/
    -10.80496,  /*f(min49)=-10.80496*/
    -12.16739,  /*f(min50)=-12.16739*/
    -11.14192,  /*f(min51)=-11.14192*/
    -10.06221,  /*f(min52)=-10.06221*/
    -9.518639,  /*f(min53)=-9.518639*/
    -9.404589,  /*f(min54)=-9.404589*/
    -10.19342,  /*f(min55)=-10.19342*/
    -9.009641,  /*f(min56)=-9.009641*/
    -10.02504,  /*f(min57)=-10.02504*/
    -10.95753,  /*f(min58)=-10.95753*/
    -9.850361,  /*f(min59)=-9.850361*/
    -11.74782,  /*f(min60)=-11.74782*/
    -10.73907,  /*f(min61)=-10.73907*/
    -9.060382,  /*f(min62)=-9.060382*/
    -10.74852,  /*f(min63)=-10.74852*/
    -8.455091,  /*f(min64)=-8.455091*/
    -11.29555,  /*f(min65)=-11.29555*/
    -10.47617,  /*f(min66)=-10.47617*/
    -9.285640,  /*f(min67)=-9.285640*/
    -9.592745,  /*f(min68)=-9.592745*/
    -10.82136,  /*f(min69)=-10.82136*/
    -9.353659,  /*f(min70)=-9.353659*/
    -9.803835,  /*f(min71)=-9.803835*/
    -10.86034,  /*f(min72)=-10.86034*/
    -10.98085,  /*f(min73)=-10.98085*/
    -8.966131,  /*f(min74)=-8.966131*/
    -11.00382,  /*f(min75)=-11.00382*/
    -9.959365,  /*f(min76)=-9.959365*/
    -11.38632,  /*f(min77)=-11.38632*/
    -11.44656,  /*f(min78)=-11.44656*/
    -14.17049,  /*f(min79)=-14.17049*/
    -10.61985,  /*f(min80)=-10.61985*/
    -11.31970,  /*f(min81)=-11.31970*/
    -10.95230,  /*f(min82)=-10.95230*/
    -11.41650,  /*f(min83)=-11.41650*/
    -12.05415,  /*f(min84)=-12.05415*/
    -10.45893,  /*f(min85)=-10.45893*/
    -9.429290,  /*f(min86)=-9.429290*/
    -10.45261,  /*f(min87)=-10.45261*/
    -10.75306,  /*f(min88)=-10.75306*/
    -9.816527,  /*f(min89)=-9.816527*/
    -12.20688,  /*f(min90)=-12.20688*/
    -10.39147,  /*f(min91)=-10.39147*/
    -9.6890304, /*f(min92)=-9.689030*/
    -11.91809,  /*f(min93)=-11.91809*/
    -10.42941,  /*f(min94)=-10.42941*/
    -10.34945,  /*f(min95)=-10.34945*/
    -9.013732,  /*f(min96)=-9.013732*/
    -8.916823,  /*f(min97)=-8.916823*/
    -8.747803,  /*f(min98)=-8.747803*/
    -12.76981,  /*f(min99)=-12.76981*/
    -11.43793   /*f(min100)=-11.43793*/
};

vagrish::GrishaginFunction grishaginFunc;

double GrishaginFunction(double* args, const size_t len) {
  return grishaginFunc.Calculate(args);
}

void newPointCalculator(const Point leftPoint, const Point rightPoint,
                        double (*fnc)(double*, size_t), const size_t N,
                        const double m, const double* a, const double* b,
                        Point& newPoint) {
  TEvolvent evolvent;
  evolvent.SetBounds(a, b);
  double newX = ((rightPoint.x + leftPoint.x) / 2) -
                ((rightPoint.z - leftPoint.z) / (2 * m));
  double* newY = new double[N];
  evolvent.GetImage(newX, newY);

  double newZ = fnc(newY, N);

  newPoint.setPoint(newX, newZ);

  delete[] newY;
}

void newMCalculator(const Point leftPoint, const Point rightPoint, double& M) {
  M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));
}

void newRCalculator(const Point leftPoint, const Point rightPoint,
                    const double m, double& R) {
  R = (m * (rightPoint.x - leftPoint.x)) +
      (((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z)) /
       (m * rightPoint.x - leftPoint.x)) -
      (2 * rightPoint.z + leftPoint.z);
}

void allParallelCalculator(const Point leftPoint, const Point rightPoint,
                           const double r, double (*fnc)(double*, size_t),
                           const size_t N, bool recountM, double maxM,
                           const double* a, const double* b, Point& newPoint,
                           double& M1, double& M2, double& R1, double& R2,
                           ForR& ForR1, ForR& ForR2, double& vm) {
  newPointCalculator(leftPoint, rightPoint, fnc, N, vm, a, b, newPoint);
  newMCalculator(leftPoint, newPoint, M1);
  newMCalculator(newPoint, rightPoint, M2);
  double M;
  if (recountM) {
    M = (M1 > M2) ? M1 : M2;
  } else {
    M = (M1 > M2) ? M1 : M2;
    M = (maxM > M) ? maxM : M;
  }
  vm = (M > 0) ? r * M : 1;
  newRCalculator(leftPoint, newPoint, vm, R1);
  ForR1.setForR(R1, newPoint, leftPoint);
  newRCalculator(newPoint, rightPoint, vm, R2);
  ForR2.setForR(R2, rightPoint, newPoint);
}

void parallelMultidimensionalMinValue(size_t numThread, double* a, double* b,
                                      size_t N, double r, double eps,
                                      size_t nMax,
                                      double (*fnc)(double*, size_t),
                                      double* resultY, double& resultF,
                                      std::chrono::duration<double>& runTime,
                                      size_t& iterCount, double& accuracy) {
  std::vector<std::thread> vThread;
  std::vector<Point> vNewPoint(numThread);
  std::vector<double> vM1(numThread);
  std::vector<double> vM2(numThread);
  std::vector<double> vR1(numThread);
  std::vector<double> vR2(numThread);
  std::vector<ForR> vForR1(numThread);
  std::vector<ForR> vForR2(numThread);
  std::vector<double> vm(numThread);

  auto start = std::chrono::system_clock::now();
  TEvolvent evolvent;
  evolvent.SetBounds(a, b);

  bool changeM = false;
  bool exitFlag = false;
  std::map<Point, double> set;
  auto* rQueue =
      new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

  double minX;
  double* newY = new double[N];
  double* y1 = new double[N];
  double* y2 = new double[N];
  evolvent.GetImage(0, y1);
  evolvent.GetImage(1, y2);

  (fnc(y1, N) < fnc(y2, N)) ? minX = 0 : minX = 1;

  double min;
  (fnc(y1, N) < fnc(y2, N)) ? min = fnc(y1, N) : min = fnc(y2, N);

  Point rightPoint(1, fnc(y2, N));
  Point leftPoint(0, fnc(y1, N));
  Point newPoint(1, fnc(y2, N));

  set.emplace(leftPoint, 0);
  set.emplace(rightPoint, 0);

  double M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));

  double m = (M > 0) ? r * M : 1;

  double R = (m * (rightPoint.x - leftPoint.x)) +
             ((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z) /
              (m * (rightPoint.x - leftPoint.x))) -
             2 * (rightPoint.z + leftPoint.z);

  ForR tmpR(R, rightPoint, leftPoint);
  (*rQueue).push(tmpR);

  size_t count = 0;

  while (count < nMax) {
    if (numThread >= (*rQueue).size()) {
      rightPoint = (*rQueue).top().rightPoint;
      leftPoint = (*rQueue).top().leftPoint;

      if (rightPoint.x - leftPoint.x <= eps) {
        break;
      }

      newPointCalculator(leftPoint, rightPoint, fnc, N, m, a, b, newPoint);
      set.emplace(newPoint, 0);

      double tmpM = 0;
      double prevM = M;

      newMCalculator(leftPoint, rightPoint, tmpM);
      if (M == tmpM) {
        M = -1;
        for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end();
             ++cur, ++prev) {
          newMCalculator((*prev).first, (*cur).first, tmpM);
          M = (tmpM > M) ? tmpM : M;
        }
      } else {
        newMCalculator(leftPoint, newPoint, tmpM);
        M = (tmpM > M) ? tmpM : M;
        newMCalculator(newPoint, rightPoint, tmpM);
        M = (tmpM > M) ? tmpM : M;
      }

      changeM = (M == prevM) ? false : true;
      m = (M > 0) ? r * M : 1;
      for (auto& iter : vm) {
        iter = m;
      }

      if (!changeM) {
        (*rQueue).pop();
        newRCalculator(leftPoint, newPoint, m, R);
        tmpR.setForR(R, newPoint, leftPoint);
        (*rQueue).push(tmpR);
        newRCalculator(newPoint, rightPoint, m, R);
        tmpR.setForR(R, rightPoint, newPoint);
        (*rQueue).push(tmpR);
      } else {
        delete rQueue;
        rQueue =
            new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

        for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end();
             ++cur, ++prev) {
          newRCalculator((*prev).first, (*cur).first, m, R);
          tmpR.setForR(R, (*cur).first, (*prev).first);
          (*rQueue).push(tmpR);
        }
      }
      ++count;
      if (newPoint.z < min) {
        min = newPoint.z;
        minX = newPoint.x;
        evolvent.GetImage(minX, newY);
        for (size_t i = 0; i < N; ++i) {
          resultY[i] = newY[i];
        }
      }
    } else {
      double prevM = M;
      double tmpM = 0;
      bool recountM = false;
      for (size_t i = 0; i < numThread; ++i) {
        rightPoint = (*rQueue).top().rightPoint;
        leftPoint = (*rQueue).top().leftPoint;

        newMCalculator(leftPoint, rightPoint, tmpM);

        if (M == tmpM) {
          recountM = true;
        }

        if (rightPoint.x - leftPoint.x <= eps) {
          exitFlag = true;
        }
        vThread.push_back(std::thread(
            allParallelCalculator, leftPoint, rightPoint, r, fnc, N, recountM,
            M, a, b, std::ref(vNewPoint[i]), std::ref(vM1[i]), std::ref(vM2[i]),
            std::ref(vR1[i]), std::ref(vR2[i]), std::ref(vForR1[i]),
            std::ref(vForR2[i]), std::ref(vm[i])));
        (*rQueue).pop();
      }
      for (size_t i = 0; i < numThread; ++i) {
        vThread[i].join();
      }
      vThread.clear();
      if (exitFlag) {
        break;
      }

      if (recountM) {
        for (size_t i = 0; i < numThread; ++i) {
          set.emplace(vNewPoint[i], 0);
          (*rQueue).push(vForR1[i]);
          (*rQueue).push(vForR2[i]);
        }
        M = -1;
        for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end();
             ++cur, ++prev) {
          newMCalculator((*prev).first, (*cur).first, tmpM);
          M = (tmpM > M) ? tmpM : M;
        }
        recountM = false;
      } else {
        for (size_t i = 0; i < numThread; ++i) {
          set.emplace(vNewPoint[i], 0);
          M = (vM1[i] > M) ? vM1[i] : M;
          M = (vM2[i] > M) ? vM2[i] : M;
          (*rQueue).push(vForR1[i]);
          (*rQueue).push(vForR2[i]);
        }
      }

      for (auto& iter : vm) {
        iter = (M > 0) ? r * M : 1;
      }
      if (prevM != M) {
        delete rQueue;
        rQueue =
            new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

        for (auto cur = (++set.begin()), prev = set.begin(); cur != set.end();
             ++cur, ++prev) {
          newRCalculator((*prev).first, (*cur).first, vm[0], R);
          tmpR.setForR(R, (*cur).first, (*prev).first);
          (*rQueue).push(tmpR);
        }
      }
      count += numThread;
      for (const auto iter : vNewPoint) {
        if (iter.z < min) {
          min = iter.z;
          minX = iter.x;
          evolvent.GetImage(minX, newY);
          for (size_t i = 0; i < N; ++i) {
            resultY[i] = newY[i];
          }
        }
      }
    }
  }
  accuracy = rightPoint.x - leftPoint.x;
  iterCount = count;
  resultF = min;

  delete[] y1;
  delete[] y2;
  delete[] newY;

  auto end = std::chrono::system_clock::now();
  runTime = end - start;
}

void multidimensionalMinValue(double* a, double* b, size_t N, double r,
                              double eps, size_t nMax,
                              double (*fnc)(double*, size_t), double* resultY,
                              double& resultF,
                              std::chrono::duration<double>& runTime,
                              size_t& iterCount, double& accuracy) {
  auto start = std::chrono::system_clock::now();
  TEvolvent evolvent;
  evolvent.SetBounds(a, b);

  std::map<Point, double> set;
  auto* rQueue =
      new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;
  std::priority_queue<double, std::vector<double>, std::less<double>> mQueue;

  double minX;
  double* newY = new double[N];
  double* y1 = new double[N];
  double* y2 = new double[N];
  evolvent.GetImage(0, y1);
  evolvent.GetImage(1, y2);

  (fnc(y1, N) < fnc(y2, N)) ? minX = 0 : minX = 1;

  double min;
  (fnc(y1, N) < fnc(y2, N)) ? min = fnc(y1, N) : min = fnc(y2, N);

  Point rightPoint(1, fnc(y2, N));
  Point leftPoint(0, fnc(y1, N));
  Point newPoint(1, fnc(y2, N));
  set.emplace(leftPoint, 0);
  auto curIter = set.emplace(rightPoint, 0);

  double M = abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x));
  double prevM = M;

  mQueue.push(M);

  double m = (M > 0) ? r * M : 1;

  double R = (m * (rightPoint.x - leftPoint.x)) +
             ((rightPoint.z - leftPoint.z) * (rightPoint.z - leftPoint.z) /
              (m * (rightPoint.x - leftPoint.x))) -
             2 * (rightPoint.z + leftPoint.z);

  ForR newR(R, rightPoint, leftPoint);
  (*rQueue).push(newR);

  size_t count = 0;

  while ((count < nMax) && (rightPoint.x - leftPoint.x > eps)) {
    double newX = ((rightPoint.x + leftPoint.x) / 2) -
                  ((rightPoint.z - leftPoint.z) / (2 * m));

    newY = new double[N];
    evolvent.GetImage(newX, newY);
    double newZ = fnc(newY, N);  // newZ

    newPoint.setPoint(newX, newZ);

    if (newZ < min) {
      minX = newX;
      min = newZ;
      for (size_t i = 0; i < N; ++i) {
        resultY[i] = newY[i];
      }
    }

    curIter = set.emplace(newPoint, 0);

    auto curTmp = curIter.first;
    auto prevTmp = std::prev(curTmp);

    if (mQueue.top() ==
        abs((rightPoint.z - leftPoint.z) / (rightPoint.x - leftPoint.x))) {
      mQueue.pop();
    }

    for (size_t i = 0; i < 2; ++i, ++curTmp, ++prevTmp) {
      M = abs(((*curTmp).first.z - (*prevTmp).first.z) /
              ((*curTmp).first.x - (*prevTmp).first.x));

      mQueue.push(M);
    }

    if (mQueue.top() > 0) {
      m = r * mQueue.top();
    } else {
      m = 1;
    }

    if (mQueue.top() == prevM) {
      curTmp = curIter.first;
      prevTmp = std::prev(curTmp);

      for (size_t i = 0; i < 2; ++i, ++curTmp, ++prevTmp) {
        R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
            ((((*curTmp).first.z - (*prevTmp).first.z) *
              ((*curTmp).first.z - (*prevTmp).first.z)) /
             (m * (*curTmp).first.x - (*prevTmp).first.x)) -
            (2 * (*curTmp).first.z + (*prevTmp).first.z);

        newR.setForR(R, (*curTmp).first, (*prevTmp).first);
        if (i == 0) {
          (*rQueue).pop();
        }
        (*rQueue).push(newR);
      }
    } else {
      delete rQueue;
      rQueue =
          new std::priority_queue<ForR, std::vector<ForR>, std::less<ForR>>;

      for (curTmp = (++set.begin()), prevTmp = set.begin(); curTmp != set.end();
           ++curTmp, ++prevTmp) {
        R = (m * ((*curTmp).first.x - (*prevTmp).first.x)) +
            ((((*curTmp).first.z - (*prevTmp).first.z) *
              ((*curTmp).first.z - (*prevTmp).first.z)) /
             (m * (*curTmp).first.x - (*prevTmp).first.x)) -
            (2 * (*curTmp).first.z + (*prevTmp).first.z);

        newR.setForR(R, (*curTmp).first, (*prevTmp).first);
        (*rQueue).push(newR);
      }
    }
    rightPoint = (*rQueue).top().rightPoint;
    leftPoint = (*rQueue).top().leftPoint;
    prevM = mQueue.top();
    ++count;
  }

  accuracy = rightPoint.x - leftPoint.x;
  iterCount = count;
  resultF = min;
  auto end = std::chrono::system_clock::now();
  runTime = end - start;
}

int main() {
  setlocale(LC_ALL, "ru");

  // double a[2] = {0.0, 0.0};
  // double b[2] = {8.0, 8.0};
  // size_t N = 2;
  // double r = 2.0;
  // double eps = 0.001;
  // size_t nMax = 10000;
  // double resultY[2] = {0.0, 0.0};
  // double resultF = 0.0;
  // std::chrono::duration<double> runTime;
  // size_t iterCount = 0;
  // double accuracy = 0.0;

  // multidimensionalMinValue(a, b, N, r, eps, nMax, twoDimentionFunction1,
  //                         resultY, resultF, runTime, iterCount, accuracy);
  // std::cout << "Значение функции = " << resultF << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultY[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCount << std::endl;
  // std::cout << "Точность алгоритма = " << accuracy << std::endl;
  // std::cout << "Время работы алгоритма = " << runTime.count() << std::endl;
  // std::cout << std::endl;

  // grishaginFunc.SetFunctionNumber(3);
  // double a1[2];
  // double b1[2];
  // grishaginFunc.GetBounds(a1, b1);
  // size_t N1 = grishaginFunc.GetDimension();
  // double r1 = 2.0;
  // double eps1 = 0.001;
  // size_t nMax1 = 10000;
  // double resultY1[2] = {0.0, 0.0};
  // double resultF1 = 0.0;
  // std::chrono::duration<double> runTime1;
  // size_t iterCount1 = 0;
  // double accuracy1 = 0.0;

  // multidimensionalMinValue(a1, b1, N1, r1, eps1, nMax1, GrishaginFunction,
  //                         resultY1, resultF1, runTime1, iterCount1,
  //                         accuracy1);
  // std::cout << "Значение функции = " << resultF1 << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultY1[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCount1 << std::endl;
  // std::cout << "Точность алгоритма = " << accuracy1 << std::endl;
  // std::cout << "Время работы алгоритма = " << runTime1.count() << std::endl;
  // std::cout << std::endl;

  // grishaginFunc.SetFunctionNumber(100);
  // double a100[2];
  // double b100[2];
  // grishaginFunc.GetBounds(a100, b100);
  // size_t N100 = grishaginFunc.GetDimension();
  // double r100 = 2.0;
  // double eps100 = 0.001;
  // size_t nMax100 = 10000;
  // double resultY100[2] = {0.0, 0.0};
  // double resultF100 = 0.0;
  // std::chrono::duration<double> runTime100;
  // size_t iterCount100 = 0;
  // double accuracy100 = 0.0;

  // multidimensionalMinValue(a100, b100, N100, r100, eps100, nMax100,
  //                         GrishaginFunction, resultY100, resultF100,
  //                         runTime100, iterCount100, accuracy100);
  // std::cout << "Значение функции = " << resultF100 << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultY100[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCount100 << std::endl;
  // std::cout << "Точность алгоритма = " << accuracy100 << std::endl;
  // std::cout << "Время работы алгоритма = " << runTime100.count() <<
  // std::endl; std::cout << std::endl;

  ////===========================================================================================================================================
  // std::cout << "Параллельные версии:" << std::endl << std::endl;
  ////===========================================================================================================================================

  // size_t numThreads = 12;

  // double ap[2] = {0.0, 0.0};
  // double bp[2] = {8.0, 8.0};
  // size_t Np = 2;
  // double rp = 2.0;
  // double epsp = 0.001;
  // size_t nMaxp = 10000;
  // double resultYp[2] = {0.0, 0.0};
  // double resultFp = 0.0;
  // std::chrono::duration<double> runTimep;
  // size_t iterCountp = 0;
  // double accuracyp = 0.0;

  // parallelMultidimensionalMinValue(numThreads, ap, bp, Np, rp, epsp, nMaxp,
  //                                 twoDimentionFunction1, resultYp, resultFp,
  //                                 runTimep, iterCountp, accuracyp);
  // std::cout << "Значение функции = " << resultFp << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultYp[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCountp << std::endl;
  // std::cout << "Точность алгоритма = " << accuracyp << std::endl;
  // std::cout << "Время работы алгоритма = " << runTimep.count() << std::endl;
  // std::cout << std::endl;

  // grishaginFunc.SetFunctionNumber(3);
  // double a1p[2];
  // double b1p[2];
  // grishaginFunc.GetBounds(a1p, b1p);
  // size_t N1p = grishaginFunc.GetDimension();
  // double r1p = 2.0;
  // double eps1p = 0.001;
  // size_t nMax1p = 10000;
  // double resultY1p[2] = {0.0, 0.0};
  // double resultF1p = 0.0;
  // std::chrono::duration<double> runTime1p;
  // size_t iterCount1p = 0;
  // double accuracy1p = 0.0;

  // parallelMultidimensionalMinValue(
  //    numThreads, a1p, b1p, N1p, r1p, eps1p, nMax1p, GrishaginFunction,
  //    resultY1p, resultF1p, runTime1p, iterCount1p, accuracy1p);
  // std::cout << "Значение функции = " << resultF1p << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultY1p[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCount1p << std::endl;
  // std::cout << "Точность алгоритма = " << accuracy1p << std::endl;
  // std::cout << "Время работы алгоритма = " << runTime1p.count() << std::endl;
  // std::cout << std::endl;

  // grishaginFunc.SetFunctionNumber(100);
  // double a100p[2];
  // double b100p[2];
  // grishaginFunc.GetBounds(a100p, b100p);
  // size_t N100p = grishaginFunc.GetDimension();
  // double r100p = 2.0;
  // double eps100p = 0.001;
  // size_t nMax100p = 10000;
  // double resultY100p[2] = {0.0, 0.0};
  // double resultF100p = 0.0;
  // std::chrono::duration<double> runTime100p;
  // size_t iterCount100p = 0;
  // double accuracy100p = 0.0;

  // parallelMultidimensionalMinValue(numThreads, a100p, b100p, N100p, r100p,
  //                                 eps100p, nMax100p, GrishaginFunction,
  //                                 resultY100p, resultF100p, runTime100p,
  //                                 iterCount100p, accuracy100p);
  // std::cout << "Значение функции = " << resultF100p << std::endl;
  // std::cout << "Точка минимума: " << std::endl;
  // for (int i = 0; i < N; i++) {
  //  std::cout << resultY100p[i] << " ";
  //}
  // std::cout << std::endl;
  // std::cout << "Число итераций алгоритма = " << iterCount100p << std::endl;
  // std::cout << "Точность алгоритма = " << accuracy100p << std::endl;
  // std::cout << "Время работы алгоритма = " << runTime100p.count() <<
  // std::endl; std::cout << std::endl;

  //===========================================================================================================================================
  std::cout << "Проверка алгоритма по функциям Гришагина:" << std::endl
            << std::endl;
  //===========================================================================================================================================

  const int numThreadsTest = 12;

  std::ofstream out;                   // поток для записи
  out.open("..\\grishagin_test_seq.txt");  // окрываем файл для записи

  for (size_t i = 0; i < 100; ++i) {
    double aTest[2];
    double bTest[2];
    size_t NTest = 2;
    double rTest = 2.0;
    double epsTest = 0.00001;
    size_t nMaxTest = 10000;
    double resultYTest[2] = {0.0, 0.0};
    double resultFTest = 0.0;
    std::chrono::duration<double> runTimeTest;
    size_t iterCountTest = 0;
    double accuracyTest = 0.0;
    grishaginFunc.SetFunctionNumber(i + 1);
    grishaginFunc.GetBounds(aTest, bTest);
    NTest = grishaginFunc.GetDimension();

    //parallelMultidimensionalMinValue(numThreadsTest, aTest, bTest, NTest, rTest,
    //                                 epsTest, nMaxTest, GrishaginFunction,
    //                                 resultYTest, resultFTest, runTimeTest,
    //                                 iterCountTest, accuracyTest);

     multidimensionalMinValue(aTest, bTest, NTest, rTest,
                             epsTest, nMaxTest, GrishaginFunction,
                             resultYTest, resultFTest, runTimeTest,
                             iterCountTest, accuracyTest);

    std::string flag = "\t0";
    if ((abs(rand_minimums[2 * i] - resultYTest[0]) > 0.01) ||
        (abs(rand_minimums[2 * i + 1] - resultYTest[1] > 0.01))) {
      flag = "\t-1";
    } else {
      flag = "\t1";
    }

    // out << std::fixed << std::setprecision(8) << "Grichagin function: " << i
    // + 1
    //    << "\t\tTrue coordinates: " << rand_minimums[2 * i] << " "
    //    << rand_minimums[2 * i + 1]
    //    << "\t\tComputed coordinates: " << resultYTest[0] << " " <<
    //    resultYTest[1]
    //    << "\tTrue values: " << resultFTest
    //    << "\tComputed values: " << rand_minimums_values[i] << flag <<
    //    std::endl;

    out << std::fixed << std::setprecision(8) << "Grichagin function: " << i + 1
        << "\t\tTrue coordinates: " << rand_minimums[2 * i] << " "
        << rand_minimums[2 * i + 1]
        << "\t\tComputed coordinates: " << resultYTest[0] << " "
        << resultYTest[1] << "\tTrue values: " << resultFTest
        << "\tComputed values: " << rand_minimums_values[i]
        << "\tRun time: " << runTimeTest.count()
        << "\tNumber of iterations: " << iterCountTest << std::endl;
  }

  return 0;
}
