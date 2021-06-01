#include <fstream>
#include <iostream>
#include <cstdlib>
#include <random>
#include <vector>
#include "DanPshen.h"

using namespace std;

void SimpleRandomSearch(ofstream& fout, double p, double eps, vector<double>& x, vector<double>& y)
{
   double p_e, Func;
   double Fmin = numeric_limits<double>::infinity();
   vector<double> minPoint(2), point(2);

   p_e = eps * eps / ((y[1] - y[0]) * (x[1] - x[0]));
   int N = ceil(log(1 - p) / log(1 - p_e));

   for (int i = 0; i < N; i++)
   {
      static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
      point[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
      point[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);
      Func = F(point);
      if (Func < Fmin)
      {
         Fmin = Func;
         minPoint = point;
      }
   }
   fout << scientific << eps << "\t" << p << "\t" << N << "\t (" << minPoint[0] << ", " << minPoint[1] << ")\t" << -Fmin << endl;
}

void SimpleRandomSearchTest(vector<double>& p, vector<double>& eps, vector<double>& x, vector<double>& y)
{
   ofstream fout;
   fout.open("SimpleRandomSearch.txt");
   fout << "Eps             P               N            (x*, y*)                       F(x*,y*)" << endl;
   for (size_t i = 0; i < p.size(); i++)
      for (size_t j = 0; j < eps.size(); j++)
         SimpleRandomSearch(fout, p[i], eps[j], x, y);
   fout.close();
}

void GlobalSearch_1(double eps, vector<double>& x, vector<double>& y, vector<int>& m, ofstream& fout, int& num)
{
   mt19937 gen(num);
   uniform_real_distribution<> distr(x[0], x[1]);
   DanPshen dp(2);
   fout << "Algorithm 1" << endl;
   fout << "     m   f_count     (x*, y*)                       F(x*,y*)" << endl;

   vector<double> minPoint(2), point(2);

   for (int i = 0; i < m.size(); i++)
   {
      double Fmin = numeric_limits<double>::infinity(), Func;
      int f_calc_count = 0;
      for (int j = 0; j < m[i]; j++)
      {
         /*static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
         point[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
         point[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);*/
         
         point[0] = distr(gen);
         point[1] = distr(gen);

         dp.FindExtremum(F, FindMinArgGolden, point, eps, eps);

         Func = F(dp.xk1);
         f_calc_count++;
         if (Func < Fmin)
         {
            Fmin = Func;
            minPoint = dp.xk1;
         }
         f_calc_count += dp.f_calc_count;
      }
      fout << setw(6) << m[i] << "\t" << setw(8) << f_calc_count << scientific << "\t (" << minPoint[0] << ", " << minPoint[1] << ")\t" << -Fmin << endl;
   }
}

void GlobalSearch_2(double eps, vector<double>& x, vector<double>& y, vector<int>& m, ofstream& fout, int& num)
{
   mt19937 gen(num);
   uniform_real_distribution<> distr(x[0], x[1]);
   DanPshen dp(2);
   fout << "Algorithm 2" << endl;
   fout << "     m   f_count     (x*, y*)                       F(x*,y*)" << endl;

   vector<double> minPoint(2), point(2);
   double Fmin, Func;

   for (int j = 0; j < m.size(); j++)
   {
      /*static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
      point[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
      point[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);*/

      point[0] = distr(gen);
      point[1] = distr(gen);

      dp.FindExtremum(F, FindMinArgGolden, point, eps, eps);
      minPoint = dp.xk1;
      Fmin = F(dp.xk1);
      int f_calc_count = dp.f_calc_count + 1;
      int i;

      while (true)
      {
         for (i = 0; i < m[j]; i++)
         {
            /*static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
            point[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
            point[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);*/

            point[0] = distr(gen);
            point[1] = distr(gen);

            Func = F(point);
            f_calc_count++;
            if (Func < Fmin)
            {
               Fmin = Func;
               minPoint = point;
               break;
            }
         }

         if (i == m[j]) { break; }

         dp.FindExtremum(F, FindMinArgGolden, minPoint, eps, eps);
         minPoint = dp.xk1;
         Fmin = F(dp.xk1);
         f_calc_count += dp.f_calc_count + 1;
      }
      fout << setw(6) << i << "\t" << setw(8) << f_calc_count << scientific << "\t (" << minPoint[0] << ", " << minPoint[1] << ")\t" << -Fmin << endl;
   }
}

void GlobalSearch_3(double eps, vector<double>& x, vector<double>& y, vector<int>& m, ofstream& fout, int& num)
{
   mt19937 gen(num);
   uniform_real_distribution<> distr(x[0], x[1]);
   DanPshen dp(2);
   fout << "Algorithm 3" << endl;
   fout << "     m   f_count     (x*, y*)                       F(x*,y*)" << endl;

   vector<double> minPoint(2), point(2), direction(2);
   double Fmin, Func;

   for (int j = 0; j < m.size(); j++)
   {
     /* static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
      point[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
      point[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);*/

      point[0] = distr(gen);
      point[1] = distr(gen);

      dp.FindExtremum(F, FindMinArgGolden, point, eps, eps);
      minPoint = dp.xk1;
      point = minPoint;
      Fmin = F(dp.xk1);
      int i, f_calc_count = dp.f_calc_count + 1;
      double delta = 0.5;

      while (true)
      {
         for (i = 0; i < m[j]; i++)
         {
            /*static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
            direction[0] = static_cast<double>(rand() * fraction * (x[1] - x[0] + 1) + x[0]);
            direction[1] = static_cast<double>(rand() * fraction * (y[1] - y[0] + 1) + y[0]);*/

            direction[0] = distr(gen);
            direction[1] = distr(gen);
           
            vector<double> precPoint = point;
            do
            {
               point[0] += delta * (direction[0] / sqrt(direction[0] * direction[0] + direction[1] * direction[1]));
               point[1] += delta * (direction[1] / sqrt(direction[0] * direction[0] + direction[1] * direction[1]));
               f_calc_count += 2;
            } while (point[0] > -10 && point[0] <= 10 && point[1] > -10 && point[1] <= 10 && F(precPoint) <= F(point));

            if (point[0] > -10 && point[0] <= 10 && point[1] > -10 && point[1] <= 10) { break; }
         }

         if (i == m[j]) { break; }

         dp.FindExtremum(F, FindMinArgGolden, point, eps, eps);
         f_calc_count += dp.f_calc_count + 1;
         Func = F(dp.xk1);
         if (Func < Fmin)
         {
            Fmin = Func;
            minPoint = dp.xk1;
         }
         point = minPoint;
      }
      fout << setw(6) << i << "\t" << setw(8) << f_calc_count << scientific << "\t (" << minPoint[0] << ", " << minPoint[1] << ")\t" << -Fmin << endl;
   }
}

int main()
{
   vector<double> x = { -10, 10 };
   vector<double> y = { -10, 10 };
   vector<double> p = { 0.8, 0.9, 0.95, 0.99 };
   vector<double> eps = { 1, 0.1, 0.01, 0.001 };
   vector<int> m = { 1, 10, 50, 100, 500, 1000, 10000 };
   int num = 4957592;

   //SimpleRandomSearchTest(p, eps, x, y);

   ofstream fout;
   fout.open("GlobalSearc.txt");

   //GlobalSearch_1(1e-6, x, y, m, fout, num);
   //GlobalSearch_2(1e-6, x, y, m, fout, num);
   GlobalSearch_3(1e-6, x, y, m, fout, num);

   fout.close();
}