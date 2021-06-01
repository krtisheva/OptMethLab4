#pragma once
#include <vector>
#include "Vector.h"
using namespace std;

const double SQRT5 = sqrt(5);
const double PI = 3.14159265;

double F(const vector<double>& point)
{
   double sum = 0.;
   vector<int> c = { 2,	4,	2,	6,	2,	3 };
   vector<int> a = { -3, -6, 2, 6, -3, 8 };
   vector<int> b = { 6, -8, -8, 8, -4, -1 };
   for (size_t i = 0; i < 6; i++)
      sum += c[i] / (1.0 + (point[0] - a[i]) * (point[0] - a[i]) + (point[1] - b[i]) * (point[1] - b[i]));
   return -sum;
}

void Grad(const vector<double>& point, vector<double>& res)
{
   vector<int> c = { 2,	4,	2,	6,	2,	3 };
   vector<int> a = { -3, -6, 2, 6, -3, 8 };
   vector<int> b = { 6, -8, -8, 8, -4, -1 };
   double sum_x = 0., sum_y = 0.;
   for (size_t i = 0; i < 6; i++)
   {
      sum_x -= 2 * c[i] * (point[0] - a[i]) / ((1.0 + (point[0] - a[i]) * (point[0] - a[i]) + (point[1] - b[i]) * (point[1] - b[i])) *
         (1.0 + (point[0] - a[i]) * (point[0] - a[i]) + (point[1] - b[i]) * (point[1] - b[i])));
      sum_y -= 2 * c[i] * (point[1] - b[i]) / ((1.0 + (point[0] - a[i]) * (point[0] - a[i]) + (point[1] - b[i]) * (point[1] - b[i])) *
         (1.0 + (point[0] - a[i]) * (point[0] - a[i]) + (point[1] - b[i]) * (point[1] - b[i])));
   }
   res[0] = sum_x;
   res[1] = sum_y;
}

// Поиск отрезка с минимумом функции
int FindSegmentWithMin(const double& delta, double funct(const vector<double>&),
   const vector<double>& x, const vector<double>& Sk, double& a, double& b)
{
   double x0 = 0;
   double xk, xk1, xk_1, h = 2;
   double f = funct(x + (x0) * Sk);
   int f_calc_count = 1;

   if(f == funct(x + (x0 + delta) * Sk))
   {
      a = x0;
      b = x0 + delta;
      return 2;
   }
   else if(f == funct(x + (x0 - delta) * Sk))
   {
      a = x0 - delta;
      b = x0;
      return 3;
   }
   else
   {
      if(f > funct(x + (x0 + delta) * Sk))
      {
         xk = x0 + delta;
         h = delta;
         f_calc_count++;
      }
      else if(f > funct(x + (x0 - delta) * Sk))
      {
         xk = x0 - delta;
         h = -delta;
         f_calc_count += 2;
      }
      else
      {
         a = x0 - delta;
         b = x0 + delta;
         return f_calc_count + 2;
      }

      xk_1 = x0;
      bool exit = false;
      do
      {
         h *= 2;
         xk1 = xk + h;

         if(funct(x + (xk)*Sk) > funct(x + (xk1)*Sk))
         {
            xk_1 = xk;
            xk = xk1;
         }
         else
            exit = true;

         f_calc_count += 2;
      } while(!exit);

      a = xk_1;
      b = xk;
   }
   return f_calc_count;
}

// Поиск аргумента минимума функции вдоль направления методом золотого сечения
int FindMinArgGolden(double funct(const vector<double>&),
   const vector<double>& x, const vector<double>& Sk, const double& eps, double& result)
{
   double a = 0, b = 0;
   int f_calc_count = FindSegmentWithMin(1e-1, funct, x, Sk, a, b);

   double x1 = a + (3 - SQRT5) / 2 * (b - a);
   double x2 = a + (SQRT5 - 1) / 2 * (b - a);
   double f1 = funct(x + (x1)*Sk);
   double f2 = funct(x + (x2)*Sk);
   double a1, b1;

   int iter_count = 0;
   for(; abs(b - a) > eps; iter_count++)
   {
      a1 = a, b1 = b;
      if(f1 < f2)
      {
         b = x2;
         x2 = x1;
         x1 = a + (3 - SQRT5) / 2 * (b - a);
         f2 = f1;
         f1 = funct(x + (x1)*Sk);
      }
      else
      {
         a = x1;
         x1 = x2;
         x2 = a + (SQRT5 - 1) / 2 * (b - a);
         f1 = f2;
         f2 = funct(x + (x2)*Sk);
      }
   }
   result =  a;
   return iter_count + f_calc_count;
}