#pragma once
#include <vector>
#include "Lab2Data.h"
#include "Vector.h"

class DanPshen
{
public:

   const int size;            // Размерность вектора
   vector<double> Sk;         // Минус гардиент на текущем шаге
   vector<double> xk;         // Приближение на текущем шаге
   vector<double> xk1;        // Новое приближение
   vector<double> grad_xk1;   // Вспомогательный вектор для нахождения экстремума
   vector<double> t;          // Вспомогательный вектор для нахождения градиента
   int f_calc_count = 0;      // Число вычислений функции

   DanPshen() : size(0)
   {
   }

   DanPshen(const int& t_size) : size(t_size)
   {
      Sk.resize(size);
      xk.resize(size);
      xk1.resize(size);
      grad_xk1.resize(size);
      t.resize(size);
   }

   int FindExtremum(double funct(const vector<double>&),
                    int min_max(double f(const vector<double>&), const vector<double>&, const vector<double>&, const double&, double&),
                    const vector<double>& x0, const double& f_eps, const double& xs_eps)
   {
      f_calc_count = 0;
      // 1. Расчет градиента функции f в точке x0
      Grad(x0, Sk);
      f_calc_count += 4;
      Sk *= -1;
      xk = x0;

      bool exit_flag = true;
      int iter_count = 0;
      do
      {
         // 2. Минимизация функции f по направлению Sk
         //Function to_minimize = Function(xk, Sk);
         double lambda;
         f_calc_count += min_max(funct, xk, Sk, 1e-5, lambda) + 1;

         // Получение нового приближения
         xk1 = xk + lambda * Sk;

         // 3. Вычисление grad(f(xk1)) и весового коэффициента omega
         Grad(xk1, grad_xk1);
         f_calc_count += 4;

         //double omega = (grad_xk1 * (grad_xk1 + Sk)) / (-1 * (Sk * Sk));
         double omega = (grad_xk1 * grad_xk1) / (Sk * Sk);

         // 4. Определение новго направления Sk1
         Sk = -1 * grad_xk1 + omega * Sk;

         // Расчет изменения решения на текущей итерации
         exit_flag = true;

         for(int i = 0; i < size; i++)
            if(abs(xk[i] - xk1[i]) > xs_eps)
               exit_flag = false;

         xk = xk1;
         iter_count++;

      } while(Norm(Sk) > f_eps && iter_count < 100 && exit_flag == false);
      return iter_count;
   }
};