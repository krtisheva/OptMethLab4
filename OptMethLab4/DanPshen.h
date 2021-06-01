#pragma once
#include <vector>
#include "Lab2Data.h"
#include "Vector.h"

class DanPshen
{
public:

   const int size;            // ����������� �������
   vector<double> Sk;         // ����� �������� �� ������� ����
   vector<double> xk;         // ����������� �� ������� ����
   vector<double> xk1;        // ����� �����������
   vector<double> grad_xk1;   // ��������������� ������ ��� ���������� ����������
   vector<double> t;          // ��������������� ������ ��� ���������� ���������
   int f_calc_count = 0;      // ����� ���������� �������

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
      // 1. ������ ��������� ������� f � ����� x0
      Grad(x0, Sk);
      f_calc_count += 4;
      Sk *= -1;
      xk = x0;

      bool exit_flag = true;
      int iter_count = 0;
      do
      {
         // 2. ����������� ������� f �� ����������� Sk
         //Function to_minimize = Function(xk, Sk);
         double lambda;
         f_calc_count += min_max(funct, xk, Sk, 1e-5, lambda) + 1;

         // ��������� ������ �����������
         xk1 = xk + lambda * Sk;

         // 3. ���������� grad(f(xk1)) � �������� ������������ omega
         Grad(xk1, grad_xk1);
         f_calc_count += 4;

         //double omega = (grad_xk1 * (grad_xk1 + Sk)) / (-1 * (Sk * Sk));
         double omega = (grad_xk1 * grad_xk1) / (Sk * Sk);

         // 4. ����������� ����� ����������� Sk1
         Sk = -1 * grad_xk1 + omega * Sk;

         // ������ ��������� ������� �� ������� ��������
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