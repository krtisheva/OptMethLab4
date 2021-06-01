#pragma once
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

// ��������� ������� �� �����
vector<double> operator * (const double& val, vector<double> vec)
{
   const size_t size = vec.size();
   
   for (size_t i = 0; i < size; ++i)
      vec[i] *= val;
   return vec;
}

// ������� ������� �� �����
vector<double> operator / (const double& val, vector<double> vec)
{
   const size_t size = vec.size();

   for(size_t i = 0; i < size; ++i)
      vec[i] /= val;
   return vec;
}

vector<double>& operator *= (vector<double>& vec, const double& val)
{
   const size_t size = vec.size();

   for(size_t i = 0; i < size; ++i)
      vec[i] *= val;

   return vec;
}

// �������� ��������
vector<double> operator + (vector<double> vec1, const vector<double>& vec2)
{
   const size_t size = vec1.size();

   for (size_t i = 0; i < size; ++i)
      vec1[i] += vec2[i];

   return vec1;
}

// ��������� ��������
vector<double> operator - (vector<double>vec1, const vector<double>& vec2)
{
   const size_t size = vec1.size();

   for (size_t i = 0; i < size; ++i)
      vec1[i] -= vec2[i];
   return vec1;
}

// ��������� ������������ ��������
double operator * (const vector<double>& vec1, const vector<double>& vec2)
{
   const size_t size = vec1.size();

   double res = 0;

   for(size_t i = 0; i < size; ++i)
      res += vec1[i] * vec2[i];

   return res;
}

// ����� �������
double Norm(const vector<double>& vec)
{
   const size_t size = vec.size();

   double res = 0;
   for(int i = 0; i < size; i++)
      res += vec[i] * vec[i];

   return sqrt(res);
}