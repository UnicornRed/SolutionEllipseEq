#pragma once

#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>

// Класс, общий для всех методов решения дифференциальных уравнений второго порядка с двумя переменными
class SolutionEq
{
private:
    // Коэффициент p
    std::function<double(double, double)> pFunc;
    // Коэффициент q
    std::function<double(double, double)> qFunc;
    // Правая часть уравнения
    std::function<double(double, double)> fFunc;
    // Точное решение (граничные условия)
    std::function<double(double, double)> uFunc;

    // Векторная норма Чебышёва
    double VectorNorm(const std::vector<double>& v) const;
protected:
    double maxX, maxY, c1, c2, d1, d2;
    double hX, hY, epsilon, diffUUPrev, roK;
    size_t N, M, k, m;
    std::vector<double> u;
    std::vector<double> u_k;
    std::vector<double> u0;
    std::vector<double> uExSol;
    std::vector<double> p;
    std::vector<double> q;
    std::vector<double> f;

    // Установка сетки
    virtual void setStep();

    // Максимальное количество итераций
    virtual size_t culc_m() = 0;

    // Численое значение оператора
    double valueLu(size_t i, size_t j, const std::vector<double>& v) const;

    // Вычисление параметров перед каждой итерацией
    virtual void setParamK() {};

    // Следующее приближение
    virtual double nextUk(size_t i, size_t j) = 0;
public:
    SolutionEq(double _maxX, double _maxY, size_t _N, size_t _M,
               std::function<double(double, double)> _pFunc,
               std::function<double(double, double)> _qFunc,
               std::function<double(double, double)> _fFunc,
               std::function<double(double, double)> _uFunc,
               double _epsilon = 0.001);

    inline const size_t getN() const { return N; }
    inline const size_t getM() const { return M; }
    inline const size_t getK() const { return k; }
    inline const size_t get_m() const { return m; }
    inline const std::vector<double> getU() const { return u; }
    inline const std::vector<double> getUExSol() const { return uExSol; }
    inline const std::vector<double> getU0() const { return u0; }
    inline const std::vector<double> getP() const { return p; }
    inline const std::vector<double> getQ() const { return q; }

    virtual const std::string getNameMethod() const = 0;

    // Подготовка состояния класса перед вычислениями
    virtual void preparing();

    // Мера апроксимации
    double ApproxMeasure(const std::vector<double>& v) const;

    // Норма погрешности
    double NormOfError(const std::vector<double>& v1, const std::vector<double>& v2) const;

    void NewStep(size_t _N, size_t _M);

    virtual void printHeader() const;

    virtual void printInfoLine() const = 0;

    void printGrid(const std::vector<double>& v) const;

    bool Next();
};