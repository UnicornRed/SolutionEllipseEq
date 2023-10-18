#pragma once

#include "SolutionEq.h"

// Общий класс для всех численных методов решения эллиптических уравнений
class SolutionEllipseEq : public SolutionEq
{
private:
    void setStep() override;
protected:
    double delta, Delta, xi, ro;

    virtual double culc_Delta() = 0;

    virtual double culc_xi() = 0;
public:
    SolutionEllipseEq(double _maxX, double _maxY, size_t _N, size_t _M,
                      std::function<double(double, double)> _pFunc,
                      std::function<double(double, double)> _qFunc,
                      std::function<double(double, double)> _fFunc,
                      std::function<double(double, double)> _uFunc,
                      double _epsilon = 0.001);

    inline const double getRo() const { return ro; }

    void printHeader() const override;

    void printInfoLine() const override;
};

// Класс общий для простых итерационных методов, Зейделя, верхней релаксации и метода с чебышевским набором параметров
class IOSUMethodSolElEq : public SolutionEllipseEq
{
private:
    double culc_Delta() override;

    double culc_xi() override;
public:
    IOSUMethodSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                      std::function<double(double, double)> _pFunc,
                      std::function<double(double, double)> _qFunc,
                      std::function<double(double, double)> _fFunc,
                      std::function<double(double, double)> _uFunc,
                      double _epsilon = 0.001);
};

// Класс общий для итерационных методов: простого и с оптимальным параметром
class IterMethodSolElEq : public IOSUMethodSolElEq
{
protected:
    size_t culc_m() override;
public:
    IterMethodSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                     std::function<double(double, double)> _pFunc,
                     std::function<double(double, double)> _qFunc,
                     std::function<double(double, double)> _fFunc,
                     std::function<double(double, double)> _uFunc,
                     double _epsilon = 0.001);
};

// Класс простого итерационного метода
class SimpIterSolElEq final : public IterMethodSolElEq
{
protected:
    double nextUk(size_t i, size_t j) override;
public:
    SimpIterSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                    std::function<double(double, double)> _pFunc,
                    std::function<double(double, double)> _qFunc,
                    std::function<double(double, double)> _fFunc,
                    std::function<double(double, double)> _uFunc,
                    double _epsilon = 0.001);

    const std::string getNameMethod() const override;
};

// Класс итерационного метода с оптимальным параметром
class OptimParamSolElEq final : public IterMethodSolElEq
{
private:
    double tau;
protected:
    double nextUk(size_t i, size_t j) override;
public:
    OptimParamSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                      std::function<double(double, double)> _pFunc,
                      std::function<double(double, double)> _qFunc,
                      std::function<double(double, double)> _fFunc,
                      std::function<double(double, double)> _uFunc,
                      double _epsilon = 0.001);

    inline double getTau() const { return tau; }

    const std::string getNameMethod() const override;

    void preparing() override;

    void culc_tau();
};

// Класс метода Зейделя
class SeidelSolElEq final : public IOSUMethodSolElEq
{
protected:
    double nextUk(size_t i, size_t j) override;

    size_t culc_m() override;
public:
    SeidelSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                  std::function<double(double, double)> _pFunc,
                  std::function<double(double, double)> _qFunc,
                  std::function<double(double, double)> _fFunc,
                  std::function<double(double, double)> _uFunc,
                  double _epsilon = 0.001);

    const std::string getNameMethod() const override;
};

// Класс метода верхней релаксации
class UpRelaxSolElEq final : public IOSUMethodSolElEq
{
private:
    double omega;
protected:
    double nextUk(size_t i, size_t j) override;

    size_t culc_m() override;
public:
    UpRelaxSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                   std::function<double(double, double)> _pFunc,
                   std::function<double(double, double)> _qFunc,
                   std::function<double(double, double)> _fFunc,
                   std::function<double(double, double)> _uFunc,
                   double _epsilon = 0.001);

    inline double getOmega() const { return omega; }

    const std::string getNameMethod() const override;

    void preparing() override;
};

// Класс чебышевских параметров
class ChebParam
{
protected:
    size_t n;
    std::vector<size_t> theta;
public:
    ChebParam(size_t _p);
};

// Класс метода итерации с чебышевским набором параметров 
class ChebParamSolElEq final : public IOSUMethodSolElEq, public ChebParam
{
private:
    double tau_k;
protected:
    void setParamK() override;

    double nextUk(size_t i, size_t j) override;

    size_t culc_m() override;
public:
    ChebParamSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                     std::function<double(double, double)> _pFunc,
                     std::function<double(double, double)> _qFunc,
                     std::function<double(double, double)> _fFunc,
                     std::function<double(double, double)> _uFunc,
                     size_t _p, double _epsilon = 0.001);

    const std::string getNameMethod() const override;
};

// Класс общий для попеременно-треугольных методов
class GeneralTriangSolElEq : public SolutionEllipseEq
{
protected:
    double omega, tau, gamma_1, gamma_2, kappa_1, kappa_2, eta;

    double culc_Delta() override;

    double culc_xi() override;

    void setParamK() override;

    std::vector<double> wSet, wk;
public:
    GeneralTriangSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                         std::function<double(double, double)> _pFunc,
                         std::function<double(double, double)> _qFunc,
                         std::function<double(double, double)> _fFunc,
                         std::function<double(double, double)> _uFunc,
                         double _epsilon = 0.001);

    void preparing() override;
};

// Класс попеременно-треугольного метода
class AlterTriangSolElEq final : public GeneralTriangSolElEq
{
private:

protected:
    double nextUk(size_t i, size_t j) override;

    size_t culc_m() override;
public:
    AlterTriangSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                       std::function<double(double, double)> _pFunc,
                       std::function<double(double, double)> _qFunc,
                       std::function<double(double, double)> _fFunc,
                       std::function<double(double, double)> _uFunc,
                       double _epsilon = 0.001);

    const std::string getNameMethod() const override;

    void preparing() override;
};

// Класс попеременно-треугольного метода с чебышевским набором параметров
class ChebAltTriSolElEq final : public GeneralTriangSolElEq, public ChebParam
{
private:

protected:
    double nextUk(size_t i, size_t j) override;

    size_t culc_m() override;

    void setParamK() override;
public:
    ChebAltTriSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                      std::function<double(double, double)> _pFunc,
                      std::function<double(double, double)> _qFunc,
                      std::function<double(double, double)> _fFunc,
                      std::function<double(double, double)> _uFunc,
                      size_t _p, double _epsilon = 0.001);

    const std::string getNameMethod() const override;

    void preparing() override;
};