#pragma once

#include "SolutionEq.h"

// Класс переменных направлений
class SolutionPoisEq : public SolutionEq
{
private:
    
protected:
    double delta, Delta, tau;

    std::vector<double> uMiddle;
    std::vector<double> uNext;
    std::vector<double> s1;
    std::vector<double> t1;
    std::vector<double> s2;
    std::vector<double> t2;
public:
    SolutionPoisEq(double _maxX, double _maxY, size_t _N, size_t _M,
                   std::function<double(double, double)> _pFunc,
                   std::function<double(double, double)> _qFunc,
                   std::function<double(double, double)> _fFunc,
                   std::function<double(double, double)> _uFunc,
                   double _epsilon = 0.001);

    void setStep() override;

    size_t culc_m() override;

    void setParamK() override;

    double nextUk(size_t i, size_t j) override;

    const std::string getNameMethod() const override;

    void preparing() override;

    void printHeader() const override;

    void printInfoLine() const override;
};
