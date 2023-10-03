#include <iostream>
#include "SolutionEllipseEq.h"

double pCoeff(double x, double y)
{
    return 3 * x + 2;
}

double qCoeff(double x, double y)
{
    return 1;
}

double fParam(double x, double y)
{
    return -x * x * (2 + 6 * y) - (12 * x + 4) * y * y * (1 + y);
}

double uExactSolution(double x, double y)
{
    return x * x * y * y * (1 + y);
}

void WorkWithMethod(SolutionEllipseEq& MyEq)
{
    MyEq.preparing();

    std::cout << MyEq.getNameMethod() << ". Variant 14.\n" << std::endl;

    std::cout << "1. Measure of approximation: ||F-AU_*|| = " << MyEq.ApproxMeasure(MyEq.getUExSol()) << "\n" << std::endl;

    std::cout << "2. Discrepancy norm for U^0:  ||F-AU^0|| = " << MyEq.ApproxMeasure(MyEq.getU0()) << "\n" << std::endl;

    std::cout << "3. Number of iterations:  m = " << MyEq.get_m() << "\n" << std::endl;

    std::cout << "4. Spectral radius:  rho(H) = " << MyEq.getRo() << "\n" << std::endl;

    std::cout << " ------------------------------------------------------------------------------------------------------------------------------\n" <<
                 " |  5.1  |      5.2      |       5.3      |       5.4      |       5.5      |       5.6       |      5.7      |       5.8     |\n" <<
                 " ------------------------------------------------------------------------------------------------------------------------------\n" <<
                 " |   k   |   ||F-AU^k||  |     rel.d.     |  ||U^k - U_*|| |    rel.error   | ||U^k-U^(k-1)|| |   apost.est.  |    sp.rad._k  |\n" <<
                 " ------------------------------------------------------------------------------------------------------------------------------\n";

    while (!MyEq.Next())
    {
        if (MyEq.getK() < 5 || MyEq.getK() % MyEq.getN() == 0)
            MyEq.printInfoLine(" | %5lu | %13.6f | %14.6f | %14.6f | %14.6f | %15.6f | %13.6f | %13.6f |\n");
    }

    MyEq.printInfoLine(" | %5lu | %13.6f | %14.6f | %14.6f | %14.6f | %15.6f | %13.6f | %13.6f |\n");

    std::cout << "\n\nApproximate solution:" << std::endl;

    MyEq.printGrid(MyEq.getU());

    std::cout << "\n\nExact solution:" << std::endl;

    MyEq.printGrid(MyEq.getUExSol());

    std::cout << std::endl << std::endl;
}

int main()
{
    size_t N, M;
    double maxX, maxY, epsilon;

    std::cin >> maxX >> maxY >> epsilon;
    std::cin >> N >> M;

    // Простая итерация

    SimpIterSolElEq MyEqSI(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, epsilon);

    WorkWithMethod(MyEqSI);

    // Оптимальный параметр

    OptimParamSolElEq MyEqOP(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, epsilon);

    WorkWithMethod(MyEqOP);

    std::cout << "Tau: " << MyEqOP.getTau() << std::endl;

    // Метод Зейделя

    SeidelSolElEq MyEqSM(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, epsilon);

    WorkWithMethod(MyEqSM);

    // Метод верхней релаксации

    UpRelaxSolElEq MyEqUR(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, epsilon);

    WorkWithMethod(MyEqUR);

    std::cout << "Omega: " << MyEqUR.getOmega() << std::endl;

    // Метод с чебышевским набором параметров

    ChebParamSolElEq MyEqCP(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, 5, epsilon);

    WorkWithMethod(MyEqCP);

    // Попеременно-треугольный итерационный метод

    AlterTriangSolElEq MyEqAT(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, epsilon);

    WorkWithMethod(MyEqAT);

    // Попеременно-треугольный итерационный метод с чебышевским набором параметров

    ChebAltTriSolElEq MyEqCA(maxX, maxY, N, M, pCoeff, qCoeff, fParam, uExactSolution, 3, epsilon);

    WorkWithMethod(MyEqCA);

    return 0;
}