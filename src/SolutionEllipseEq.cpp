#include "SolutionEllipseEq.h"

SolutionEllipseEq::SolutionEllipseEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     double _epsilon) : SolutionEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc,
                                                                   _fFunc, _uFunc, _epsilon)
{
    
}

void SolutionEllipseEq::setStep()
{
    SolutionEq::setStep();

    double sincosHelp = sin(M_PI * hX / (2 * maxX));
    delta = c1 * 4 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = sin(M_PI * hY / (2 * maxY));
    delta += d1 * 4 / (hY * hY) * (sincosHelp * sincosHelp);

    Delta = culc_Delta();

    xi = culc_xi();
    ro = (1 - xi) / (1 + xi);

    m = culc_m();
}

void SolutionEllipseEq::printHeader() const
{
    SolutionEq::printHeader();

    std::cout << "3. Number of iterations:  m = " << get_m() << "\n" << std::endl;

    std::cout << "4. Spectral radius:  rho(H) = " << getRo() << "\n" << std::endl;

    std::cout << " ------------------------------------------------------------------------------------------------------------------------------\n" <<
                 " |  5.1  |      5.2      |       5.3      |       5.4      |       5.5      |       5.6       |      5.7      |       5.8     |\n" <<
                 " ------------------------------------------------------------------------------------------------------------------------------\n" <<
                 " |   k   |   ||F-AU^k||  |     rel.d.     |  ||U^k - U_*|| |    rel.error   | ||U^k-U^(k-1)|| |   apost.est.  |    sp.rad._k  |\n" <<
                 " ------------------------------------------------------------------------------------------------------------------------------\n";
}

void SolutionEllipseEq::printInfoLine() const
{
    printf(" | %5lu | %13.6f | %14.6f | %14.6f | %14.6f | %15.6f | %13.6f | %13.6f |\n", getK(), ApproxMeasure(getU()),
           ApproxMeasure(getU()) / ApproxMeasure(getU0()), NormOfError(getU(), getUExSol()),
           NormOfError(getU(), getUExSol()) / NormOfError(getU0(), getUExSol()),
           diffUUPrev,
           getRo() * diffUUPrev / (1 - getRo()), k == 1 ? 0 : roK);
}

IOSUMethodSolElEq::IOSUMethodSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     double _epsilon) : SolutionEllipseEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

double IOSUMethodSolElEq::culc_Delta()
{
    double sincosHelp, _Delta;

    sincosHelp = cos(M_PI * hX / (2 * maxX));
    _Delta = c2 * 4.0 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = cos(M_PI * hY / (2 * maxY));
    _Delta += d2 * 4.0 / (hY * hY) * (sincosHelp * sincosHelp);

    return _Delta;
}

double IOSUMethodSolElEq::culc_xi()
{
    return delta / Delta;
}

IterMethodSolElEq::IterMethodSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     double _epsilon) : IOSUMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

size_t IterMethodSolElEq::culc_m()
{
    return ceil(log(1.0 / epsilon) / (2 * xi));
}

SimpIterSolElEq::SimpIterSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                 std::function<double(double, double)> _pFunc,
                                 std::function<double(double, double)> _qFunc,
                                 std::function<double(double, double)> _fFunc,
                                 std::function<double(double, double)> _uFunc,
                                 double _epsilon) : IterMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

const std::string SimpIterSolElEq::getNameMethod() const
{
    return "Simple iteration method";
}

double SimpIterSolElEq::nextUk(size_t i, size_t j)
{
    return ((p[i - 1 + (j - 1) * N] * u[i - 1 + j * (N + 1)]) / (hX * hX) +
            (p[i + (j - 1) * N] * u[i + 1 + j * (N + 1)]) / (hX * hX) +
            (q[i - 1 + (j - 1) * (N - 1)] * u[i + (j - 1) * (N + 1)]) / (hY * hY) +
            (q[i - 1 + j * (N - 1)] * u[i + (j + 1) * (N + 1)]) / (hY * hY) +
             f[i - 1 + (j - 1) * (N - 1)]) /
           ( p[i - 1 + (j - 1) * N] / (hX * hX) + p[i + (j - 1) * N] / (hX * hX) +
             q[i - 1 + (j - 1) * (N - 1)] / (hY * hY) + q[i - 1 + j * (N - 1)] / (hY * hY));
}

OptimParamSolElEq::OptimParamSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     double _epsilon) : IterMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

void OptimParamSolElEq::preparing()
{
    SolutionEllipseEq::preparing();
    culc_tau();
}

void OptimParamSolElEq::culc_tau()
{
    tau = 2.0 / (Delta + delta);
}

const std::string OptimParamSolElEq::getNameMethod() const
{
    return "Iteration method with optimal parameter";
}

double OptimParamSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + tau *
           ((p[i + (j - 1) * N] * (u[i + 1 + j * (N + 1)] - u[i + j * (N + 1)])) / (hX * hX) -
            (p[i - 1 + (j - 1) * N] * (u[i + j * (N + 1)] - u[i - 1 + j * (N + 1)])) / (hX * hX) +
            (q[i - 1 + j * (N - 1)] * (u[i + (j + 1) * (N + 1)] - u[i + j * (N + 1)])) / (hY * hY) -
            (q[i - 1 + (j - 1) * (N - 1)] * (u[i + j * (N + 1)] - u[i + (j - 1) * (N + 1)])) / (hY * hY) +
             f[i - 1 + (j - 1) * (N - 1)]);
}

SeidelSolElEq::SeidelSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                             std::function<double(double, double)> _pFunc,
                             std::function<double(double, double)> _qFunc,
                             std::function<double(double, double)> _fFunc,
                             std::function<double(double, double)> _uFunc,
                             double _epsilon) : IOSUMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

size_t SeidelSolElEq::culc_m()
{
    return ceil(log(1.0 / epsilon) / (4 * xi));
}

const std::string SeidelSolElEq::getNameMethod() const
{
    return "Seidel's method";
}

double SeidelSolElEq::nextUk(size_t i, size_t j)
{
    return ((p[i - 1 + (j - 1) * N] * u_k[i - 1 + j * (N + 1)]) / (hX * hX) +
            (p[i + (j - 1) * N] * u[i + 1 + j * (N + 1)]) / (hX * hX) +
            (q[i - 1 + (j - 1) * (N - 1)] * u_k[i + (j - 1) * (N + 1)]) / (hY * hY) +
            (q[i - 1 + j * (N - 1)] * u[i + (j + 1) * (N + 1)]) / (hY * hY) +
             f[i - 1 + (j - 1) * (N - 1)]) /
           ( p[i - 1 + (j - 1) * N] / (hX * hX) + p[i + (j - 1) * N] / (hX * hX) +
             q[i - 1 + (j - 1) * (N - 1)] / (hY * hY) + q[i - 1 + j * (N - 1)] / (hY * hY));
}

UpRelaxSolElEq::UpRelaxSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                               std::function<double(double, double)> _pFunc,
                               std::function<double(double, double)> _qFunc,
                               std::function<double(double, double)> _fFunc,
                               std::function<double(double, double)> _uFunc,
                               double _epsilon) : IOSUMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{
    omega = 2 / (1 + std::sqrt(1 - getRo() * getRo()));
}

void UpRelaxSolElEq::preparing()
{
    SolutionEllipseEq::preparing();
    omega = 2 / (1 + std::sqrt(1 - getRo() * getRo()));
}

size_t UpRelaxSolElEq::culc_m()
{
    return ceil(log(1.0 / epsilon) / std::sqrt(xi));
}

const std::string UpRelaxSolElEq::getNameMethod() const
{
    return "Upper relaxation method";
}

double UpRelaxSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + omega *
           ((p[i + (j - 1) * N] * (u[i + 1 + j * (N + 1)] - u[i + j * (N + 1)])) / (hX * hX) -
            (p[i - 1 + (j - 1) * N] * (u[i + j * (N + 1)] - u_k[i - 1 + j * (N + 1)])) / (hX * hX) +
            (q[i - 1 + j * (N - 1)] * (u[i + (j + 1) * (N + 1)] - u[i + j * (N + 1)])) / (hY * hY) -
            (q[i - 1 + (j - 1) * (N - 1)] * (u[i + j * (N + 1)] - u_k[i + (j - 1) * (N + 1)])) / (hY * hY) +
             f[i - 1 + (j - 1) * (N - 1)]) /
           ( p[i - 1 + (j - 1) * N] / (hX * hX) + p[i + (j - 1) * N] / (hX * hX) +
             q[i - 1 + (j - 1) * (N - 1)] / (hY * hY) + q[i - 1 + j * (N - 1)] / (hY * hY));
}

ChebParam::ChebParam(size_t _p) : n(1)
{
    for (size_t i{0}; i < _p; ++i)
        n *= 2;

    theta.resize(n);

    std::vector<size_t> th = {1};

    for (size_t m{2}; m <= n; m *= 2)
    {
        std::copy(th.begin(), th.end(), theta.begin());

        th.resize(m);

        for (size_t i{1}; i < m; ++i)
            if (i % 2)
                th[i] = 2 * m - th[i - 1];
            else
                th[i] = theta[i / 2];
    }

    std::copy(th.begin(), th.end(), theta.begin());
}

ChebParamSolElEq::ChebParamSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                   std::function<double(double, double)> _pFunc,
                                   std::function<double(double, double)> _qFunc,
                                   std::function<double(double, double)> _fFunc,
                                   std::function<double(double, double)> _uFunc,
                                   size_t _p, double _epsilon) : IOSUMethodSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon), ChebParam(_p)
{
    
}

size_t ChebParamSolElEq::culc_m()
{
    return ceil(log(2.0 / epsilon) / (2 * std::sqrt(xi)));
}

const std::string ChebParamSolElEq::getNameMethod() const
{
    return "Method with an optimal Chebyshev set of parameters";
}

void ChebParamSolElEq::setParamK()
{
    tau_k = 2.0 / (Delta + delta + (Delta - delta) * std::cos(double(theta[(k - 1) % n]) / (2.0 * n) * M_PI));
}

double ChebParamSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + tau_k *
           ((p[i + (j - 1) * N] * (u[i + 1 + j * (N + 1)] - u[i + j * (N + 1)])) / (hX * hX) -
            (p[i - 1 + (j - 1) * N] * (u[i + j * (N + 1)] - u[i - 1 + j * (N + 1)])) / (hX * hX) +
            (q[i - 1 + j * (N - 1)] * (u[i + (j + 1) * (N + 1)] - u[i + j * (N + 1)])) / (hY * hY) -
            (q[i - 1 + (j - 1) * (N - 1)] * (u[i + j * (N + 1)] - u[i + (j - 1) * (N + 1)])) / (hY * hY) +
             f[i - 1 + (j - 1) * (N - 1)]);
}

GeneralTriangSolElEq::GeneralTriangSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                           std::function<double(double, double)> _pFunc,
                                           std::function<double(double, double)> _qFunc,
                                           std::function<double(double, double)> _fFunc,
                                           std::function<double(double, double)> _uFunc,
                                           double _epsilon) : SolutionEllipseEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

double GeneralTriangSolElEq::culc_Delta()
{
    return c2 * 4.0 / (hX * hX) + d2 * 4.0 / (hY * hY);
}

double GeneralTriangSolElEq::culc_xi()
{
    eta = delta / Delta;

    gamma_1 = delta / (2 + 2 * std::sqrt(eta));
    gamma_2 = delta / (4 * std::sqrt(eta));

    return gamma_1 / gamma_2;
}

void GeneralTriangSolElEq::preparing()
{
    SolutionEllipseEq::preparing();
    omega = 2.0 / std::sqrt(delta * Delta);
    kappa_1 = omega / (hX * hX);
    kappa_2 = omega / (hY * hY);

    wSet.resize((N + 1) * (M + 1));
    wk.resize((N + 1) * (M + 1));
}

void GeneralTriangSolElEq::setParamK()
{
    for (size_t j{1}; j < M; ++j)
        for (size_t i{1}; i < N; ++i)
            wSet[i + j * (N + 1)] = (kappa_1 * p[i - 1 + (j - 1) * N] * wSet[i - 1 + j * (N + 1)] +
                                     kappa_2 * q[i - 1 + (j - 1) * (N - 1)] * wSet[i + (j - 1) * (N + 1)] +
                                     valueLu(i, j, u) + f[i - 1 + (j - 1) * (N - 1)]) / 
                                    (1 + kappa_1 * p[i - 1 + (j - 1) * N] + kappa_2 * q[i - 1 + (j - 1) * (N - 1)]);

    for (size_t j{M - 1}; j >= 1; --j)
        for (size_t i{N - 1}; i >= 1; --i)
            wk[i + j * (N + 1)] = (kappa_1 * p[i + (j - 1) * N] * wk[i + 1 + j * (N + 1)] +
                                   kappa_2 * q[i - 1 + j * (N - 1)] * wk[i + (j + 1) * (N + 1)] +
                                   wSet[i + j * (N + 1)]) / 
                                  (1 + kappa_1 * p[i + (j - 1) * N] + kappa_2 * q[i - 1 + j * (N - 1)]);
}

AlterTriangSolElEq::AlterTriangSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                       std::function<double(double, double)> _pFunc,
                                       std::function<double(double, double)> _qFunc,
                                       std::function<double(double, double)> _fFunc,
                                       std::function<double(double, double)> _uFunc,
                                       double _epsilon) : GeneralTriangSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon)
{

}

void AlterTriangSolElEq::preparing()
{
    GeneralTriangSolElEq::preparing();
    tau = 2.0 / (gamma_1 + gamma_2);
}

size_t AlterTriangSolElEq::culc_m()
{
    return ceil(log(1.0 / epsilon) / log(1.0 / ro));
}

const std::string AlterTriangSolElEq::getNameMethod() const
{
    return "Alternately triangular iterative method";
}

double AlterTriangSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + tau * wk[i + j * (N + 1)];
}

ChebAltTriSolElEq::ChebAltTriSolElEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     size_t _p, double _epsilon) : GeneralTriangSolElEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc, _fFunc, _uFunc, _epsilon), ChebParam(_p)
{

}

size_t ChebAltTriSolElEq::culc_m()
{
    return ceil(log(2.0 / epsilon) / (2 * std::sqrt(2 * std::sqrt(eta))));
}

const std::string ChebAltTriSolElEq::getNameMethod() const
{
    return "Alternately triangular iterative method with Chebyshev set of parameters";
}

void ChebAltTriSolElEq::preparing()
{
    GeneralTriangSolElEq::preparing();
}

void ChebAltTriSolElEq::setParamK()
{
    GeneralTriangSolElEq::setParamK();
#ifdef DEBUG
    std::cout << "Theta " << k << ": " << theta[(k - 1) % n] << "\n";
#endif
    tau = 2.0 / (gamma_2 + gamma_1 + (gamma_2 - gamma_1) * std::cos(double(theta[(k - 1) % n]) / (2.0 * n) * M_PI));
    // tau = 2.0 / (gamma_2 + gamma_1 + (gamma_2 - gamma_1) * std::cos(double(2 * k - 1) / (2.0 * n) * M_PI));
}

double ChebAltTriSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + tau * wk[i + j * (N + 1)];
}