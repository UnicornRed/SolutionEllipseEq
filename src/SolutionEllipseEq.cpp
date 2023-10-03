#include <algorithm>
#include <cmath>
#include "SolutionEllipseEq.h"

// #define DEBUG

SolutionEllipseEq::SolutionEllipseEq(double _maxX, double _maxY, size_t _N, size_t _M,
                                     std::function<double(double, double)> _pFunc,
                                     std::function<double(double, double)> _qFunc,
                                     std::function<double(double, double)> _fFunc,
                                     std::function<double(double, double)> _uFunc,
                                     double _epsilon) : pFunc(_pFunc), qFunc(_qFunc), fFunc(_fFunc), uFunc(_uFunc),
                                                                                     maxX(_maxX), maxY(_maxY), epsilon(_epsilon),
                                                                                     N(_N), M(_M), k(0)
{
    
}

double SolutionEllipseEq::valueLu(size_t i, size_t j, const std::vector<double>& v)
{
    double value = p[i + (j - 1) * N] * (v[i + 1 + j * (N + 1)] - v[i + j * (N + 1)]) / (hX * hX) - 
                   p[i - 1 + (j - 1) * N] * (v[i + j * (N + 1)] - v[i - 1 + j * (N + 1)]) / (hX * hX) +
                   q[i - 1 + j * (N - 1)] * (v[i + (j + 1) * (N + 1)] - v[i + j * (N + 1)]) / (hY * hY) - 
                   q[i - 1 + (j - 1) * (N - 1)] * (v[i + j * (N + 1)] - v[i + (j - 1) * (N + 1)]) / (hY * hY);

    #ifdef DEBUG
        std::cout << "valueLu: (" << i << ", " << j << ") - " << value << "\nv: " << std::endl;

        for (const auto& i : v)
            std::cout << i << " ";

        std::cout << "\np: " << std::endl;

        for (const auto& i : p)
            std::cout << i << " ";

        std::cout << "\nq: " << std::endl;

        for (const auto& i : q)
            std::cout << i << " ";

        std::cout << std::endl;
    #endif

    return value;
}

double SolutionEllipseEq::VectorNorm(const std::vector<double>& v)
{
    double norm = 0;

    std::for_each(v.begin(), v.end(), [&norm](const double& a)
    {
        double abs_a = std::abs(a);

        if (abs_a > norm)
            norm = abs_a;
    });

    return norm;
}

double SolutionEllipseEq::ApproxMeasure(const std::vector<double>& v)
{
    std::vector<double> ApMeas((N - 1) * (M - 1));

    for (size_t j{}; j < M - 1; ++j)
        for (size_t i{}; i < N - 1; ++i)
            ApMeas[i + j * (N - 1)] = valueLu(i + 1, j + 1, v) + f[i + j * (N - 1)];

    return VectorNorm(ApMeas);
}

double SolutionEllipseEq::NormOfError(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> NormEr(v1.begin(), v1.end());

    for (size_t i{}; i < v1.size(); ++i)
        NormEr[i] -= v2[i];

    return VectorNorm(NormEr);
}

void SolutionEllipseEq::setStep()
{
    hX = maxX / N;
    hY = maxY / M;

    u.resize((N + 1) * (M + 1));
    u_k.resize((N + 1) * (M + 1));
    u0.resize((N + 1) * (M + 1));
    uExSol.resize((N + 1) * (M + 1));
    p.resize(N * (M - 1));
    q.resize((N - 1) * M);
    f.resize((N - 1) * (M - 1));

    for (size_t i{}; i < N + 1; ++i)
    {
        u[i] = uFunc(hX * i, 0);
        u[i + M * (N + 1)] = uFunc(hX * i, maxY);
    }

    for (size_t j{1}; j < M; ++j)
    {
        u[j * (N + 1)] = uFunc(0, hY * j);
        u[N + j * (N + 1)] = uFunc(maxX, hY * j);
    }

    std::copy(u.begin(), u.end(), u0.begin());

    for (size_t j{}; j < M + 1; ++j)
        for (size_t i{}; i < N + 1; ++i)
            uExSol[i + j * (N + 1)] = uFunc(hX * i, hY * j);

    for (size_t j{1}; j < M; ++j)
        for (size_t i{}; i < N; ++i)
            p[i + (j - 1) * N] = pFunc(hX * i + hX / 2, hY * j);
    
    for (size_t j{}; j < M; ++j)
        for (size_t i{1}; i < N; ++i)
            q[(i - 1) + j * (N - 1)] = qFunc(hX * i, hY * j + hY / 2);

    for (size_t j{}; j < M - 1; ++j)
        for (size_t i{}; i < N - 1; ++i)
            f[i + j * (N - 1)] = fFunc(hX * (i + 1), hY * (j + 1));

    #ifdef DEBUG
        std::cout << "u:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << u[i + j * (N + 1)] << ";";
            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "uExSol:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << uExSol[i + j * (N + 1)] << ";";
                            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "p:\n";
        for (size_t j{1}; j < M; ++j)
        {
            for (size_t i{}; i < N; ++i)
                std::cout << "(" << 0.5 + i << ", " << j << ") - " << p[i + (j - 1) * N] << ";";
                            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "q:\n";
        for (size_t j{}; j < M; ++j)
        {
            for (size_t i{1}; i < N; ++i)
                std::cout << "(" << i << ", " << 0.5 + j << ") - " << q[(i - 1) + j * (N - 1)] << ";";
                            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "f:\n";
        for (size_t j{}; j < M - 1; ++j)
        {
            for (size_t i{}; i < N - 1; ++i)
                std::cout << "(" << i + 1 << ", " << j + 1 << ") - " << f[i + j * (N - 1)] << ";";
                            
            std::cout << std::endl;
        }
        std::cout << std::endl;
    #endif

    c1 = pFunc(0, 0);
    c2 = pFunc(0, 0);
    d1 = qFunc(0, 0);
    d2 = qFunc(0, 0);

    for (double x = 0; x <= maxX + epsilon; x += hX)
    {
        for (double y = 0; y <= maxY + epsilon; y += hY)
        {
            c1 = std::min(c1, pFunc(x, y));
            c2 = std::max(c2, pFunc(x, y));

            d1 = std::min(d1, qFunc(x, y));
            d2 = std::min(d2, qFunc(x, y));
        }
    }

    double sincosHelp = sin(M_PI * hX / (2 * maxX));
    delta = c1 * 4 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = sin(M_PI * hY / (2 * maxY));
    delta += d1 * 4 / (hY * hY) * (sincosHelp * sincosHelp);

    Delta = culc_Delta();

    xi = culc_xi();
    ro = (1 - xi) / (1 + xi);

    m = culc_m();
}

void SolutionEllipseEq::preparing()
{
    setStep();
}

void SolutionEllipseEq::printInfoLine(const std::string& format)
{
    printf(format.c_str(), getK(), ApproxMeasure(getU()),
           ApproxMeasure(getU()) / ApproxMeasure(getU0()), NormOfError(getU(), getUExSol()),
           NormOfError(getU(), getUExSol()) / NormOfError(getU0(), getUExSol()),
           diffUUPrev,
           getRo() * diffUUPrev / (1 - getRo()), k == 1 ? 0 : roK);
}

void SolutionEllipseEq::printGrid(const std::vector<double>& v)
{
    if (N % 5 != 0 && M % 5 != 0)
    {
        std::cout << "Grids error! N and M must be % 5 == 0.\n";

        return;
    }

    size_t multi = N / 5;

    printf("-----------");

    for (size_t i{}; i <= N / multi; ++i)
        printf("--------------");

    printf("\n|    x\\y  |");

    for (size_t i{}; i <= N / multi; ++i)
        printf(" %11.2f |", hX * (i * multi));

    printf("\n");

    for (size_t j{}; j <= M / multi; ++j)
    {
        printf("-----------");

        for (size_t i{}; i <= N / multi; ++i)
            printf("--------------");

        printf("\n| %7.2f |", hY * (j * multi));

        for (size_t i{}; i <= N / multi; ++i)
            printf(" %11.5f |", v[i * multi + (j * multi) * (N + 1)]);

        printf("\n");
    }

    printf("-----------");

    for (size_t i{}; i <= N / multi; ++i)
        printf("--------------");
}

bool SolutionEllipseEq::Next()
{
    if (k >= m)
        return true;
    
    ++k;

    std::copy(u.begin(), u.end(), u_k.begin());

    setParamK();

    for (size_t j{1}; j < M; ++j)
        for (size_t i{1}; i < N; ++i)
            u_k[i + j * (N + 1)] = nextUk(i, j);

    #ifdef DEBUG
        std::cout << "u_k:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << u_k[i + j * (N + 1)] << ";";
            
            std::cout << std::endl;
        }
        std::cout << std::endl;
    #endif

    roK = NormOfError(u, u_k) / diffUUPrev;

    diffUUPrev = NormOfError(u, u_k);

    std::copy(u_k.begin(), u_k.end(), u.begin());

    if (NormOfError(u, uExSol) / NormOfError(u0, uExSol) < epsilon)
        return true;

    #ifdef DEBUG
        std::cout << "u:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << u[i + j * (N + 1)] << ";";
            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "uPrev:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << uPrev[i + j * (N + 1)] << ";";
            
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "uPrevPrev:\n";
        for (size_t j{}; j < M + 1; ++j)
        {
            for (size_t i{}; i < N + 1; ++i)
                std::cout << "(" << i << ", " << j << ") - " << uPrevPrev[i + j * (N + 1)] << ";";
            
            std::cout << std::endl;
        }
        std::cout << std::endl;
    #endif

    return false;
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
    _Delta = c2 * 4 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = cos(M_PI * hY / (2 * maxY));
    _Delta += d2 * 4 / (hY * hY) * (sincosHelp * sincosHelp);

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
    return c2 * 4 / (hX * hX) + d2 * 4 / (hY * hY);
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
}

void GeneralTriangSolElEq::setParamK()
{
    for (size_t j{1}; j < M; ++j)
        for (size_t i{1}; i < N; ++i)
            wSet[i + j * (N + 1)] = (kappa_1 * p[i - 1 + (j - 1) * N] * wSet[i - 1 + j * (N + 1)] +
                                     kappa_2 * q[i - 1 + (j - 1) * (N - 1)] * wSet[i + (j - 1) * (N + 1)] +
                                     valueLu(i, j, u) + f[i - 1 + (j - 1) * (N - 1)]) / 
                                    (1 + kappa_1 * p[i - 1 + (j - 1) * N] + kappa_2 * q[i - 1 + (j - 1) * (N - 1)]);
    
    for (size_t j{M}; j >= 1; --j)
        for (size_t i{N}; i >= 1; --i)
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

    wSet.resize((N + 1) * (M + 1));
    wk.resize((N + 1) * (M + 1));
}

size_t AlterTriangSolElEq::culc_m()
{
    return ceil(log(2.0 / epsilon) / (2 * std::sqrt(xi)));
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

void ChebAltTriSolElEq::preparing()
{
    GeneralTriangSolElEq::preparing();

    wSet.resize((N + 1) * (M + 1));
    wk.resize((N + 1) * (M + 1));
}

size_t ChebAltTriSolElEq::culc_m()
{
    return ceil(log(2.0 / epsilon) / (2 * std::sqrt(2 * std::sqrt(eta))));
}

const std::string ChebAltTriSolElEq::getNameMethod() const
{
    return "Alternately triangular iterative method with Chebyshev set of parameters";
}

void ChebAltTriSolElEq::setParamK()
{
    GeneralTriangSolElEq::setParamK();

    tau = 2.0 / (gamma_2 + gamma_1 + (gamma_2 - gamma_1) * std::cos(double(theta[(k - 1) % n]) / (2.0 * n) * M_PI));
}

double ChebAltTriSolElEq::nextUk(size_t i, size_t j)
{
    return u[i + j * (N + 1)] + tau * wk[i + j * (N + 1)];
}