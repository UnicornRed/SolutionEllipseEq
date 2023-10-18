#include "SolutionEq.h"

SolutionEq::SolutionEq(double _maxX, double _maxY, size_t _N, size_t _M,
                       std::function<double(double, double)> _pFunc,
                       std::function<double(double, double)> _qFunc,
                       std::function<double(double, double)> _fFunc,
                       std::function<double(double, double)> _uFunc,
                       double _epsilon) : pFunc(_pFunc), qFunc(_qFunc), fFunc(_fFunc), uFunc(_uFunc),
                                          maxX(_maxX), maxY(_maxY), epsilon(_epsilon),
                                          N(_N), M(_M), k(0)
{
    
}

double SolutionEq::valueLu(size_t i, size_t j, const std::vector<double>& v) const
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

double SolutionEq::VectorNorm(const std::vector<double>& v) const
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

double SolutionEq::ApproxMeasure(const std::vector<double>& v) const
{
    std::vector<double> ApMeas((N - 1) * (M - 1));

    for (size_t j{}; j < M - 1; ++j)
        for (size_t i{}; i < N - 1; ++i)
            ApMeas[i + j * (N - 1)] = valueLu(i + 1, j + 1, v) + f[i + j * (N - 1)];

    return VectorNorm(ApMeas);
}

double SolutionEq::NormOfError(const std::vector<double>& v1, const std::vector<double>& v2) const
{
    std::vector<double> NormEr(v1.begin(), v1.end());

    for (size_t i{}; i < v1.size(); ++i)
        NormEr[i] -= v2[i];

    return VectorNorm(NormEr);
}

void SolutionEq::NewStep(size_t _N, size_t _M)
{
    N = _N;
    M = _M;
    k = 0;
    u.clear();
    u_k.clear();
    u0.clear();
    uExSol.clear();
    p.clear();
    q.clear();
    f.clear();

    setStep();
}

void SolutionEq::setStep()
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
}

void SolutionEq::preparing()
{
    setStep();
}

void SolutionEq::printHeader() const
{
    std::cout << getNameMethod() << ". Variant 14.\n" << std::endl;

    std::cout << "1. Measure of approximation: ||F-AU_*|| = " << ApproxMeasure(getUExSol()) << "\n" << std::endl;

    std::cout << "2. Discrepancy norm for U^0:  ||F-AU^0|| = " << ApproxMeasure(getU0()) << "\n" << std::endl;
}

void SolutionEq::printGrid(const std::vector<double>& v) const
{
    if (N % 5 != 0 && M % 5 != 0)
    {
        std::cout << "Grids error! N and M must be % 5 == 0.\n";

        return;
    }

    size_t multiX = N / 5, multiY = M / 5;

    printf("-----------");

    for (size_t i{}; i <= N / multiX; ++i)
        printf("--------------");

    printf("\n|    x\\y  |");

    for (size_t i{}; i <= N / multiX; ++i)
        printf(" %11.2f |", hX * (i * multiX));

    printf("\n");

    for (size_t j{}; j <= M / multiY; ++j)
    {
        printf("-----------");

        for (size_t i{}; i <= N / multiX; ++i)
            printf("--------------");

        printf("\n| %7.2f |", hY * (j * multiY));

        for (size_t i{}; i <= N / multiX; ++i)
            printf(" %11.5f |", v[i * multiX + (j * multiY) * (N + 1)]);

        printf("\n");
    }

    printf("-----------");

    for (size_t i{}; i <= N / multiX; ++i)
        printf("--------------");
}

bool SolutionEq::Next()
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
            std::cout << "(" << i << ", " << j << ") - " << u_k[i + j * (N + 1)] << ";";
        
        std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    return false;
}