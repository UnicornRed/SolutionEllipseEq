#include "SolutionPoisEq.h"

SolutionPoisEq::SolutionPoisEq(double _maxX, double _maxY, size_t _N, size_t _M,
                               std::function<double(double, double)> _pFunc,
                               std::function<double(double, double)> _qFunc,
                               std::function<double(double, double)> _fFunc,
                               std::function<double(double, double)> _uFunc,
                               double _epsilon) : SolutionEq(_maxX, _maxY, _N, _M, _pFunc, _qFunc,
                                                                     _fFunc, _uFunc, _epsilon)
{

}

void SolutionPoisEq::setStep()
{
    SolutionEq::setStep();

    double sincosHelp = sin(M_PI * hX / (2 * maxX));
    delta = c1 * 4 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = sin(M_PI * hY / (2 * maxY));
    delta = std::min(delta, d1 * 4 / (hY * hY) * (sincosHelp * sincosHelp));

    sincosHelp = cos(M_PI * hX / (2 * maxX));
    Delta = c2 * 4 / (hX * hX) * (sincosHelp * sincosHelp);
    sincosHelp = cos(M_PI * hY / (2 * maxY));
    Delta = std::max(Delta, d2 * 4 / (hY * hY) * (sincosHelp * sincosHelp));

    m = culc_m();

    for (size_t j{1}; j < M; ++j)
        for (size_t i{1}; i < N; ++i)
            u[i + j * (N + 1)] = 1;
}

size_t SolutionPoisEq::culc_m()
{
    return N / (2 * M_PI) * std::log(1 / epsilon);
}

void SolutionPoisEq::setParamK()
{
    double A, B, C, G, gamma;

    for (size_t j{1}; j < M; ++j)
    {
        A = 0.0;
        B = -1.0;
        C = 0.0;
        G = uExSol[j * (N + 1)];

        s1[0] = C / B;
        t1[0] = - G / B;

        for (size_t i{1}; i < N; ++i)
        {
            A = tau * p[i - 1 + (j - 1) * N] / (2.0 * hX * hX);
            C = tau * p[i + (j - 1) * N] / (2.0 * hX * hX);
            B = A + C + 1.0;
            G = -u[i + j * (N + 1)] - tau / 2.0 *
                (q[i - 1 + j * (N - 1)] * (u[i + (j + 1) * (N + 1)] - u[i + j * (N + 1)]) / (hY * hY) - 
                 q[i - 1 + (j - 1) * (N - 1)] * (u[i + j * (N + 1)] - u[i + (j - 1) * (N + 1)]) / (hY * hY) +
                 f[i - 1 + (j - 1) * (N - 1)]);

            gamma = B - A * s1[i - 1];
            s1[i] = C / gamma;
            t1[i] = (A * t1[i - 1] - G) / gamma;
        }

        A = 0.0;
        B = -1.0;
        C = 0.0;
        G = uExSol[N + j * (N + 1)];

        gamma = B - A * s1[N - 1];
        s1[N] = C / gamma;
        t1[N] = (A * t1[N - 1] - G) / gamma;

        uMiddle[N + j * (N + 1)] = t1[N];

        for (int i{static_cast<int>(N) - 1}; i >= 0; --i)
            uMiddle[i + j * (N + 1)] = s1[i] * uMiddle[i + 1 + j * (N + 1)] + t1[i];
    }

    for (size_t i{1}; i < N; ++i)
    {
        A = 0.0;
        B = -1.0;
        C = 0.0;
        G = uExSol[i];

        s2[0] = C / B;
        t2[0] = - G / B;

        for (size_t j{1}; j < M; ++j)
        {
            A = tau * q[i - 1 + (j - 1) * (N - 1)] / (2.0 * hY * hY);
            C = tau * q[i - 1 + j * (N - 1)] / (2.0 * hY * hY);
            B = A + C + 1.0;
            G = -uMiddle[i + j * (N + 1)] - tau / 2.0 *
                (p[i + (j - 1) * N] * (uMiddle[i + 1 + j * (N + 1)] - uMiddle[i + j * (N + 1)]) / (hX * hX) - 
                 p[i - 1 + (j - 1) * N] * (uMiddle[i + j * (N + 1)] - uMiddle[i - 1 + j * (N + 1)]) / (hX * hX) +
                 f[i - 1 + (j - 1) * (N - 1)]);

            gamma = B - A * s2[j - 1];
            s2[j] = C / gamma;
            t2[j] = (A * t2[j - 1] - G) / gamma;
        }

        A = 0.0;
        B = -1.0;
        C = 0.0;
        G = uExSol[i + M * (N + 1)];

        gamma = B - A * s2[M - 1];
        s2[M] = C / gamma;
        t2[M] = (A * t2[M - 1] - G) / gamma;

        uNext[i + M * (N + 1)] = t2[M];

        for (int j{static_cast<int>(M) - 1}; j >= 0; --j)
            uNext[i + j * (N + 1)] = s2[j] * uNext[i + (j + 1) * (N + 1)] + t2[j];
    }
}

double SolutionPoisEq::nextUk(size_t i, size_t j)
{
    return uNext[i + j * (N + 1)];
}

const std::string SolutionPoisEq::getNameMethod() const
{
    return "The method of variable directions";
}

void SolutionPoisEq::preparing()
{
    SolutionEq::preparing();

    tau = 2.0 / std::sqrt(delta * Delta);

    uMiddle.resize((N + 1) * (M + 1));
    uNext.resize((N + 1) * (M + 1));
    s1.resize((N + 1));
    t1.resize((N + 1));
    s2.resize((M + 1));
    t2.resize((M + 1));

    for (size_t i{}; i < N + 1; ++i)
    {
        uMiddle[i] = uExSol[i];
        uMiddle[i + M * (N + 1)] = uExSol[i + M * (N + 1)];
    }

    for (size_t j{}; j < M + 1; ++j)
    {
        uNext[j * (N + 1)] = uExSol[j * (N + 1)];
        uNext[N + j * (N + 1)] = uExSol[N + j * (N + 1)];
    }
}

void SolutionPoisEq::printHeader() const
{
    SolutionEq::printHeader();

    std::cout << " ----------------------------------------------------------------------------------------------\n" <<
                 " |  5.1  |      5.2      |       5.3      |       5.4      |       5.5      |       5.6       |\n" <<
                 " ----------------------------------------------------------------------------------------------\n" <<
                 " |   k   |   ||F-AU^k||  |     rel.d.     |  ||U^k - U_*|| |    rel.error   | ||U^k-U^(k-1)|| |\n" <<
                 " ----------------------------------------------------------------------------------------------\n";
}

void SolutionPoisEq::printInfoLine() const
{
    printf(" | %5lu | %13.6f | %14.6f | %14.6f | %14.6f | %15.6f |\n", getK(), ApproxMeasure(getU()),
           ApproxMeasure(getU()) / ApproxMeasure(getU0()), NormOfError(getU(), getUExSol()),
           NormOfError(getU(), getUExSol()) / NormOfError(getU0(), getUExSol()),
           diffUUPrev);
}