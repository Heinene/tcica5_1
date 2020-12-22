#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>

const double a = 0.0; const double b = 3.0; const double c = -1.0; const double d = 3.0; const size_t N = 10; const double A = 3.0; const double step = (double)((b - a) / (N - 1));

struct Point {
    double x;
    double y;
};

double f(const double& x) {
    return c * x + d;
}

std::vector<Point> random(const size_t N, const double shu) {
    std::vector<Point> points(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> error(-0.5, 0.5);
    for (size_t i = 0; i < N; ++i) {
        points[i].x = a + i * step;
        points[i].y = f(a + i * step) + shu * error(gen);
    }
    return points;
}

std::vector<Point> edge(const double niz, const double vver,
    const size_t num, const double shu) {
    std::vector<Point> points(num);
    const double step = (vver - niz) / static_cast<double>(num - 1);
    for (size_t i = 0; i < num; ++i) {
        points[i].x = niz + i * step;
        points[i].y = f(points[i].x) - shu / 2 + rand() * 1. / RAND_MAX * (shu);
    }
    return points;
}

double error(const std::vector<Point>& pra, const double c, const double d) {
    double sum = 0.;
    for (auto point : pra) {
        sum += pow(point.y - (c * f(point.x)), 2);
    }
    return sum;
}

double golden_ratio(std::vector<Point>& p, double Cm, double C) {
    double l = std::abs(C - Cm);
    std::swap(Cm, C);
    Cm = std::fabs(Cm);
    C = std::fabs(C);
    const double e = 0.1;
    const double t = (std::sqrt(5) + 1) / 2;
    double c1 = Cm + (1 - 1 / t) * C;
    double c2 = Cm + C / t;
    double f1 = f(-c1);
    double f2 = f(-c2);
    while (l > e) {
        if (f1 < f2) {
            C = c2;
            c2 = Cm + C - c1;
            f2 = f(-c2);
        }
        else {
            Cm = c1;
            c1 = Cm + C - c2;
            f1 = f(-c1);
        }
        if (c1 > c2) {
            std::swap(c1, c2);
            std::swap(f1, f2);
        }
        l = std::abs(C - Cm);
    }
    return -((C + Cm) / 2);
}

int F(int f)
{
    if (f == 1) return 1;
    else if (f == 2) return 1;
    else if (f > 2)
        return (F(f - 1) + F(f - 2));
}

double Fibon(std::vector<Point>& p, double Dm, double D) {
    double a = Dm, b = D, x1, x2, y1, y2;
    int G = 10;
    x1 = a + (double)F(G - 2) / F(G) * (b - a);
    x2 = a + (double)F(G - 1) / F(G) * (b - a);
    y1 = error(p, x1, 0);
    y2 = error(p, x2, 0);
    for (int i = G; i >= 1; --i) {
        if (y1 > y2) {
            a = x1;
            x1 = x2;
            x2 = b - (x2 - a);
            y1 = y2;
            y2 = error(p, x2, 0);
        }
        else {
            b = x2;
            x2 = x1;
            x1 = a + (b - x2);
            y2 = y1;
            y1 = error(p, x1, 0);
        }
    }
    return (x1 + x2) / 2;
}

void print(const double noise)
{
    std::vector<Point> p = random(N, noise);
    double Cm, C, Dm, D;
    edge(Cm, C, Dm, D);
    std::cout << "Cm = " << Cm << "\nC = " << C << "\nDm = " << Dm << "\nD = " << D << std::endl;
    double w1 = golden_ratio(p, Cm, C);
    double w0 = Fibon(p, Dm, D);
    std::cout << "w1 = " << w1 << std::endl;
    std::cout << "w0 = " << w0 << std::endl;
}


int main() {
    std::cout << "bezshu" << std::endl;
    print(0.0);
    std::cout << "s shu " << A << std::endl;
    print(A);
    return 0;
}