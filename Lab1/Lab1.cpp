#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

struct TableRow {
    double X, T, U;
};

std::vector<TableRow> table;

std::vector<TableRow> readTable(const std::string& filename);
double interpolate(const std::vector<TableRow>& table, double x, char col);
double T(double x);
double U(double x);

double Gold(double x, double y);
double Glr(double x, double y);
double Srz(double x, double y, double z);
double Grs(double x, double y);
double fun1(double x, double y, double z);

double Gold2(double x, double y);
double Glr2(double x, double y);
double Srz2(double x, double y, double z);
double Grs2(double x, double y);
double fun2(double x, double y, double z);

double fun3(double x, double y, double z) {
    return 1.3498 * x + 2.2362 * y - 2.348 * x * y;
}

std::vector<TableRow> readTable(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) throw std::runtime_error("Cannot open file: " + filename);

    std::vector<TableRow> table;
    std::string line;
    while (std::getline(fin, line)) {
        std::istringstream iss(line);
        TableRow row;
        if (!(iss >> row.X >> row.T >> row.U)) continue;
        table.push_back(row);
    }
    if (table.empty()) throw std::runtime_error("File is empty: " + filename);
    return table;
}

double interpolate(const std::vector<TableRow>& table, double x, char col) {
    for (size_t i = 1; i < table.size(); ++i) {
        double x0 = table[i - 1].X, x1 = table[i].X;
        double y0 = (col == 'T') ? table[i - 1].T : table[i - 1].U;
        double y1 = (col == 'T') ? table[i].T : table[i].U;
        if ((x0 <= x && x <= x1) || (x1 <= x && x <= x0)) {
            if (x1 == x0) return y0;
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        }
    }
    if (x < table.front().X)
        return (col == 'T') ? table.front().T : table.front().U;
    if (x > table.back().X)
        return (col == 'T') ? table.back().T : table.back().U;
    throw std::runtime_error("Interpolation error: x out of bounds");
}

double T(double x) {
    return interpolate(table, x, 'T');
}
double U(double x) {
    return interpolate(table, x, 'U');
}

void selectTable(double x) {
    std::string fname;
    if (std::abs(x) <= 1.0) fname = "dat_X_1_1.dat";
    else if (x < -1.0) fname = "dat_X00_1.dat";
    else if (x > 1.0) fname = "dat_X1_00.dat";
    else throw std::runtime_error("Cannot select table for x");
    table = readTable(fname);
}

double Srz(double x, double y, double z) {
    double denom = (T(y) + U(y) - U(z));
    if (fabs(denom) < 1e-12) throw std::runtime_error("Srz: denominator is zero");
    return (T(x) + U(z) - T(y)) / denom;
}
double Gold(double x, double y) {
    if (x > y && y == 0) throw std::runtime_error("Gold: y=0, x>y, needs recalculation");
    if (x < y && x == 0) throw std::runtime_error("Gold: x=0, x<y, needs recalculation");
    if (x > y) return x * y;
    if (x < y) return x / y;
    return 0.0;
}
double Glr(double x, double y) {
    if (fabs(x) < 1.0) return x;
    if (fabs(y) < 1.0) return y;
    double rad = x * x + y * y - 4.0;
    if (fabs(x) >= 1.0 && fabs(y) >= 1.0 && rad > 0.1) return sqrt(rad);
    throw std::runtime_error("Glr: arguments not in allowed range");
}
double Grs(double x, double y) {
    return 0.1389 * Srz(x + y, Gold(x, y), Glr(x, x * y)) +
        1.8389 * Srz(x - y, Gold(y, x / 5.0), Glr(5.0 * x, x * y)) +
        0.33 * Srz(x - 0.9, Glr(y, x / 5.0), Gold(5.0 * y, y));
}
double fun1(double x, double y, double z) {
    return x * x * Grs(y, z) + y * y * Grs(x, z) + 0.33 * x * y * Grs(x, z);
}

double Srz2(double x, double y, double z) {
    double denom = (T(y) + U(y) - U(z));
    if (fabs(denom) < 1e-12) throw std::runtime_error("Srz2: denominator is zero");
    return (T(x) + U(z) - T(y)) / denom;
}
double Gold2(double x, double y) {
    if (y > 0.1 && x > y) return x / y;
    if (y > 0.1 && x <= y) return 0.15;
    if (y == 0) return 0.1;
    return 0.0;
}
double Glr2(double x, double y) {
    if (fabs(x) < 1.0) return x;
    if (fabs(y) < 1.0) return y;
    return 0.0;
}
double Grs2(double x, double y) {
    return 0.14 * Srz2(x + y, Gold2(x, y), Glr2(x, x * y)) +
        1.83 * Srz2(x - y, Gold2(y, x / 5.0), Glr2(4.0 * x, x * y)) +
        0.33 * Srz2(x, Glr2(y, x / 4.0), Gold2(4.0 * y, y));
}
double fun2(double x, double y, double z) {
    return x * Grs2(x, y) + y * Grs2(y, z) + z * Grs2(z, x);
}

double fun(double x, double y, double z) {
    try {
        selectTable(x);
    }
    catch (const std::exception& e) {
        std::cerr << "File not opened: " << e.what() << "\n";
        return fun3(x, y, z);
    }

    try {
        return fun1(x, y, z);
    }
    catch (const std::exception& e) {
        std::cerr << "fun1 failed: " << e.what() << "\n";
        try {
            return fun2(x, y, z);
        }
        catch (const std::exception& e2) {
            std::cerr << "fun2 failed: " << e2.what() << "\n";
            return fun3(x, y, z);
        }
    }
}

int main() {
    double x, y, z;
    std::cout << "Enter x, y, z: ";
    std::cin >> x >> y >> z;

    try {
        double res = fun(x, y, z);
        std::cout << "fun(x, y, z) = " << std::setprecision(10) << res << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}