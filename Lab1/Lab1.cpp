#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>

using namespace std;

struct TableEntry {
    double X, T, U;
};

vector<TableEntry> readTable(const string& filename) {
    vector<TableEntry> table;
    ifstream fin(filename);
    if (!fin.is_open()) return table;
    string line;
    if (!getline(fin, line)) return table;
    while (getline(fin, line)) {
        for (auto& c : line) {
            if (c == ',') c = '.';
        }
        if (line.empty() || (!isdigit(line.front()) && line.front() != '-' && line.front() != '+'))
            continue;
        istringstream iss(line);
        TableEntry entry;
        if (iss >> entry.X >> entry.T >> entry.U)
            table.push_back(entry);
    }
    return table;
}

static double interpolate(double x, double x1, double y1, double x2, double y2) {
    if (fabs(x2 - x1) < 1e-10) return y1;
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

pair<double, double> get_TU(double x) {
    vector<TableEntry> table;
    string fname;
    if (fabs(x) <= 1.0)
        fname = "dat_X_1_1.dat";
    else if (x < -1.0)
        fname = "dat_X00_1.dat";
    else if (x > 1.0)
        fname = "dat_X1_00.dat";
    else
        return { NAN, NAN };

    table = readTable(fname);
    if (table.empty()) {
        //cout << "[ERROR] Table is empty!" << endl;
        return { NAN, NAN };
    }

    for (size_t i = 0; i + 1 < table.size(); ++i) {
        if ((x >= table[i].X && x <= table[i + 1].X) ||
            (x <= table[i].X && x >= table[i + 1].X)) {

            double T = interpolate(x, table[i].X, table[i].T, table[i + 1].X, table[i + 1].T);
            double U = interpolate(x, table[i].X, table[i].U, table[i + 1].X, table[i + 1].U);
            //cout << "[INFO] Interpolation: x = " << x << ", T = " << T << ", U = " << U << endl;
            return { T, U };
        }
    }

    //cout << "[WARN] x = " << x << " is out of table range" << endl;
    return { NAN, NAN };
}



double Srz(double x, double y, double z) {
    if (x > y) {
        auto Tx = get_TU(x).first;
        auto Uz = get_TU(z).second;
        auto Ty = get_TU(y).first;
        if (isnan(Tx) || isnan(Uz) || isnan(Ty)) return NAN;
        return Tx + Uz - Ty;
    }
    else {
        auto Ty = get_TU(y).first;
        auto Uy = get_TU(y).second; 
        auto Uz = get_TU(z).second;
        if (isnan(Ty) || isnan(Uz) || isnan(Uy)) return NAN;
        return Ty + Uy - Uz;     
    }
}
double Gold(double x, double y) {
    if (y != 0) {
        if (x > y) return x / y;
        if (x < y) return y / x;
    }
    return NAN;
}

double Glr(double x, double y) {
    if (fabs(x) < 1) return x;
    if (fabs(y) >= 1 && fabs(y) < 1) return y;
    double val = y / sqrt(x * x + y * y - 4);
    if (fabs(y) >= 1 && (x * x + y * y - 4) > 0.1) return val;
    if ((x * x + y * y - 4) < 0.1) return NAN;
    return NAN;
}

double Grs(double x, double y) {
    double s1 = 0.1389 * Srz(x + y, Gold(x, y), Glr(x, x * y));
    double s2 = 1.8389 * Srz(x - y, Gold(y, x / 5), Glr(5 * x, Gold(5 * y, y)));
    return s1 + s2 * 0.83;
}

double fun2(double x, double y, double z) {
    auto Grsl = [](double a, double b) {
        double s1 = 0.14 * Srz(a + b, Gold(a, b), Glr(a, a * b));
        double s2 = 1.83 * Srz(a - b, Gold(b, a / 5), Glr(4 * a, Gold(4 * b, b)));
        return s1 + s2 * 0.83;
        };
    double result = x * Grsl(x, y) + y * Grsl(y, z) + z * Grsl(z, x);
    return result;
}

double fun3(double x, double y, double z) {
    return 1.3498 * z + 2.2362 * y - 2.348 * x * y;
}

double fun1(double x, double y, double z) {
    auto Grs1 = [&](double a, double b) {
        double gold = Gold(a, b);
        if (isnan(gold)) return fun2(x, y, z);
        double glr = Glr(a, b);
        if (isnan(glr)) return fun2(x, y, z);
        return Grs(a, b);
        };
    double val = x * x * Grs1(y, z) + y * y * Grs1(x, z) + 0.33 * x * y * Grs1(x, z);
    return val;
}

double fun(double x, double y, double z) {
    double result = fun1(x, y, z);
    if (!isnan(result)) {
        //cout << "[DEBUG] Using fun1()" << endl;
        return result;
    }
    result = fun2(x, y, z);
    if (!isnan(result)) {
        //cout << "[DEBUG] Using fun2()" << endl;
        return result;
    }
    //cout << "[DEBUG] Using fun3()" << endl;
    return fun3(x, y, z);
}


int main() {
    double x, y, z;
    cout << "Enter x, y, z: ";
    cin >> x >> y >> z;

    double res = fun(x, y, z);
    cout << "fun(x, y, z) = " << res << endl;
    return 0;
}
