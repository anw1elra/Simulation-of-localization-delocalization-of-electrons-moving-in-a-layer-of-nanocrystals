#include <iostream>
#include <random>
#include <fstream>
#include <chrono>
#include <string>

//Reciprocal electron mass in nm^2/fs^2/eV
const double RME = 0.176;
const std::string DIR_NAME = "results";

//Input parameters
double me, tau, u0, a, b, l, dt, ee, eb;
unsigned nn, ne, nc1;

//Auxillary parameters
double tf, cm, ctau, ca, cb, cwx, cwy, dc, w;
unsigned nc, nt;

int execute(const std::string &str) {
    return system(str.c_str());
}

double runtime_in_milliseconds(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

int main() {
    std::fstream fin("input.txt");
    fin >> me >> tau >> u0 >> a >> b >> l;
    fin >> nn >> dt >> ee >> eb >> ne >> nc1;

    auto start = std::chrono::steady_clock::now();

    nc = nc1 * nc1;
    tf = 6 * tau;
    nt = tf / dt;
    cm = RME / me;
    ctau = 1 / tau;
    ca = 0.5 / a / a;
    cb = 0.5 / b / b;
    cwx = 2 * cm * ca * u0;
    cwy = 2 * cm * cb * u0;
    dc = l / nc1;
    w = 3 * a;

    if (!execute("mkdir " + DIR_NAME)) {
        std::cout << "Unable to create folder." << std::endl;
        return 1;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed: " << runtime_in_milliseconds(start, end);
    return 0;
}