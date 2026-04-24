#include "program.h"
#include "circuit.h"
#include "gatedata.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>
#include <string>

namespace sim = sdim;

// helper
bool close_enough(double observed, double expected, double tol) {
    return std::abs(observed - expected) < tol;
}

// ------------------ BIT FLIP TEST ------------------
void test_n1_bit_flip() {
    const int shots = 100000;
    const double p = 0.3;
    const double tol = 0.02;

    sim::Circuit circuit(1, 2);

    circuit.add_gate("N1", 0, std::nullopt, sim::ParameterMap{
        {"channel", std::string("f")},
        {"prob", p}
    });

    circuit.add_gate("M", 0);

    sim::Program program(circuit, 11111);
    auto results = program.simulate_shots(shots, false);

    int count_1 = 0;
    for (const auto& r : results[0][0]) {
        if (r.measurement_value == 1) count_1++;
    }

    double f1 = static_cast<double>(count_1) / shots;

    std::cout << "\n[Bit Flip Test]\n";
    std::cout << "Observed P(1) = " << f1 << ", expected = " << p << "\n";

    assert(close_enough(f1, p, tol));
}

// ------------------ PHASE FLIP TEST ------------------
void test_n1_phase_flip() {
    const int shots = 100000;
    const double p = 0.3;
    const double tol = 0.02;

    sim::Circuit circuit(1, 2);

    // Convert phase error into measurable bit flip
    circuit.add_gate("H", 0);

    circuit.add_gate("N1", 0, std::nullopt, sim::ParameterMap{
        {"channel", std::string("p")},
        {"prob", p}
    });

    circuit.add_gate("H", 0);
    circuit.add_gate("M", 0);

    sim::Program program(circuit, 22222);
    auto results = program.simulate_shots(shots, false);

    int count_1 = 0;
    for (const auto& r : results[0][0]) {
        if (r.measurement_value == 1) count_1++;
    }

    double f1 = static_cast<double>(count_1) / shots;

    std::cout << "\n[Phase Flip Test]\n";
    std::cout << "Observed P(1) = " << f1 << ", expected = " << p << "\n";

    assert(close_enough(f1, p, tol));
}

// ------------------ DEPOLARIZING TEST ------------------
void test_n1_depolarizing() {
    const int shots = 100000;
    const double p = 0.3;
    const double tol = 0.02;

    sim::Circuit circuit(1, 2);

    circuit.add_gate("N1", 0, std::nullopt, sim::ParameterMap{
        {"channel", std::string("d")},
        {"prob", p}
    });

    circuit.add_gate("M", 0);

    sim::Program program(circuit, 33333);
    auto results = program.simulate_shots(shots, false);

    int count_1 = 0;
    for (const auto& r : results[0][0]) {
        if (r.measurement_value == 1) count_1++;
    }

    double f1 = static_cast<double>(count_1) / shots;

    double expected = 2.0 * p / 3.0;

    std::cout << "\n[Depolarizing Test]\n";
    std::cout << "Observed P(1) = " << f1 << ", expected = " << expected << "\n";

    assert(close_enough(f1, expected, tol));
}

// ------------------ MAIN ------------------
int main() {
    test_n1_bit_flip();
    test_n1_phase_flip();
    test_n1_depolarizing();

    std::cout << "\nAll N1 channel tests passed!\n";
    return 0;
}
