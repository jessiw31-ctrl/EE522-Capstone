#include "program.h"
#include "circuit.h"
#include "gatedata.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <optional>

namespace sim = sdim;

bool close_enough(double observed, double expected, double tol) {
    return std::abs(observed - expected) < tol;
}

void test_depolarizing_1qubit() {
    const int shots = 10000;
    const double p = 0.1;
    const double tol = 0.02;

    sim::Circuit circuit(1, 2);

    circuit.add_gate(
        "N1",
        0,
        std::nullopt,
        sim::ParameterMap{
            {"channel", std::string("d")},
            {"prob", p}
        }
    );

    circuit.add_gate("M", 0);

    sim::Program program(circuit, 12345);
    auto results = program.simulate_shots(shots, false);

    int count_0 = 0;
    int count_1 = 0;

    for (const auto& result : results[0][0]) {
        if (result.measurement_value == 0) count_0++;
        else if (result.measurement_value == 1) count_1++;
    }

    double f0 = static_cast<double>(count_0) / shots;
    double f1 = static_cast<double>(count_1) / shots;

    double expected_1 = 2.0 * p / 3.0;
    double expected_0 = 1.0 - expected_1;

    std::cout << "\nDepolarizing noise test\n";
    std::cout << "Observed P(0) = " << f0 << ", expected = " << expected_0 << "\n";
    std::cout << "Observed P(1) = " << f1 << ", expected = " << expected_1 << "\n";

    assert(close_enough(f0, expected_0, tol));
    assert(close_enough(f1, expected_1, tol));
}

void test_pauli_channel_1qubit() {
    const int shots = 100000;
    const double prob_x = 0.20;
    const double prob_z = 0.10;
    const double prob_xz = 0.05;
    const double tol = 0.02;

    sim::Circuit circuit(1, 2);

    circuit.add_gate(
        "PAULI1",
        0,
        std::nullopt,
        sim::ParameterMap{
            {"prob_x", prob_x},
            {"prob_z", prob_z},
            {"prob_xz", prob_xz}
        }
    );

    circuit.add_gate("M", 0);

    sim::Program program(circuit, 54321);
    auto results = program.simulate_shots(shots, false);

    int count_0 = 0;
    int count_1 = 0;

    for (const auto& result : results[0][0]) {
        if (result.measurement_value == 0) count_0++;
        else if (result.measurement_value == 1) count_1++;
    }

    double f0 = static_cast<double>(count_0) / shots;
    double f1 = static_cast<double>(count_1) / shots;

    double expected_1 = prob_x + prob_xz;
    double expected_0 = 1.0 - expected_1;

    std::cout << "\nPauli-channel noise test\n";
    std::cout << "Observed P(0) = " << f0 << ", expected = " << expected_0 << "\n";
    std::cout << "Observed P(1) = " << f1 << ", expected = " << expected_1 << "\n";

    assert(close_enough(f0, expected_0, tol));
    assert(close_enough(f1, expected_1, tol));
}

int main() {
    test_depolarizing_1qubit();
    test_pauli_channel_1qubit();

    std::cout << "\nAll noise tests passed!\n";
    return 0;
}
