#include "program.h"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using exstabsim::BuildIrResult;
using exstabsim::Circuit;
using exstabsim::FrameResults3D;
using exstabsim::MeasurementArray2D;
using exstabsim::MeasurementRecord;
using exstabsim::Program;

constexpr int kDimension = 3;

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

int mod_qutrit(int value) {
    int r = value % kDimension;
    return r < 0 ? r + kDimension : r;
}

struct Options {
    int shots = 1000;
    int distance = 3;
    int rounds = 5;
    int num_p_points = 11;
    double p_min = 0.001;
    double p_max = 0.1;
    int logical_value = 0;
    std::uint64_t seed = 0;
    std::string syndrome_output = "syndrome.csv";
    std::string final_data_output = "final_data.csv";
    std::string summary_output = "summary.csv";
    bool export_all_syndromes = false;
};

Options parse_options(int argc, char **argv) {
    Options opt;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        auto val = [&](const std::string &name) {
            if (i + 1 >= argc) throw std::runtime_error("Missing value after " + name);
            return argv[++i];
        };

        if (arg == "--shots") opt.shots = std::stoi(val(arg));
        else if (arg == "--distance") opt.distance = std::stoi(val(arg));
        else if (arg == "--rounds") opt.rounds = std::stoi(val(arg));
        else if (arg == "--num-p-points") opt.num_p_points = std::stoi(val(arg));
        else if (arg == "--p-min") opt.p_min = std::stod(val(arg));
        else if (arg == "--p-max") opt.p_max = std::stod(val(arg));
        else if (arg == "--logical") opt.logical_value = std::stoi(val(arg));
        else if (arg == "--seed") opt.seed = std::stoull(val(arg));
        else if (arg == "--syndrome-output") opt.syndrome_output = val(arg);
        else if (arg == "--final-data-output") opt.final_data_output = val(arg);
        else if (arg == "--summary-output") opt.summary_output = val(arg);
        else if (arg == "--export-all-syndromes") opt.export_all_syndromes = true;
        else throw std::runtime_error("Unknown argument: " + arg);
    }

    require(opt.shots > 0, "shots must be positive");
    require(opt.distance >= 2, "distance must be at least 2");
    require(opt.rounds > 0, "rounds must be positive");
    require(opt.num_p_points > 0, "num-p-points must be positive");
    require(opt.p_min > 0, "p-min must be positive");
    require(opt.p_max >= opt.p_min, "p-max must be >= p-min");
    require(opt.logical_value >= 0 && opt.logical_value < kDimension,
            "logical must be 0, 1, or 2");

    return opt;
}

std::vector<double> make_p_grid(const Options &opt) {
    std::vector<double> ps;

    if (opt.num_p_points == 1) {
        ps.push_back(opt.p_min);
        return ps;
    }

    double log_min = std::log(opt.p_min);
    double log_max = std::log(opt.p_max);

    for (int i = 0; i < opt.num_p_points; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(opt.num_p_points - 1);
        ps.push_back(std::exp(log_min + t * (log_max - log_min)));
    }

    return ps;
}

int ancilla_index(int distance, int edge) {
    return distance + edge;
}

std::size_t ancilla_syndrome_measurement_round(int circuit_round) {
    return static_cast<std::size_t>(2 * circuit_round + 1);
}

Circuit make_circuit(int distance, int rounds, int logical_value, double p_x) {
    Circuit circuit(distance + distance - 1, kDimension);

    for (int q = 0; q < distance; ++q) {
        for (int k = 0; k < logical_value; ++k) {
            circuit.add_gate("X", q);
        }
    }

    const exstabsim::ParameterMap noise_params{
        {"prob_x", p_x},
        {"prob_z", 0.0},
        {"prob_xz", 0.0},
    };

    for (int round = 0; round < rounds; ++round) {
        for (int q = 0; q < distance; ++q) {
            circuit.add_gate("PAULI1", q, std::nullopt, noise_params);
        }

        for (int edge = 0; edge < distance - 1; ++edge) {
            int ancilla = ancilla_index(distance, edge);

            circuit.add_gate("RESET", ancilla);
            circuit.add_gate("CNOT", edge, ancilla);
            circuit.add_gate("CNOT_INV", edge + 1, ancilla);
            circuit.add_gate("M", ancilla);
        }
    }

    std::vector<int> data_qudits;
    for (int q = 0; q < distance; ++q) {
        data_qudits.push_back(q);
    }

    circuit.add_gate("M", data_qudits);

    return circuit;
}

const MeasurementRecord &measurement_at(
        const FrameResults3D &measurements,
        int qudit,
        std::size_t round,
        std::size_t shot) {

    return measurements.at(static_cast<std::size_t>(qudit), round, shot);
}

int measurement_value_at(
        const FrameResults3D &measurements,
        int qudit,
        std::size_t round,
        std::size_t shot) {

    return mod_qutrit(
        static_cast<int>(
            measurement_at(measurements, qudit, round, shot).measurement_value));
}

int final_data_check_value(
        const FrameResults3D &measurements,
        int edge,
        std::size_t shot) {

    int left = measurement_value_at(measurements, edge, 0, shot);
    int right = measurement_value_at(measurements, edge + 1, 0, shot);

    return mod_qutrit(left - right);
}

void export_syndrome_defects(
        std::ofstream &syndrome_out,
        const FrameResults3D &measurements,
        double p,
        int distance,
        int rounds,
        int shot,
        bool export_all_syndromes) {

    for (int edge = 0; edge < distance - 1; ++edge) {
        int previous_syndrome = 0;

        for (int round = 0; round < rounds; ++round) {
            int ancilla = ancilla_index(distance, edge);
            std::size_t meas_round = ancilla_syndrome_measurement_round(round);

            int syndrome = measurement_value_at(
                measurements,
                ancilla,
                meas_round,
                static_cast<std::size_t>(shot));

            int defect = mod_qutrit(syndrome - previous_syndrome);

            if (export_all_syndromes || defect != 0) {
                syndrome_out << p << ','
                             << distance << ','
                             << rounds << ','
                             << shot << ','
                             << round << ','
                             << edge << ','
                             << 0 << ','
                             << syndrome << ','
                             << previous_syndrome << ','
                             << defect << '\n';
            }

            previous_syndrome = syndrome;
        }

        int final_syndrome = final_data_check_value(
            measurements,
            edge,
            static_cast<std::size_t>(shot));

        int final_defect = mod_qutrit(final_syndrome - previous_syndrome);

        if (export_all_syndromes || final_defect != 0) {
            syndrome_out << p << ','
                         << distance << ','
                         << rounds << ','
                         << shot << ','
                         << rounds << ','
                         << edge << ','
                         << 1 << ','
                         << final_syndrome << ','
                         << previous_syndrome << ','
                         << final_defect << '\n';
        }
    }
}

void export_final_data(
        std::ofstream &final_data_out,
        const FrameResults3D &measurements,
        double p,
        int distance,
        int rounds,
        int shot) {

    for (int q = 0; q < distance; ++q) {
        int value = measurement_value_at(
            measurements,
            q,
            0,
            static_cast<std::size_t>(shot));

        final_data_out << p << ','
                       << distance << ','
                       << rounds << ','
                       << shot << ','
                       << q << ','
                       << value << '\n';
    }
}

int majority_vote(
        const FrameResults3D &measurements,
        int distance,
        int shot) {

    std::vector<int> counts(kDimension, 0);

    for (int q = 0; q < distance; ++q) {
        int value = measurement_value_at(
            measurements,
            q,
            0,
            static_cast<std::size_t>(shot));

        counts[value]++;
    }

    int best_symbol = -1;
    int best_count = -1;
    bool tied = false;

    for (int symbol = 0; symbol < kDimension; ++symbol) {
        if (counts[symbol] > best_count) {
            best_symbol = symbol;
            best_count = counts[symbol];
            tied = false;
        } else if (counts[symbol] == best_count) {
            tied = true;
        }
    }

    return tied ? -1 : best_symbol;
}
}

int main(int argc, char **argv) {
    try {
        Options opt = parse_options(argc, argv);
        std::vector<double> ps = make_p_grid(opt);

        std::ofstream syndrome_out(opt.syndrome_output);
        std::ofstream final_data_out(opt.final_data_output);
        std::ofstream summary_out(opt.summary_output);

        syndrome_out << "p,distance,rounds,shot,round,edge,"
                     << "is_final_boundary,syndrome,previous_syndrome,defect\n";

        final_data_out << "p,distance,rounds,shot,q,value\n";

        summary_out << "p,distance,rounds,shots,raw_majority_failures,"
                    << "raw_majority_logical_error_rate\n";

        syndrome_out << std::setprecision(17);
        final_data_out << std::setprecision(17);
        summary_out << std::setprecision(17);

        for (double p : ps) {
            Circuit circuit = make_circuit(
                opt.distance,
                opt.rounds,
                opt.logical_value,
                p);

            Program program(
                circuit,
                opt.seed + static_cast<std::uint64_t>(std::llround(1e9 * p)));

            (void) program.simulate();

            MeasurementArray2D reference_results =
                Program::results_to_array(program.measurement_results());

            BuildIrResult ir_result = program.build_ir(opt.shots);

            FrameResults3D measurements =
                Program::simulate_frame(
                    ir_result.ir,
                    reference_results,
                    circuit.num_qudits(),
                    circuit.dimension(),
                    opt.shots,
                    ir_result.noise,
                    0);

            int raw_failures = 0;

            for (int shot = 0; shot < opt.shots; ++shot) {
                export_syndrome_defects(
                    syndrome_out,
                    measurements,
                    p,
                    opt.distance,
                    opt.rounds,
                    shot,
                    opt.export_all_syndromes);

                export_final_data(
                    final_data_out,
                    measurements,
                    p,
                    opt.distance,
                    opt.rounds,
                    shot);

                int logical_out = majority_vote(
                    measurements,
                    opt.distance,
                    shot);

                if (logical_out < 0 || logical_out != opt.logical_value) {
                    raw_failures++;
                }
            }

            double raw_logical_error_rate =
                static_cast<double>(raw_failures) /
                static_cast<double>(opt.shots);

            summary_out << p << ','
                        << opt.distance << ','
                        << opt.rounds << ','
                        << opt.shots << ','
                        << raw_failures << ','
                        << raw_logical_error_rate << '\n';

            std::cout << "[pass] p=" << p
                      << " d=" << opt.distance
                      << " rounds=" << opt.rounds
                      << " raw logical error rate="
                      << raw_logical_error_rate << "\n";
        }

        std::cout << "[pass] wrote " << opt.syndrome_output << "\n";
        std::cout << "[pass] wrote " << opt.final_data_output << "\n";
        std::cout << "[pass] wrote " << opt.summary_output << "\n";

        return EXIT_SUCCESS;

    } catch (const std::exception &ex) {
        std::cerr << "[fail] " << ex.what() << "\n";
        return EXIT_FAILURE;
    }
}
