
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
    using sdim::BuildIrResult;
    using sdim::Circuit;
    using sdim::FrameResults3D;
    using sdim::MeasurementArray2D;
    using sdim::MeasurementRecord;
    using sdim::Program;

    constexpr int kDimension = 3;
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

    void require(bool condition, const std::string &message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    int mod_qutrit(int value) {
        int r = value % kDimension;
        return r < 0 ? r + kDimension : r;
    }

    double seconds_between(TimePoint start, TimePoint end) {
        return std::chrono::duration<double>(end - start).count();
    }

    struct TimingSummary {
        double build_seconds = 0.0;
        double reference_seconds = 0.0;
        double ir_seconds = 0.0;
        double frame_seconds = 0.0;
        double postprocess_seconds = 0.0;
        double total_seconds = 0.0;
    };

    struct BenchmarkOptions {
        int shots = 10000;
        int distance_min = 3;
        int distance_max = 9;
        int num_p_points = 11;
        double p_min = 0.001;
        double p_max = 0.1;
        int logical_value = 0;
        std::uint64_t seed = 0;
        std::string output = "qutrit_repetition_memory_benchmark.csv";
    };

    BenchmarkOptions parse_options(int argc, char **argv) {
        BenchmarkOptions options;

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto require_value = [&](const std::string &name) -> const char * {
                if (i + 1 >= argc) {
                    throw std::runtime_error("Missing value after " + name);
                }
                return argv[++i];
            };

            if (arg == "--shots") {
                options.shots = std::stoi(require_value(arg));
            } else if (arg == "--distance-min" || arg == "--rounds-min") {
                options.distance_min = std::stoi(require_value(arg));
            } else if (arg == "--distance-max" || arg == "--rounds-max") {
                options.distance_max = std::stoi(require_value(arg));
            } else if (arg == "--num-p-points") {
                options.num_p_points = std::stoi(require_value(arg));
            } else if (arg == "--p-min") {
                options.p_min = std::stod(require_value(arg));
            } else if (arg == "--p-max") {
                options.p_max = std::stod(require_value(arg));
            } else if (arg == "--logical") {
                options.logical_value = std::stoi(require_value(arg));
            } else if (arg == "--seed") {
                options.seed = static_cast<std::uint64_t>(std::stoull(require_value(arg)));
            } else if (arg == "--output") {
                options.output = require_value(arg);
            } else {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        require(options.shots > 0, "shots must be positive");
        require(options.distance_min >= 2, "distance-min must be at least 2");
        require(options.distance_max >= options.distance_min, "distance-max must be >= distance-min");
        require(options.num_p_points > 0, "num-p-points must be positive");
        require(options.p_min > 0.0, "p-min must be positive for logarithmic spacing");
        require(options.p_max >= options.p_min, "p-max must be >= p-min");
        require(options.logical_value >= 0 && options.logical_value < kDimension, "logical must be 0, 1, or 2");

        return options;
    }

    std::vector<double> make_probability_grid(const BenchmarkOptions &options) {
        std::vector<double> values;
        values.reserve(static_cast<std::size_t>(options.num_p_points));

        if (options.num_p_points == 1) {
            values.push_back(options.p_min);
            return values;
        }

        const double log_min = std::log(options.p_min);
        const double log_max = std::log(options.p_max);
        for (int i = 0; i < options.num_p_points; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(options.num_p_points - 1);
            values.push_back(std::exp(log_min + t * (log_max - log_min)));
        }
        return values;
    }

    int first_odd_distance(int distance_min) {
        return distance_min % 2 == 0 ? distance_min + 1 : distance_min;
    }

    int ancilla_index(int distance, int edge) {
        return distance + edge;
    }

    Circuit make_memory_benchmark_circuit(int distance, int logical_value, double p_x) {
        const int rounds = distance;
        Circuit circuit(distance + distance - 1, kDimension);

        // Layout: data qudits occupy [0, distance), syndrome ancillas occupy
        // [distance, 2 * distance - 1). Edge e measures data[e] - data[e + 1].
        for (int q = 0; q < distance; ++q) {
            for (int k = 0; k < logical_value; ++k) {
                circuit.add_gate("X", q);
            }
        }

        const sdim::ParameterMap noise_params{
                {"prob_x", p_x},
                {"prob_z", 0.0},
                {"prob_xz", 0.0},
        };

        for (int round = 0; round < rounds; ++round) {
            for (int q = 0; q < distance; ++q) {
                circuit.add_gate("PAULI1", q, std::nullopt, noise_params);
            }

            for (int edge = 0; edge < distance - 1; ++edge) {
                const int ancilla = ancilla_index(distance, edge);
                circuit.add_gate("RESET", ancilla);
                circuit.add_gate("CNOT", edge, ancilla);
                circuit.add_gate("CNOT_INV", edge + 1, ancilla);
                circuit.add_gate("M", ancilla);
            }
        }

        std::vector<int> data_qudits;
        data_qudits.reserve(static_cast<std::size_t>(distance));
        for (int q = 0; q < distance; ++q) {
            data_qudits.push_back(q);
        }
        circuit.add_gate("M", data_qudits);
        return circuit;
    }

    const MeasurementRecord &measurement_at(const FrameResults3D &measurements,
                                            int qudit,
                                            std::size_t round,
                                            std::size_t shot) {
        require(qudit >= 0 && static_cast<std::size_t>(qudit) < measurements.qudits,
                "measurement qudit index out of range");
        require(round < measurements.rounds, "measurement round index out of range");
        require(shot < measurements.shots, "measurement shot index out of range");
        return measurements.at(static_cast<std::size_t>(qudit), round, shot);
    }

    int decode_by_maximum_likelihood(const FrameResults3D &measurements, int distance, std::size_t shot_index) {
        std::vector<int> counts(static_cast<std::size_t>(kDimension), 0);
        for (int q = 0; q < distance; ++q) {
            const MeasurementRecord &measurement = measurement_at(measurements, q, 0, shot_index);
            counts[static_cast<std::size_t>(mod_qutrit(static_cast<int>(measurement.measurement_value)))] += 1;
        }

        int best_symbol = -1;
        int best_count = -1;
        bool tied = false;
        for (int symbol = 0; symbol < kDimension; ++symbol) {
            const int count = counts[static_cast<std::size_t>(symbol)];
            if (count > best_count) {
                best_symbol = symbol;
                best_count = count;
                tied = false;
            } else if (count == best_count) {
                tied = true;
            }
        }
        return tied ? -1 : best_symbol;
    }

    void require_measurement_shape(const FrameResults3D &measurements, int distance) {
        require(measurements.qudits == static_cast<std::size_t>(distance + distance - 1),
                "measurement qudit count mismatch");

        const std::size_t measured_rounds = measurements.rounds;
        require(measured_rounds % 2 == 0, "ancilla measurements should be reset/syndrome pairs");
    }

    void print_usage() {
        std::cout
                << "Options:\n"
                << "  --shots <int>\n"
                << "  --distance-min <int>  (even values advance to the next odd distance)\n"
                << "  --distance-max <int>  (only odd distances are simulated)\n"
                << "  --num-p-points <int>\n"
                << "  --p-min <double>  (positive; probabilities are log-spaced)\n"
                << "  --p-max <double>\n"
                << "  --logical <0|1|2>\n"
                << "  --seed <uint64>\n"
                << "  --output <path>\n";
    }
} // namespace

int main(int argc, char **argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
                print_usage();
                return EXIT_SUCCESS;
            }
        }

        const BenchmarkOptions options = parse_options(argc, argv);
        const std::vector<double> probabilities = make_probability_grid(options);
        const std::filesystem::path output_path = std::filesystem::absolute(options.output);

        std::cout << "[info] qutrit repetition-code memory benchmark\n";
        std::cout << "[info] shots = " << options.shots
                  << ", logical = " << options.logical_value
                  << ", distance = [" << options.distance_min << ", " << options.distance_max << "]"
                  << " odd only"
                  << ", rounds = 1"
                  << ", log-spaced p in [" << options.p_min << ", " << options.p_max << "]"
                  << ", points = " << options.num_p_points << "\n";
        std::cout << "[info] writing CSV to " << output_path << "\n";

        std::ofstream out(output_path);
        require(static_cast<bool>(out), "failed to open output file: " + output_path.string());

        out << "p,distance,rounds,shots,failures,logical_error_rate,avg_ambiguous_decodes,avg_corrections,invalid_final_votes,"
            << "build_seconds,reference_seconds,ir_seconds,frame_seconds,postprocess_seconds,total_seconds,shots_per_second\n";
        out << std::setprecision(17);

        for (double p : probabilities) {
            for (int distance = first_odd_distance(options.distance_min);
                 distance <= options.distance_max;
                 distance += 2) {
                const int rounds = 1;
                TimingSummary timing;

                const TimePoint total_start = Clock::now();
                const TimePoint build_start = Clock::now();
                const Circuit circuit = make_memory_benchmark_circuit(distance, options.logical_value, p);
                Program program(circuit,
                                options.seed
                                + static_cast<std::uint64_t>(1000 * distance + std::llround(1e9 * p)));
                const TimePoint build_end = Clock::now();
                timing.build_seconds = seconds_between(build_start, build_end);

                // Use the frame simulator for noisy shots. The reference shot supplies
                // deterministic baseline measurements; build_ir samples PAULI1 errors,
                // and simulate_frame propagates those errors through the circuit IR.
                const TimePoint reference_start = Clock::now();
                (void) program.simulate();
                const MeasurementArray2D reference_results = Program::results_to_array(program.measurement_results());
                const TimePoint reference_end = Clock::now();
                timing.reference_seconds = seconds_between(reference_start, reference_end);

                const TimePoint ir_start = Clock::now();
                const BuildIrResult ir_result = program.build_ir(options.shots);
                const TimePoint ir_end = Clock::now();
                timing.ir_seconds = seconds_between(ir_start, ir_end);

                const TimePoint frame_start = Clock::now();
                const FrameResults3D measurements = Program::simulate_frame(
                        ir_result.ir,
                        reference_results,
                        circuit.num_qudits(),
                        circuit.dimension(),
                        options.shots,
                        ir_result.noise,
                        0);
                const TimePoint frame_end = Clock::now();
                timing.frame_seconds = seconds_between(frame_start, frame_end);

                int failures = 0;
                int invalid_final_votes = 0;
                long long total_ambiguous_decodes = 0;
                long long total_corrections = 0;

                const TimePoint postprocess_start = Clock::now();
                require_measurement_shape(measurements, distance);
                for (int shot = 0; shot < options.shots; ++shot) {
                    const int logical_out = decode_by_maximum_likelihood(
                            measurements,
                            distance,
                            static_cast<std::size_t>(shot));
                    const bool invalid_final_vote = logical_out < 0;
                    if (invalid_final_vote || logical_out != options.logical_value) {
                        failures += 1;
                    }
                    if (invalid_final_vote) {
                        invalid_final_votes += 1;
                        total_ambiguous_decodes += 1;
                    }
                }
                const TimePoint postprocess_end = Clock::now();
                timing.postprocess_seconds = seconds_between(postprocess_start, postprocess_end);
                timing.total_seconds = seconds_between(total_start, postprocess_end);

                const double logical_error_rate =
                        static_cast<double>(failures) / static_cast<double>(options.shots);
                const double avg_ambiguous_decodes =
                        static_cast<double>(total_ambiguous_decodes) / static_cast<double>(options.shots);
                const double avg_corrections =
                        static_cast<double>(total_corrections) / static_cast<double>(options.shots);
                const double shots_per_second =
                        timing.total_seconds > 0.0 ? static_cast<double>(options.shots) / timing.total_seconds : 0.0;

                out << p << ','
                    << distance << ','
                    << rounds << ','
                    << options.shots << ','
                    << failures << ','
                    << logical_error_rate << ','
                    << avg_ambiguous_decodes << ','
                    << avg_corrections << ','
                    << invalid_final_votes << ','
                    << timing.build_seconds << ','
                    << timing.reference_seconds << ','
                    << timing.ir_seconds << ','
                    << timing.frame_seconds << ','
                    << timing.postprocess_seconds << ','
                    << timing.total_seconds << ','
                    << shots_per_second << '\n';

                std::cout << "[pass] p=" << std::fixed << std::setprecision(6) << p
                          << " distance=" << distance
                          << " rounds=" << rounds
                          << " logical_error_rate=" << std::setprecision(8) << logical_error_rate
                          << " failures=" << failures << '/' << options.shots
                          << " total_seconds=" << std::setprecision(6) << timing.total_seconds
                          << " shots_per_second=" << std::setprecision(2) << shots_per_second << "\n";
            }
        }

        std::cout << "[pass] wrote memory benchmark CSV\n";
        return EXIT_SUCCESS;
    } catch (const std::exception &ex) {
        std::cerr << "[fail] repetition_memory_benchmark: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
