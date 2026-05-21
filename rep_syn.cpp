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
        int rounds = -1; // if negative, use rounds = distance
        int num_p_points = 11;
        double p_min = 0.001;
        double p_max = 0.1;
        int logical_value = 0;
        std::uint64_t seed = 0;
        std::string output = "qutrit_repetition_memory_benchmark.csv";
        std::string syndrome_output = "qutrit_repetition_syndrome_defects.csv";
        bool export_all_syndromes = false;
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
            } else if (arg == "--distance-min") {
                options.distance_min = std::stoi(require_value(arg));
            } else if (arg == "--distance-max") {
                options.distance_max = std::stoi(require_value(arg));
            } else if (arg == "--rounds") {
                options.rounds = std::stoi(require_value(arg));
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
            } else if (arg == "--syndrome-output") {
                options.syndrome_output = require_value(arg);
            } else if (arg == "--export-all-syndromes") {
                options.export_all_syndromes = true;
            } else {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        require(options.shots > 0, "shots must be positive");
        require(options.distance_min >= 2, "distance-min must be at least 2");
        require(options.distance_max >= options.distance_min, "distance-max must be >= distance-min");
        require(options.rounds != 0, "rounds must be positive, or omitted to use rounds = distance");
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

    Circuit make_memory_benchmark_circuit(int distance, int rounds, int logical_value, double p_x) {
        Circuit circuit(distance + distance - 1, kDimension);

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

    int measurement_value_at(const FrameResults3D &measurements,
                             int qudit,
                             std::size_t round,
                             std::size_t shot) {
        const MeasurementRecord &measurement = measurement_at(measurements, qudit, round, shot);
        return mod_qutrit(static_cast<int>(measurement.measurement_value));
    }

    std::size_t ancilla_syndrome_measurement_round(int circuit_round) {
        return static_cast<std::size_t>(2 * circuit_round + 1);
    }

    int decode_by_maximum_likelihood(const FrameResults3D &measurements, int distance, std::size_t shot_index) {
        std::vector<int> counts(static_cast<std::size_t>(kDimension), 0);
        for (int q = 0; q < distance; ++q) {
            const int value = measurement_value_at(measurements, q, 0, shot_index);
            counts[static_cast<std::size_t>(value)] += 1;
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

    void require_measurement_shape(const FrameResults3D &measurements, int distance, int rounds) {
        require(measurements.qudits == static_cast<std::size_t>(distance + distance - 1),
                "measurement qudit count mismatch");

        const std::size_t last_needed_ancilla_round =
                ancilla_syndrome_measurement_round(rounds - 1);
        require(measurements.rounds > last_needed_ancilla_round,
                "not enough measurement rounds to contain all ancilla syndromes");
    }

    int final_data_check_value(const FrameResults3D &measurements,
                               int edge,
                               std::size_t shot) {
        const int left = measurement_value_at(measurements, edge, 0, shot);
        const int right = measurement_value_at(measurements, edge + 1, 0, shot);
        return mod_qutrit(left - right);
    }

    void export_syndrome_defects(std::ofstream &syndrome_out,
                                 const FrameResults3D &measurements,
                                 double p,
                                 int distance,
                                 int rounds,
                                 int shot,
                                 bool export_all_syndromes) {
        for (int edge = 0; edge < distance - 1; ++edge) {
            int previous_syndrome = 0;

            for (int round = 0; round < rounds; ++round) {
                const int ancilla = ancilla_index(distance, edge);
                const std::size_t meas_round = ancilla_syndrome_measurement_round(round);

                const int syndrome = measurement_value_at(
                        measurements,
                        ancilla,
                        meas_round,
                        static_cast<std::size_t>(shot));

                const int defect = mod_qutrit(syndrome - previous_syndrome);

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

            const int final_syndrome =
                    final_data_check_value(
                            measurements,
                            edge,
                            static_cast<std::size_t>(shot));

            const int final_defect =
                    mod_qutrit(final_syndrome - previous_syndrome);

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

    void print_usage() {
        std::cout
                << "Options:\n"
                << "  --shots <int>\n"
                << "  --distance-min <int>\n"
                << "  --distance-max <int>\n"
                << "  --rounds <int>        optional; default is rounds = distance\n"
                << "  --num-p-points <int>\n"
                << "  --p-min <double>\n"
                << "  --p-max <double>\n"
                << "  --logical <0|1|2>\n"
                << "  --seed <uint64>\n"
                << "  --output <path>\n"
                << "  --syndrome-output <path>\n"
                << "  --export-all-syndromes\n";
    }
}

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

        const std::filesystem::path output_path =
                std::filesystem::absolute(options.output);

        const std::filesystem::path syndrome_output_path =
                std::filesystem::absolute(options.syndrome_output);

        std::cout << "[info] qutrit repetition-code memory benchmark + syndrome export\n";
        std::cout << "[info] shots = " << options.shots
                  << ", logical = " << options.logical_value
                  << ", distance = [" << options.distance_min << ", " << options.distance_max << "]"
                  << " odd only"
                  << ", rounds = "
                  << (options.rounds > 0 ? std::to_string(options.rounds) : std::string("distance"))
                  << ", p in [" << options.p_min << ", " << options.p_max << "]"
                  << ", points = " << options.num_p_points << "\n";

        std::cout << "[info] writing summary CSV to " << output_path << "\n";
        std::cout << "[info] writing syndrome-defect CSV to " << syndrome_output_path << "\n";

        std::ofstream out(output_path);
        require(static_cast<bool>(out),
                "failed to open output file: " + output_path.string());

        std::ofstream syndrome_out(syndrome_output_path);
        require(static_cast<bool>(syndrome_out),
                "failed to open syndrome output file: " + syndrome_output_path.string());

        out << "p,distance,rounds,shots,failures,logical_error_rate,"
            << "avg_ambiguous_decodes,avg_corrections,invalid_final_votes,"
            << "build_seconds,reference_seconds,ir_seconds,frame_seconds,"
            << "postprocess_seconds,total_seconds,shots_per_second\n";
        out << std::setprecision(17);

        syndrome_out << "p,distance,rounds,shot,round,edge,"
                     << "is_final_boundary,syndrome,previous_syndrome,defect\n";
        syndrome_out << std::setprecision(17);

        for (double p : probabilities) {
            for (int distance = first_odd_distance(options.distance_min);
                 distance <= options.distance_max;
                 distance += 2) {
                const int rounds =
                        options.rounds > 0 ? options.rounds : distance;

                TimingSummary timing;

                const TimePoint total_start = Clock::now();

                const TimePoint build_start = Clock::now();
                const Circuit circuit =
                        make_memory_benchmark_circuit(
                                distance,
                                rounds,
                                options.logical_value,
                                p);

                Program program(
                        circuit,
                        options.seed
                        + static_cast<std::uint64_t>(
                                1000 * distance
                                + 100000 * rounds
                                + std::llround(1e9 * p)));

                const TimePoint build_end = Clock::now();
                timing.build_seconds =
                        seconds_between(build_start, build_end);

                const TimePoint reference_start = Clock::now();
                (void) program.simulate();

                const MeasurementArray2D reference_results =
                        Program::results_to_array(
                                program.measurement_results());

                const TimePoint reference_end = Clock::now();
                timing.reference_seconds =
                        seconds_between(reference_start, reference_end);

                const TimePoint ir_start = Clock::now();
                const BuildIrResult ir_result =
                        program.build_ir(options.shots);
                const TimePoint ir_end = Clock::now();
                timing.ir_seconds =
                        seconds_between(ir_start, ir_end);

                const TimePoint frame_start = Clock::now();
                const FrameResults3D measurements =
                        Program::simulate_frame(
                                ir_result.ir,
                                reference_results,
                                circuit.num_qudits(),
                                circuit.dimension(),
                                options.shots,
                                ir_result.noise,
                                0);
                const TimePoint frame_end = Clock::now();
                timing.frame_seconds =
                        seconds_between(frame_start, frame_end);

                int failures = 0;
                int invalid_final_votes = 0;
                long long total_ambiguous_decodes = 0;
                long long total_corrections = 0;

                const TimePoint postprocess_start = Clock::now();

                require_measurement_shape(
                        measurements,
                        distance,
                        rounds);

                for (int shot = 0; shot < options.shots; ++shot) {
                    export_syndrome_defects(
                            syndrome_out,
                            measurements,
                            p,
                            distance,
                            rounds,
                            shot,
                            options.export_all_syndromes);

                    const int logical_out =
                            decode_by_maximum_likelihood(
                                    measurements,
                                    distance,
                                    static_cast<std::size_t>(shot));

                    const bool invalid_final_vote =
                            logical_out < 0;

                    if (invalid_final_vote
                        || logical_out != options.logical_value) {
                        failures += 1;
                    }

                    if (invalid_final_vote) {
                        invalid_final_votes += 1;
                        total_ambiguous_decodes += 1;
                    }
                }

                const TimePoint postprocess_end = Clock::now();

                timing.postprocess_seconds =
                        seconds_between(postprocess_start, postprocess_end);
                timing.total_seconds =
                        seconds_between(total_start, postprocess_end);

                const double logical_error_rate =
                        static_cast<double>(failures)
                        / static_cast<double>(options.shots);

                const double avg_ambiguous_decodes =
                        static_cast<double>(total_ambiguous_decodes)
                        / static_cast<double>(options.shots);

                const double avg_corrections =
                        static_cast<double>(total_corrections)
                        / static_cast<double>(options.shots);

                const double shots_per_second =
                        timing.total_seconds > 0.0
                        ? static_cast<double>(options.shots)
                          / timing.total_seconds
                        : 0.0;

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
                          << " logical_error_rate="
                          << std::setprecision(8) << logical_error_rate
                          << " failures=" << failures << '/'
                          << options.shots
                          << " total_seconds="
                          << std::setprecision(6) << timing.total_seconds
                          << " shots_per_second="
                          << std::setprecision(2) << shots_per_second
                          << "\n";
            }
        }

        std::cout << "[pass] wrote summary benchmark CSV\n";
        std::cout << "[pass] wrote syndrome-defect CSV\n";

        return EXIT_SUCCESS;
    } catch (const std::exception &ex) {
        std::cerr << "[fail] repetition_memory_benchmark: "
                  << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
