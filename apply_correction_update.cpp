#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

constexpr int Q = 3;

int mod_q(int x) {
    int r = x % Q;
    return r < 0 ? r + Q : r;
}

std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
}

struct DataKey {
    std::string p;
    int shot;
    int q;

    bool operator<(const DataKey &o) const {
        return std::tie(p, shot, q) < std::tie(o.p, o.shot, o.q);
    }
};

struct GroupKey {
    std::string p;
    int distance;
    int rounds;

    bool operator<(const GroupKey &o) const {
        return std::tie(p, distance, rounds) < std::tie(o.p, o.distance, o.rounds);
    }
};

struct FinalDataRow {
    std::string p;
    int distance;
    int rounds;
    int shot;
    int q;
    int value;
    int corrected_value;
};

struct DataCorrection {
    std::string p;
    int shot;
    int cluster_id;
    int charge;
    int data_qutrit;
    int t;
};

struct RawSummary {
    int shots = 0;
    int raw_majority_failures = 0;
    double raw_majority_logical_error_rate = 0.0;
};

std::map<DataKey, FinalDataRow> read_final_data(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Could not open " + filename);

    std::string line;
    std::getline(in, line);

    std::map<DataKey, FinalDataRow> data;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto c = split_csv_line(line);
        if (c.size() < 6) continue;

        FinalDataRow row;
        row.p = c[0];
        row.distance = std::stoi(c[1]);
        row.rounds = std::stoi(c[2]);
        row.shot = std::stoi(c[3]);
        row.q = std::stoi(c[4]);
        row.value = std::stoi(c[5]);
        row.corrected_value = row.value;

        data[{row.p, row.shot, row.q}] = row;
    }

    return data;
}

std::vector<DataCorrection> read_data_corrections(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Could not open " + filename);

    std::string line;
    std::getline(in, line);

    std::vector<DataCorrection> corrections;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto c = split_csv_line(line);
        if (c.size() < 10) continue;

        DataCorrection corr;
        corr.p = c[0];
        corr.shot = std::stoi(c[1]);
        corr.cluster_id = std::stoi(c[2]);
        corr.charge = std::stoi(c[3]);
        corr.data_qutrit = std::stoi(c[4]);
        corr.t = std::stoi(c[5]);

        corrections.push_back(corr);
    }

    return corrections;
}

std::map<GroupKey, RawSummary> read_summary(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Could not open " + filename);

    std::string line;
    std::getline(in, line);

    std::map<GroupKey, RawSummary> summary;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto c = split_csv_line(line);
        if (c.size() < 6) continue;

        GroupKey key{c[0], std::stoi(c[1]), std::stoi(c[2])};

        RawSummary s;
        s.shots = std::stoi(c[3]);
        s.raw_majority_failures = std::stoi(c[4]);
        s.raw_majority_logical_error_rate = std::stod(c[5]);

        summary[key] = s;
    }

    return summary;
}

void apply_data_corrections(
        std::map<DataKey, FinalDataRow> &data,
        const std::vector<DataCorrection> &corrections) {

    int applied = 0;
    int missed = 0;

    for (const DataCorrection &corr : corrections) {
        DataKey key{corr.p, corr.shot, corr.data_qutrit};

        auto it = data.find(key);
        if (it == data.end()) {
            missed++;
            continue;
        }

        it->second.corrected_value =
            mod_q(it->second.corrected_value - corr.charge);

        applied++;
    }

    std::cout << "applied data corrections = " << applied << "\n";
    std::cout << "missed data corrections = " << missed << "\n";
}

int majority_vote(const std::vector<int> &values) {
    std::vector<int> counts(Q, 0);

    for (int v : values) counts[mod_q(v)]++;

    int best_symbol = -1;
    int best_count = -1;
    bool tied = false;

    for (int s = 0; s < Q; ++s) {
        if (counts[s] > best_count) {
            best_symbol = s;
            best_count = counts[s];
            tied = false;
        } else if (counts[s] == best_count) {
            tied = true;
        }
    }

    return tied ? -1 : best_symbol;
}

int main(int argc, char **argv) {
    try {
        int logical_value = 0;
        if (argc >= 2) logical_value = std::stoi(argv[1]);

        auto data = read_final_data("final_data.csv");
        auto corrections = read_data_corrections("data_corrections.csv");
        auto raw_summary = read_summary("summary.csv");

        apply_data_corrections(data, corrections);

        std::ofstream corrected_data_out("corrected_final_data.csv");
        corrected_data_out << "p,distance,rounds,shot,q,value,corrected_value\n";

        std::map<GroupKey, std::map<int, std::vector<int>>> grouped;

        for (const auto &entry : data) {
            const FinalDataRow &r = entry.second;

            corrected_data_out << r.p << ','
                               << r.distance << ','
                               << r.rounds << ','
                               << r.shot << ','
                               << r.q << ','
                               << r.value << ','
                               << r.corrected_value << '\n';

            grouped[{r.p, r.distance, r.rounds}][r.shot]
                .push_back(r.corrected_value);
        }

        std::ofstream out("corrected_logical_summary.csv");

        out << "p,distance,rounds,shots,"
            << "raw_majority_failures,"
            << "raw_majority_logical_error_rate,"
            << "corrected_logical_failures,"
            << "corrected_logical_error_rate\n";

        for (const auto &entry : grouped) {
            const GroupKey &key = entry.first;
            const auto &shots = entry.second;

            int corrected_failures = 0;

            for (const auto &shot_entry : shots) {
                int logical_out = majority_vote(shot_entry.second);

                if (logical_out < 0 || logical_out != logical_value) {
                    corrected_failures++;
                }
            }

            int num_shots = static_cast<int>(shots.size());
            double corrected_rate =
                static_cast<double>(corrected_failures) /
                static_cast<double>(num_shots);

            RawSummary raw;
            auto it = raw_summary.find(key);
            if (it != raw_summary.end()) raw = it->second;

            out << key.p << ','
                << key.distance << ','
                << key.rounds << ','
                << num_shots << ','
                << raw.raw_majority_failures << ','
                << raw.raw_majority_logical_error_rate << ','
                << corrected_failures << ','
                << corrected_rate << '\n';
        }

        std::cout << "wrote corrected_final_data.csv\n";
        std::cout << "wrote corrected_logical_summary.csv\n";

        return 0;

    } catch (const std::exception &ex) {
        std::cerr << "[fail] " << ex.what() << "\n";
        return 1;
    }
}
