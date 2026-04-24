
#include <algorithm>
#include <cctype>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <optional>
#include <variant>

#include "gatedata.h"

namespace sdim {
namespace {

ParameterMap make_n1_defaults() {
    return ParameterMap{{"channel", std::string("d")}, {"prob", 0.01}};
}

ParameterMap make_n2_defaults(int d) {
    const std::size_t size = static_cast<std::size_t>(d) * d * d * d;
    std::vector<double> dist(size, 1.0 / static_cast<double>(size));
    return ParameterMap{{"prob_dist", std::move(dist)}};
}

ParameterMap make_dep2_defaults() {
    return ParameterMap{{"prob", 0.01}};
}

ParameterMap make_pauli1_defaults() {
    return ParameterMap{
        {"prob_x", 0.01},
        {"prob_z", 0.01},
        {"prob_xz", 0.01}
    };
}

}  // namespace

GateData::GateData(int dimension) : dimension_(dimension) {
    if (dimension_ < 2) {
        throw std::invalid_argument("Dimension must be greater than 1");
    }
    add_gate_data_pauli(dimension_);
    add_gate_hada(dimension_);
    add_gate_controlled(dimension_);
    add_gate_collapsing(dimension_);
    add_gate_noise(dimension_);
}

std::string GateData::to_upper_ascii(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return s;
}

void GateData::add_gate(const std::string &name, int arg_count) {
    const std::string canonical_name = to_upper_ascii(name);
    const int gate_id = num_gates_;
    gate_map_[canonical_name] = Gate(canonical_name, arg_count, gate_id);
    gate_names_by_id_.push_back(canonical_name);
    ++num_gates_;
}

void GateData::add_gate_alias(const std::string &name, const std::vector<std::string> &aliases) {
    const std::string canonical_name = to_upper_ascii(name);
    if (gate_map_.find(canonical_name) == gate_map_.end()) {
        throw std::invalid_argument("Cannot add aliases for unknown gate: " + canonical_name);
    }
    for (const auto &alias : aliases) {
        alias_map_[to_upper_ascii(alias)] = canonical_name;
    }
}

std::optional<int> GateData::get_gate_id(const std::string &gate_name) const {
    const Gate *gate = find_gate(gate_name);
    if (!gate) {
        return std::nullopt;
    }
    return gate->gate_id;
}

const std::string &GateData::get_gate_name(int gate_id) const {
    if (gate_id < 0 || gate_id >= static_cast<int>(gate_names_by_id_.size())) {
        throw std::out_of_range("Gate ID not found: " + std::to_string(gate_id));
    }
    return gate_names_by_id_[static_cast<std::size_t>(gate_id)];
}

const Gate *GateData::find_gate(const std::string &gate_name_or_alias) const {
    std::string key = to_upper_ascii(gate_name_or_alias);

    auto gate_it = gate_map_.find(key);
    if (gate_it != gate_map_.end()) {
        return &gate_it->second;
    }

    auto alias_it = alias_map_.find(key);
    if (alias_it == alias_map_.end()) {
        return nullptr;
    }

    gate_it = gate_map_.find(alias_it->second);
    if (gate_it == gate_map_.end()) {
        return nullptr;
    }
    return &gate_it->second;
}

const Gate &GateData::get_gate(const std::string &gate_name_or_alias) const {
    const Gate *gate = find_gate(gate_name_or_alias);
    if (!gate) {
        throw std::invalid_argument("Gate not found: " + gate_name_or_alias);
    }
    return *gate;
}

void GateData::add_gate_data_pauli(int) {
    add_gate("I", 1);
    add_gate("X", 1);
    add_gate("X_INV", 1);
    add_gate("Z", 1);
    add_gate("Z_INV", 1);
}

void GateData::add_gate_hada(int) {
    add_gate("H", 1);
    add_gate_alias("H", {"R", "DFT"});

    add_gate("H_INV", 1);
    add_gate_alias("H_INV", {"R_INV", "DFT_INV", "H_DAG", "R_DAG", "DFT_DAG"});

    add_gate("P", 1);
    add_gate_alias("P", {"PHASE", "S"});

    add_gate("P_INV", 1);
    add_gate_alias("P_INV", {"PHASE_INV", "S_INV"});
}

void GateData::add_gate_controlled(int) {
    add_gate("CNOT", 2);
    add_gate_alias("CNOT", {"SUM", "CX", "C"});

    add_gate("CNOT_INV", 2);
    add_gate_alias("CNOT_INV", {"SUM_INV", "CX_INV", "C_INV"});

    add_gate("CZ", 2);
    add_gate("CZ_INV", 2);
    add_gate("SWAP", 2);
}

void GateData::add_gate_collapsing(int) {
    add_gate("M", 1);
    add_gate_alias("M", {"MEASURE", "COLLAPSE", "MZ"});

    add_gate("M_X", 1);
    add_gate_alias("M_X", {"MEASURE_X", "MX"});

    add_gate("RESET", 1);
    add_gate_alias("RESET", {"MR", "MEASURE_RESET", "MEASURE_R"});
}

void GateData::add_gate_noise(int d) {
    add_gate("N1", 1);
    add_gate_alias("N1", {"NOISE1"});
    gate_map_.at("N1").defaults = make_n1_defaults();

    add_gate("N2", 2);
    add_gate_alias("N2", {"NOISE2"});
    gate_map_.at("N2").defaults = make_n2_defaults(d);

    add_gate("DEP2", 2);
    add_gate_alias("DEP2", {"DEPOLARIZING_2"});
    gate_map_.at("DEP2").defaults = make_dep2_defaults();

    add_gate("PAULI1", 1);
    add_gate_alias("PAULI1", {"PAULI_NOISE_1", "PAULI_CHANNEL_1", "PAULI_1"});
    gate_map_.at("PAULI1").defaults = make_pauli1_defaults();
}

std::ostream &operator<<(std::ostream &os, const Gate &gate) {
    os << gate.name << ' ' << gate.gate_id;
    return os;
}

std::string parameter_value_to_string(const ParameterValue &value) {
    struct Visitor {
        std::string operator()(std::monostate) const { return "null"; }
        std::string operator()(std::int64_t v) const { return std::to_string(v); }

        std::string operator()(double v) const {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        }

        std::string operator()(const std::string &v) const {
            return '"' + v + '"';
        }

        std::string operator()(const std::vector<double> &v) const {
            std::ostringstream oss;
            oss << '[';
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (i > 0) {
                    oss << ", ";
                }
                oss << v[i];
            }
            oss << ']';
            return oss.str();
        }
    };

    return std::visit(Visitor{}, value);
}

std::ostream &operator<<(std::ostream &os, const GateData &gate_data) {
    bool first = true;
    for (const auto &[name, gate] : gate_data.gate_map()) {
        if (!first) {
            os << '\n';
        }
        os << gate;
        first = false;
    }
    return os;
}

}  // namespace sdim
