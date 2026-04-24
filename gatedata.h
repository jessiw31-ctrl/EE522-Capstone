#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace sdim {

using ParameterValue = std::variant<
    std::monostate,
    std::int64_t,
    double,
    std::string,
    std::vector<double>
>;

using ParameterMap = std::map<std::string, ParameterValue>;

struct Gate {
    std::string name;
    int arg_count;
    int gate_id;
    std::optional<ParameterMap> defaults;

    Gate() = default;

    Gate(std::string name_, int arg_count_, int gate_id_)
        : name(std::move(name_)),
          arg_count(arg_count_),
          gate_id(gate_id_) {}
};

class GateData {
public:
    explicit GateData(int dimension = 2);

    void add_gate(const std::string &name, int arg_count);

    void add_gate_alias(
        const std::string &name,
        const std::vector<std::string> &aliases
    );

    std::optional<int> get_gate_id(const std::string &gate_name) const;

    const std::string &get_gate_name(int gate_id) const;

    const Gate *find_gate(const std::string &gate_name_or_alias) const;

    const Gate &get_gate(const std::string &gate_name_or_alias) const;

    int dimension() const noexcept { return dimension_; }

    int num_gates() const noexcept { return num_gates_; }

    const std::map<std::string, Gate> &gate_map() const noexcept {
        return gate_map_;
    }

    const std::map<std::string, std::string> &alias_map() const noexcept {
        return alias_map_;
    }

private:
    std::map<std::string, Gate> gate_map_;
    std::map<std::string, std::string> alias_map_;
    std::vector<std::string> gate_names_by_id_;

    int num_gates_ = 0;
    int dimension_ = 2;

    static std::string to_upper_ascii(std::string s);

    void add_gate_data_pauli(int d);
    void add_gate_hada(int d);
    void add_gate_controlled(int d);
    void add_gate_collapsing(int d);
    void add_gate_noise(int d);
};

std::ostream &operator<<(std::ostream &os, const Gate &gate);

std::ostream &operator<<(std::ostream &os, const GateData &gate_data);

std::string parameter_value_to_string(const ParameterValue &value);

}  // namespace sdim
