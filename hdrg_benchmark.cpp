#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

constexpr int Q = 3;

constexpr int WX = 1;
constexpr int WT = 4;

int mod_q(int x) {
    int r = x % Q;
    return r < 0 ? r + Q : r;
}

struct ShotKey {
    std::string p;
    int shot;

    bool operator<(const ShotKey &other) const {
        if (p != other.p) return p < other.p;
        return shot < other.shot;
    }
};

struct Defect {
    int id = -1;
    std::string p;
    int shot = -1;
    int x = 0;
    int t = 0;
    int charge = 0;
    bool final_boundary = false;
};

struct Cluster {
    int id = -1;
    std::vector<int> defect_indices;
    int total_charge_mod3 = 0;
    bool touches_left_boundary = false;
    bool touches_right_boundary = false;
    bool touches_initial_time_boundary = false;
    bool touches_final_time_boundary = false;
    bool neutral = false;
    int radius = 0;
};

struct CorrectionChain {
    int cluster_id = -1;
    int defect_a = -1;
    int defect_b = -1;
    int charge = 0;
    bool to_boundary = false;
    std::string boundary_name;
    bool is_temporal_boundary = false;
    std::vector<std::pair<int, int>> path;
};

struct HDRGResult {
    bool success = false;
    int final_radius = 0;
    std::vector<Cluster> clusters;
    std::vector<CorrectionChain> chains;
};

std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
}

std::map<ShotKey, std::vector<Defect>>
read_syndrome_csv(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Could not open " + filename);

    std::string header;
    std::getline(in, header);

    std::map<ShotKey, std::vector<Defect>> shots;
    std::string line;
    int global_id = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        auto cols = split_csv_line(line);
        if (cols.size() < 10) continue;

        std::string p = cols[0];
        int shot = std::stoi(cols[3]);
        int round = std::stoi(cols[4]);
        int edge = std::stoi(cols[5]);
        bool final_boundary = std::stoi(cols[6]) != 0;
        int defect = std::stoi(cols[9]);

        ShotKey key{p, shot};
        shots[key];

        defect = mod_q(defect);
        if (defect == 0) continue;

        Defect d;
        d.id = global_id++;
        d.p = p;
        d.shot = shot;
        d.x = edge;
        d.t = round;
        d.charge = defect;
        d.final_boundary = final_boundary;

        shots[key].push_back(d);
    }

    return shots;
}

int spacetime_distance(const Defect &a, const Defect &b) {
    return WX * std::abs(a.x - b.x)
         + WT * std::abs(a.t - b.t);
}

int distance_to_left_boundary(const Defect &d) {
    return WX * (d.x + 1);
}

int distance_to_right_boundary(const Defect &d, int distance) {
    return WX * ((distance - 1) - d.x);
}

int distance_to_initial_time_boundary(const Defect &d) {
    return WT * d.t;
}

int distance_to_final_time_boundary(const Defect &d, int rounds) {
    return WT * (rounds - d.t);
}

std::vector<std::vector<int>> connected_components(
        const std::vector<Defect> &defects,
        int radius) {
    int n = static_cast<int>(defects.size());
    std::vector<std::vector<int>> adj(n);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (spacetime_distance(defects[i], defects[j]) <= radius) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    std::vector<bool> visited(n, false);
    std::vector<std::vector<int>> comps;

    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;

        std::vector<int> comp;
        std::queue<int> q;
        q.push(i);
        visited[i] = true;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            comp.push_back(u);

            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        comps.push_back(comp);
    }

    return comps;
}

Cluster make_cluster(
        int cluster_id,
        const std::vector<int> &component,
        const std::vector<Defect> &defects,
        int radius,
        int distance,
        int rounds) {
    int charge = 0;

    bool touches_left = false;
    bool touches_right = false;
    bool touches_initial_time = false;
    bool touches_final_time = false;

    for (int idx : component) {
        const Defect &d = defects[idx];
        charge += d.charge;

        if (distance_to_left_boundary(d) <= radius) touches_left = true;
        if (distance_to_right_boundary(d, distance) <= radius) touches_right = true;
        if (distance_to_initial_time_boundary(d) <= radius) touches_initial_time = true;

        if (distance_to_final_time_boundary(d, rounds) <= radius || d.final_boundary) {
            touches_final_time = true;
        }
    }

    Cluster c;
    c.id = cluster_id;
    c.defect_indices = component;
    c.total_charge_mod3 = mod_q(charge);
    c.touches_left_boundary = touches_left;
    c.touches_right_boundary = touches_right;
    c.touches_initial_time_boundary = touches_initial_time;
    c.touches_final_time_boundary = touches_final_time;
    c.radius = radius;

    c.neutral =
        (c.total_charge_mod3 == 0) ||
        c.touches_left_boundary ||
        c.touches_right_boundary ||
        c.touches_initial_time_boundary ||
        c.touches_final_time_boundary;

    return c;
}

std::vector<std::pair<int, int>> manhattan_path_between_points(
        int x1, int t1, int x2, int t2) {
    std::vector<std::pair<int, int>> path;
    int x = x1;
    int t = t1;
    path.push_back({x, t});

    while (x != x2) {
        x += (x2 > x) ? 1 : -1;
        path.push_back({x, t});
    }

    while (t != t2) {
        t += (t2 > t) ? 1 : -1;
        path.push_back({x, t});
    }

    return path;
}

struct TreeNode {
    bool is_boundary = false;
    int defect_idx = -1;
    std::string boundary_name;
    int x = 0;
    int t = 0;
    int charge = 0;
};

int node_distance(const TreeNode &a, const TreeNode &b) {
    return WX * std::abs(a.x - b.x)
         + WT * std::abs(a.t - b.t);
}

std::vector<int> build_mst_parent(const std::vector<TreeNode> &nodes) {
    int n = static_cast<int>(nodes.size());

    std::vector<int> parent(n, -1);
    std::vector<int> best_dist(n, 1000000000);
    std::vector<bool> used(n, false);

    best_dist[0] = 0;

    for (int it = 0; it < n; ++it) {
        int v = -1;

        for (int i = 0; i < n; ++i) {
            if (!used[i] && (v == -1 || best_dist[i] < best_dist[v])) {
                v = i;
            }
        }

        if (v == -1) break;
        used[v] = true;

        for (int u = 0; u < n; ++u) {
            if (used[u]) continue;

            int dist = node_distance(nodes[v], nodes[u]);
            if (dist < best_dist[u]) {
                best_dist[u] = dist;
                parent[u] = v;
            }
        }
    }

    return parent;
}

void choose_best_boundary_node(
        const Cluster &cluster,
        const std::vector<Defect> &defects,
        int distance,
        int rounds,
        TreeNode &boundary,
        bool &found) {
    int best_dist = 1000000000;
    found = false;

    for (int idx : cluster.defect_indices) {
        const Defect &d = defects[idx];

        auto try_boundary = [&](bool allowed,
                                const std::string &name,
                                int bx,
                                int bt,
                                int dist) {
            if (allowed && dist < best_dist) {
                best_dist = dist;
                found = true;
                boundary.is_boundary = true;
                boundary.defect_idx = -1;
                boundary.boundary_name = name;
                boundary.x = bx;
                boundary.t = bt;
                boundary.charge = 0;
            }
        };

        try_boundary(cluster.touches_left_boundary,
                     "L", -1, d.t, distance_to_left_boundary(d));

        try_boundary(cluster.touches_right_boundary,
                     "R", distance - 1, d.t, distance_to_right_boundary(d, distance));

        try_boundary(cluster.touches_initial_time_boundary,
                     "T0", d.x, 0, distance_to_initial_time_boundary(d));

        try_boundary(cluster.touches_final_time_boundary,
                     "TF", d.x, rounds, distance_to_final_time_boundary(d, rounds));
    }
}

std::vector<TreeNode> make_cluster_nodes(
        const Cluster &cluster,
        const std::vector<Defect> &defects,
        int distance,
        int rounds) {
    std::vector<TreeNode> nodes;

    if (cluster.total_charge_mod3 != 0) {
        TreeNode b;
        bool found = false;
        choose_best_boundary_node(cluster, defects, distance, rounds, b, found);
        if (found) nodes.push_back(b);
    }

    for (int idx : cluster.defect_indices) {
        const Defect &d = defects[idx];

        TreeNode n;
        n.is_boundary = false;
        n.defect_idx = idx;
        n.boundary_name = "";
        n.x = d.x;
        n.t = d.t;
        n.charge = d.charge;

        nodes.push_back(n);
    }

    return nodes;
}

bool is_temporal_boundary_name(const std::string &name) {
    return name == "T0" || name == "TF";
}

int compute_subtree_charge_and_emit(
        int u,
        const std::vector<std::vector<int>> &children,
        const std::vector<TreeNode> &nodes,
        const std::vector<Defect> &defects,
        const Cluster &cluster,
        std::vector<CorrectionChain> &chains) {
    int subtree_charge = nodes[u].charge;

    for (int v : children[u]) {
        int child_charge = compute_subtree_charge_and_emit(
            v, children, nodes, defects, cluster, chains);

        subtree_charge = mod_q(subtree_charge + child_charge);

        if (child_charge != 0) {
            const TreeNode &child = nodes[v];
            const TreeNode &par = nodes[u];

            CorrectionChain chain;
            chain.cluster_id = cluster.id;
            chain.charge = child_charge;
            chain.path = manhattan_path_between_points(child.x, child.t, par.x, par.t);

            if (!child.is_boundary) {
                chain.defect_a = defects[child.defect_idx].id;
            }

            if (par.is_boundary) {
                chain.defect_b = -1;
                chain.to_boundary = true;
                chain.boundary_name = par.boundary_name;
                chain.is_temporal_boundary = is_temporal_boundary_name(par.boundary_name);
            } else {
                chain.defect_b = defects[par.defect_idx].id;
                chain.to_boundary = false;
                chain.boundary_name = "";
                chain.is_temporal_boundary = false;
            }

            chains.push_back(chain);
        }
    }

    return mod_q(subtree_charge);
}

std::vector<CorrectionChain> make_correction_chains(
        const std::vector<Cluster> &clusters,
        const std::vector<Defect> &defects,
        int distance,
        int rounds) {
    std::vector<CorrectionChain> chains;

    for (const Cluster &cluster : clusters) {
        if (!cluster.neutral) continue;
        if (cluster.defect_indices.empty()) continue;

        std::vector<TreeNode> nodes =
            make_cluster_nodes(cluster, defects, distance, rounds);

        if (nodes.empty()) continue;

        std::vector<int> parent = build_mst_parent(nodes);
        int n = static_cast<int>(nodes.size());

        std::vector<std::vector<int>> children(n);
        for (int i = 1; i < n; ++i) {
            if (parent[i] >= 0) {
                children[parent[i]].push_back(i);
            }
        }

        int leftover_charge = compute_subtree_charge_and_emit(
            0, children, nodes, defects, cluster, chains);

        if (leftover_charge != 0 && !nodes[0].is_boundary) {
            std::cerr << "Warning: neutral cluster left nonzero charge "
                      << leftover_charge << " in cluster " << cluster.id << "\n";
        }
    }

    return chains;
}

HDRGResult hdrg_decode_one_shot(
        const std::vector<Defect> &defects,
        int distance,
        int rounds) {
    HDRGResult result;

    if (defects.empty()) {
        result.success = true;
        result.final_radius = 0;
        return result;
    }

    int max_radius = WX * distance + WT * rounds + 2;

    for (int radius = 1; radius <= max_radius; ++radius) {
        auto comps = connected_components(defects, radius);

        std::vector<Cluster> clusters;
        bool all_neutral = true;

        for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
            Cluster c = make_cluster(i, comps[i], defects, radius, distance, rounds);
            if (!c.neutral) all_neutral = false;
            clusters.push_back(c);
        }

        if (all_neutral) {
            result.success = true;
            result.final_radius = radius;
            result.clusters = clusters;
            result.chains = make_correction_chains(clusters, defects, distance, rounds);
            return result;
        }
    }

    auto comps = connected_components(defects, max_radius);

    std::vector<Cluster> clusters;
    bool all_neutral = true;

    for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
        Cluster c = make_cluster(i, comps[i], defects, max_radius, distance, rounds);
        if (!c.neutral) all_neutral = false;
        clusters.push_back(c);
    }

    result.success = all_neutral;
    result.final_radius = max_radius;
    result.clusters = clusters;
    result.chains = make_correction_chains(clusters, defects, distance, rounds);

    return result;
}

void write_clusters_csv(
        std::ofstream &out,
        const ShotKey &key,
        const HDRGResult &result,
        const std::vector<Defect> &defects) {
    for (const Cluster &c : result.clusters) {
        for (int idx : c.defect_indices) {
            const Defect &d = defects[idx];

            out << key.p << ','
                << key.shot << ','
                << c.id << ','
                << c.radius << ','
                << c.total_charge_mod3 << ','
                << c.touches_left_boundary << ','
                << c.touches_right_boundary << ','
                << c.touches_initial_time_boundary << ','
                << c.touches_final_time_boundary << ','
                << c.neutral << ','
                << d.id << ','
                << d.x << ','
                << d.t << ','
                << d.charge << ','
                << d.final_boundary << '\n';
        }
    }
}

void write_correction_chains_csv(
        std::ofstream &out,
        const ShotKey &key,
        const HDRGResult &result) {
    for (const CorrectionChain &chain : result.chains) {
        for (int step = 0; step < static_cast<int>(chain.path.size()); ++step) {
            out << key.p << ','
                << key.shot << ','
                << chain.cluster_id << ','
                << chain.defect_a << ','
                << chain.defect_b << ','
                << chain.charge << ','
                << chain.to_boundary << ','
                << chain.boundary_name << ','
                << chain.is_temporal_boundary << ','
                << step << ','
                << chain.path[step].first << ','
                << chain.path[step].second << '\n';
        }
    }
}

void write_data_corrections_csv(
        std::ofstream &out,
        const ShotKey &key,
        const HDRGResult &result,
        int distance) {
    for (const CorrectionChain &chain : result.chains) {
        for (int step = 1; step < static_cast<int>(chain.path.size()); ++step) {
            int x0 = chain.path[step - 1].first;
            int t0 = chain.path[step - 1].second;
            int x1 = chain.path[step].first;
            int t1 = chain.path[step].second;

            if (t0 != t1) continue;
            if (x0 == x1) continue;

            int right = std::max(x0, x1);
            int data_qutrit = right;

            if (data_qutrit < 0 || data_qutrit >= distance) continue;

            out << key.p << ','
                << key.shot << ','
                << chain.cluster_id << ','
                << chain.charge << ','
                << data_qutrit << ','
                << t0 << ','
                << x0 << ','
                << t0 << ','
                << x1 << ','
                << t1 << '\n';
        }
    }
}

int main(int argc, char **argv) {
    std::string syndrome_file = "syndrome.csv";
    int distance = 3;
    int rounds = 3;

    if (argc >= 2) syndrome_file = argv[1];
    if (argc >= 3) distance = std::stoi(argv[2]);
    if (argc >= 4) rounds = std::stoi(argv[3]);

    auto shots = read_syndrome_csv(syndrome_file);

    std::ofstream summary("hdrg_summary.csv");
    std::ofstream clusters_out("clusters.csv");
    std::ofstream chains_out("correction_chains.csv");
    std::ofstream data_corr_out("data_corrections.csv");

    summary << "p,shot,d,rounds,num_defects,success,final_radius,num_clusters,num_chains,runtime_ms\n";

    clusters_out
        << "p,shot,cluster_id,radius,total_charge_mod3,"
        << "touches_left_boundary,touches_right_boundary,"
        << "touches_initial_time_boundary,touches_final_time_boundary,"
        << "neutral,defect_id,x,t,charge,final_boundary\n";

    chains_out
        << "p,shot,cluster_id,defect_a,defect_b,charge,"
        << "to_boundary,boundary_name,is_temporal_boundary,path_step,x,t\n";

    data_corr_out
        << "p,shot,cluster_id,charge,data_qutrit,t,x0,t0,x1,t1\n";

    int failures = 0;
    double total_runtime_ms = 0.0;

    std::map<std::string, double> runtime_by_p;
    std::map<std::string, int> shots_by_p;
    std::map<std::string, int> failures_by_p;
    std::map<std::string, int> total_defects_by_p;

    auto total_start = std::chrono::high_resolution_clock::now();

    for (const auto &entry : shots) {
        const ShotKey &key = entry.first;
        const std::vector<Defect> &defects = entry.second;

        auto start = std::chrono::high_resolution_clock::now();

        HDRGResult result = hdrg_decode_one_shot(defects, distance, rounds);

        auto end = std::chrono::high_resolution_clock::now();

        double runtime_ms =
            std::chrono::duration<double, std::milli>(end - start).count();

        total_runtime_ms += runtime_ms;

        runtime_by_p[key.p] += runtime_ms;
        shots_by_p[key.p] += 1;
        total_defects_by_p[key.p] += static_cast<int>(defects.size());

        if (!result.success) {
            failures += 1;
            failures_by_p[key.p] += 1;
        }

        summary << key.p << ','
                << key.shot << ','
                << distance << ','
                << rounds << ','
                << defects.size() << ','
                << result.success << ','
                << result.final_radius << ','
                << result.clusters.size() << ','
                << result.chains.size() << ','
                << runtime_ms << '\n';

        write_clusters_csv(clusters_out, key, result, defects);
        write_correction_chains_csv(chains_out, key, result);
        write_data_corrections_csv(data_corr_out, key, result, distance);
    }

    auto total_end = std::chrono::high_resolution_clock::now();

    double wall_runtime_ms =
        std::chrono::duration<double, std::milli>(total_end - total_start).count();

    double failure_rate =
        !shots.empty()
        ? static_cast<double>(failures) / static_cast<double>(shots.size())
        : 0.0;

    double avg_runtime_ms =
        !shots.empty()
        ? total_runtime_ms / static_cast<double>(shots.size())
        : 0.0;

    std::ofstream benchmark_out("benchmark_summary.csv");

    benchmark_out
        << "p,d,rounds,shots,avg_runtime_ms,total_runtime_ms,"
        << "avg_num_defects,failures,failure_rate\n";

    for (const auto &entry : runtime_by_p) {
        const std::string &p_value = entry.first;

        double p_total_runtime_ms = entry.second;
        int p_shots = shots_by_p[p_value];
        int p_failures = failures_by_p[p_value];
        int p_total_defects = total_defects_by_p[p_value];

        double p_avg_runtime_ms =
            p_shots > 0
            ? p_total_runtime_ms / static_cast<double>(p_shots)
            : 0.0;

        double p_avg_num_defects =
            p_shots > 0
            ? static_cast<double>(p_total_defects) / static_cast<double>(p_shots)
            : 0.0;

        double p_failure_rate =
            p_shots > 0
            ? static_cast<double>(p_failures) / static_cast<double>(p_shots)
            : 0.0;

        benchmark_out << p_value << ','
                      << distance << ','
                      << rounds << ','
                      << p_shots << ','
                      << p_avg_runtime_ms << ','
                      << p_total_runtime_ms << ','
                      << p_avg_num_defects << ','
                      << p_failures << ','
                      << p_failure_rate << '\n';
    }

    std::cout << "HDRG decode complete\n";
    std::cout << "WX = " << WX << ", WT = " << WT << "\n";
    std::cout << "decoded shots = " << shots.size() << "\n";
    std::cout << "cluster-neutralization failures = " << failures << "\n";
    std::cout << "cluster-neutralization failure rate = " << failure_rate << "\n";
    std::cout << "total decoder runtime from per-shot sum = "
              << total_runtime_ms << " ms\n";
    std::cout << "average decoder runtime per shot = "
              << avg_runtime_ms << " ms\n";
    std::cout << "total wall runtime for decoding loop = "
              << wall_runtime_ms << " ms\n";
    std::cout << "wrote hdrg_summary.csv\n";
    std::cout << "wrote benchmark_summary.csv\n";
    std::cout << "wrote clusters.csv\n";
    std::cout << "wrote correction_chains.csv\n";
    std::cout << "wrote data_corrections.csv\n";

    return 0;
}
