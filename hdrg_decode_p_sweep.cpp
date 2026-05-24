#include <algorithm>
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
    int id;
    std::string p;
    int shot;
    int x;
    int t;
    int charge;
    bool final_boundary;
};

struct Cluster {
    int id;
    std::vector<int> defect_indices;
    int total_charge_mod3;
    bool touches_left_boundary;
    bool touches_right_boundary;
    bool neutral;
    int radius;
};

struct CorrectionChain {
    int cluster_id;
    int defect_a;
    int defect_b;
    int charge;
    bool to_boundary;
    std::string boundary_name;
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

    while (std::getline(ss, item, ',')) {
        out.push_back(item);
    }

    return out;
}

std::map<ShotKey, std::vector<Defect>>
read_syndrome_csv(const std::string &filename) {
    std::ifstream in(filename);

    if (!in) {
        throw std::runtime_error("Could not open " + filename);
    }

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

        if (defect == 0) continue;

        Defect d;
        d.id = global_id++;
        d.p = p;
        d.shot = shot;
        d.x = edge;
        d.t = round;
        d.charge = defect;
        d.final_boundary = final_boundary;

        ShotKey key{p, shot};
        shots[key].push_back(d);
    }

    return shots;
}

int spacetime_distance(const Defect &a, const Defect &b) {
    return std::abs(a.x - b.x) + std::abs(a.t - b.t);
}

int distance_to_left_boundary(const Defect &d) {
    return d.x + 1;
}

int distance_to_right_boundary(const Defect &d, int distance) {
    return (distance - 1) - d.x;
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
        int distance) {

    int charge = 0;
    bool touches_left = false;
    bool touches_right = false;

    for (int idx : component) {
        const Defect &d = defects[idx];

        charge += d.charge;

        if (distance_to_left_boundary(d) <= radius) {
            touches_left = true;
        }

        if (distance_to_right_boundary(d, distance) <= radius) {
            touches_right = true;
        }

        if (d.final_boundary) {
            touches_left = true;
            touches_right = true;
        }
    }

    charge = mod_q(charge);

    Cluster c;
    c.id = cluster_id;
    c.defect_indices = component;
    c.total_charge_mod3 = charge;
    c.touches_left_boundary = touches_left;
    c.touches_right_boundary = touches_right;
    c.neutral = (charge == 0) || touches_left || touches_right;
    c.radius = radius;

    return c;
}

std::vector<std::pair<int, int>> manhattan_path_between_points(
        int x1, int t1,
        int x2, int t2) {

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

std::vector<std::pair<int, int>> manhattan_path(
        const Defect &a,
        const Defect &b) {

    return manhattan_path_between_points(a.x, a.t, b.x, b.t);
}

std::vector<std::pair<int, int>> path_to_boundary(
        const Defect &a,
        const std::string &boundary_name,
        int distance) {

    int boundary_x = (boundary_name == "L") ? -1 : distance - 1;
    return manhattan_path_between_points(a.x, a.t, boundary_x, a.t);
}

std::vector<CorrectionChain> make_correction_chains(
        const std::vector<Cluster> &clusters,
        const std::vector<Defect> &defects,
        int distance) {

    std::vector<CorrectionChain> chains;

    for (const Cluster &cluster : clusters) {
        if (!cluster.neutral) continue;

        std::vector<int> charge1;
        std::vector<int> charge2;

        for (int idx : cluster.defect_indices) {
            if (defects[idx].charge == 1) {
                charge1.push_back(idx);
            } else if (defects[idx].charge == 2) {
                charge2.push_back(idx);
            }
        }

        while (!charge1.empty() && !charge2.empty()) {
            int a_idx = charge1.back();
            charge1.pop_back();

            int best_j = 0;
            int best_dist = 1000000000;

            for (int j = 0; j < static_cast<int>(charge2.size()); ++j) {
                int dist = spacetime_distance(defects[a_idx], defects[charge2[j]]);
                if (dist < best_dist) {
                    best_dist = dist;
                    best_j = j;
                }
            }

            int b_idx = charge2[best_j];
            charge2.erase(charge2.begin() + best_j);

            CorrectionChain chain;
            chain.cluster_id = cluster.id;
            chain.defect_a = defects[a_idx].id;
            chain.defect_b = defects[b_idx].id;
            chain.charge = defects[a_idx].charge;
            chain.to_boundary = false;
            chain.boundary_name = "";
            chain.path = manhattan_path(defects[a_idx], defects[b_idx]);

            chains.push_back(chain);
        }

        std::vector<int> leftovers;
        leftovers.insert(leftovers.end(), charge1.begin(), charge1.end());
        leftovers.insert(leftovers.end(), charge2.begin(), charge2.end());

        for (int idx : leftovers) {
            const Defect &d = defects[idx];

            if (!cluster.touches_left_boundary && !cluster.touches_right_boundary) {
                continue;
            }

            std::string boundary_name;

            if (cluster.touches_left_boundary && cluster.touches_right_boundary) {
                int dl = distance_to_left_boundary(d);
                int dr = distance_to_right_boundary(d, distance);
                boundary_name = (dl <= dr) ? "L" : "R";
            } else if (cluster.touches_left_boundary) {
                boundary_name = "L";
            } else {
                boundary_name = "R";
            }

            CorrectionChain chain;
            chain.cluster_id = cluster.id;
            chain.defect_a = d.id;
            chain.defect_b = -1;
            chain.charge = d.charge;
            chain.to_boundary = true;
            chain.boundary_name = boundary_name;
            chain.path = path_to_boundary(d, boundary_name, distance);

            chains.push_back(chain);
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

    int max_radius = distance + rounds + 2;

    for (int radius = 1; radius <= max_radius; ++radius) {
        auto comps = connected_components(defects, radius);

        std::vector<Cluster> clusters;
        bool all_neutral = true;

        for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
            Cluster c = make_cluster(i, comps[i], defects, radius, distance);

            if (!c.neutral) {
                all_neutral = false;
            }

            clusters.push_back(c);
        }

        if (all_neutral) {
            result.success = true;
            result.final_radius = radius;
            result.clusters = clusters;
            result.chains = make_correction_chains(clusters, defects, distance);
            return result;
        }
    }

    auto comps = connected_components(defects, max_radius);

    std::vector<Cluster> clusters;
    bool all_neutral = true;

    for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
        Cluster c = make_cluster(i, comps[i], defects, max_radius, distance);

        if (!c.neutral) {
            all_neutral = false;
        }

        clusters.push_back(c);
    }

    result.success = all_neutral;
    result.final_radius = max_radius;
    result.clusters = clusters;
    result.chains = make_correction_chains(clusters, defects, distance);

    return result;
}

void write_graph_dot(
        const std::string &filename,
        const std::vector<Defect> &defects,
        const HDRGResult &result,
        int distance,
        int rounds) {

    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error("Could not write " + filename);
    }

    out << "graph G {\n";
    out << "  layout=neato;\n";
    out << "  overlap=false;\n";
    out << "  splines=true;\n";
    out << "  node [fixedsize=true, width=0.38, height=0.38];\n";

    for (int t = 0; t <= rounds; ++t) {
        for (int x = 0; x <= distance - 2; ++x) {
            out << "  v_" << x << "_" << t
                << " [label=\"\", pos=\"" << x << "," << -t << "!\", "
                << "style=filled, fillcolor=lightgray, shape=circle];\n";
        }
    }

    for (int t = 0; t <= rounds; ++t) {
        out << "  boundaryL_" << t
            << " [label=\"\", pos=\"-1," << -t << "!\", "
            << "style=filled, fillcolor=black, shape=square, "
            << "width=0.20, height=0.20];\n";

        out << "  boundaryR_" << t
            << " [label=\"\", pos=\"" << distance - 1 << "," << -t << "!\", "
            << "style=filled, fillcolor=black, shape=square, "
            << "width=0.20, height=0.20];\n";
    }

    for (int t = 0; t <= rounds; ++t) {
        for (int x = 0; x < distance - 2; ++x) {
            out << "  v_" << x << "_" << t
                << " -- v_" << x + 1 << "_" << t
                << " [color=gray, penwidth=1];\n";
        }
    }

    for (int t = 0; t < rounds; ++t) {
        for (int x = 0; x <= distance - 2; ++x) {
            out << "  v_" << x << "_" << t
                << " -- v_" << x << "_" << t + 1
                << " [color=gray, penwidth=1];\n";
        }
    }

    for (int i = 0; i < static_cast<int>(defects.size()); ++i) {
        const Defect &d = defects[i];
        std::string color = (d.charge == 1) ? "red" : "orange";

        out << "  d" << i
            << " [label=\"" << d.charge << "\", "
            << "pos=\"" << d.x << "," << -d.t << "!\", "
            << "style=filled, fillcolor=" << color << ", "
            << "shape=circle, width=0.60, height=0.60];\n";
    }

    for (const Cluster &c : result.clusters) {
        for (int a = 0; a < static_cast<int>(c.defect_indices.size()); ++a) {
            for (int b = a + 1; b < static_cast<int>(c.defect_indices.size()); ++b) {
                int i = c.defect_indices[a];
                int j = c.defect_indices[b];
                int dist = spacetime_distance(defects[i], defects[j]);

                if (dist <= c.radius) {
                    std::string edge_color = c.neutral ? "green" : "red";

                    out << "  d" << i << " -- d" << j
                        << " [color=" << edge_color
                        << ", penwidth=3, label=\"" << dist << "\"];\n";
                }
            }
        }

        for (int idx : c.defect_indices) {
            const Defect &d = defects[idx];

            if (c.touches_left_boundary) {
                out << "  d" << idx
                    << " -- boundaryL_" << d.t
                    << " [color=black, style=dashed, penwidth=2];\n";
            }

            if (c.touches_right_boundary) {
                out << "  d" << idx
                    << " -- boundaryR_" << d.t
                    << " [color=black, style=dashed, penwidth=2];\n";
            }
        }
    }

    int chain_counter = 0;

    for (const CorrectionChain &chain : result.chains) {
        for (int step = 0; step < static_cast<int>(chain.path.size()); ++step) {
            auto p = chain.path[step];

            std::string n =
                    "p_" + std::to_string(chain_counter)
                    + "_" + std::to_string(step);

            out << "  " << n
                << " [label=\"\", pos=\"" << p.first << "," << -p.second << "!\", "
                << "width=0.13, height=0.13, "
                << "style=filled, fillcolor=blue, shape=circle];\n";

            if (step > 0) {
                std::string prev =
                        "p_" + std::to_string(chain_counter)
                        + "_" + std::to_string(step - 1);

                out << "  " << prev << " -- " << n
                    << " [color=blue, penwidth=3];\n";
            }
        }

        chain_counter += 1;
    }

    out << "}\n";
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
                << step << ','
                << chain.path[step].first << ','
                << chain.path[step].second << '\n';
        }
    }
}

int main(int argc, char **argv) {
    std::string syndrome_file = "syndrome.csv";
    int distance = 3;
    int rounds = 5;
    int graphs_to_save = 3;

    if (argc >= 2) syndrome_file = argv[1];
    if (argc >= 3) distance = std::stoi(argv[2]);
    if (argc >= 4) rounds = std::stoi(argv[3]);
    if (argc >= 5) graphs_to_save = std::stoi(argv[4]);

    auto shots = read_syndrome_csv(syndrome_file);

    std::ofstream summary("hdrg_summary.csv");
    std::ofstream clusters_out("clusters.csv");
    std::ofstream chains_out("correction_chains.csv");

    summary
        << "p,shot,num_defects,success,final_radius,num_clusters,num_chains\n";

    clusters_out
        << "p,shot,cluster_id,radius,total_charge_mod3,"
        << "touches_left_boundary,touches_right_boundary,neutral,"
        << "defect_id,x,t,charge,final_boundary\n";

    chains_out
        << "p,shot,cluster_id,defect_a,defect_b,charge,"
        << "to_boundary,boundary_name,path_step,x,t\n";

    int failures = 0;
    int graph_count = 0;

    for (const auto &entry : shots) {
        const ShotKey &key = entry.first;
        const std::vector<Defect> &defects = entry.second;

        HDRGResult result = hdrg_decode_one_shot(defects, distance, rounds);

        if (!result.success) {
            failures += 1;
        }

        summary << key.p << ','
                << key.shot << ','
                << defects.size() << ','
                << result.success << ','
                << result.final_radius << ','
                << result.clusters.size() << ','
                << result.chains.size() << '\n';

        write_clusters_csv(clusters_out, key, result, defects);
        write_correction_chains_csv(chains_out, key, result);

        if (graph_count < graphs_to_save && !defects.empty()) {
            std::string graph_name =
                    "graph_p_" + key.p
                    + "_shot_" + std::to_string(key.shot)
                    + ".dot";

            std::replace(graph_name.begin(), graph_name.end(), '.', '_');

            graph_name += ".dot";

            write_graph_dot(graph_name, defects, result, distance, rounds);
            graph_count += 1;
        }
    }

    double failure_rate =
            !shots.empty()
            ? static_cast<double>(failures) / static_cast<double>(shots.size())
            : 0.0;

    std::cout << "HDRG decode complete\n";
    std::cout << "decoded non-empty syndrome shots = " << shots.size() << "\n";
    std::cout << "cluster-neutralization failures = " << failures << "\n";
    std::cout << "cluster-neutralization failure rate = " << failure_rate << "\n";
    std::cout << "wrote hdrg_summary.csv\n";
    std::cout << "wrote clusters.csv\n";
    std::cout << "wrote correction_chains.csv\n";
    std::cout << "wrote first " << graph_count << " graph .dot files\n";

    return 0;
}
