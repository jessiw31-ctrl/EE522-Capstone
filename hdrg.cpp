#include <algorithm>
#include <cstdlib>
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

struct Defect {
    int id;
    int shot;
    int x;       // edge index
    int t;       // round index
    int charge;  // 1 or 2
    bool final_boundary;
};

struct Cluster {
    int id;
    std::vector<int> defect_indices;
    int total_charge_mod3;
    bool neutral;
    int radius;
};

struct CorrectionChain {
    int cluster_id;
    int defect_a;
    int defect_b;
    int charge;
    std::vector<std::pair<int, int>> path;
};

struct HDRGResult {
    bool success;
    int final_radius;
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

std::map<int, std::vector<Defect>> read_syndrome_csv(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Could not open " + filename);
    }

    std::string header;
    std::getline(in, header);

    std::map<int, std::vector<Defect>> shots;
    std::string line;

    int global_id = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        auto cols = split_csv_line(line);

        int shot = std::stoi(cols[3]);
        int round = std::stoi(cols[4]);
        int edge = std::stoi(cols[5]);
        bool final_boundary = std::stoi(cols[6]) != 0;
        int defect = std::stoi(cols[9]);

        if (defect == 0) continue;

        Defect d;
        d.id = global_id++;
        d.shot = shot;
        d.x = edge;
        d.t = round;
        d.charge = defect;
        d.final_boundary = final_boundary;

        shots[shot].push_back(d);
    }

    return shots;
}

int spacetime_distance(const Defect &a, const Defect &b) {
    return std::abs(a.x - b.x) + std::abs(a.t - b.t);
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
        int radius) {

    int charge = 0;

    for (int idx : component) {
        charge += defects[idx].charge;
    }

    charge = mod_q(charge);

    Cluster c;
    c.id = cluster_id;
    c.defect_indices = component;
    c.total_charge_mod3 = charge;
    c.neutral = (charge == 0);
    c.radius = radius;

    return c;
}

std::vector<std::pair<int, int>> manhattan_path(const Defect &a, const Defect &b) {
    std::vector<std::pair<int, int>> path;

    int x = a.x;
    int t = a.t;

    path.push_back({x, t});

    while (x != b.x) {
        x += (b.x > x) ? 1 : -1;
        path.push_back({x, t});
    }

    while (t != b.t) {
        t += (b.t > t) ? 1 : -1;
        path.push_back({x, t});
    }

    return path;
}

std::vector<CorrectionChain> make_correction_chains(
        const std::vector<Cluster> &clusters,
        const std::vector<Defect> &defects) {

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
            chain.path = manhattan_path(defects[a_idx], defects[b_idx]);

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
            Cluster c = make_cluster(i, comps[i], defects, radius);
            if (!c.neutral) {
                all_neutral = false;
            }
            clusters.push_back(c);
        }

        if (all_neutral) {
            result.success = true;
            result.final_radius = radius;
            result.clusters = clusters;
            result.chains = make_correction_chains(clusters, defects);
            return result;
        }
    }

    auto comps = connected_components(defects, max_radius);
    std::vector<Cluster> clusters;

    for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
        clusters.push_back(make_cluster(i, comps[i], defects, max_radius));
    }

    result.success = false;
    result.final_radius = max_radius;
    result.clusters = clusters;
    result.chains = make_correction_chains(clusters, defects);
    return result;
}

void write_graph_dot(
        const std::string &filename,
        const std::vector<Defect> &defects,
        int radius) {

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Could not write " + filename);
    }

    out << "graph G {\n";
    out << "  layout=neato;\n";
    out << "  overlap=false;\n";
    out << "  splines=true;\n";

    for (int i = 0; i < static_cast<int>(defects.size()); ++i) {
        const Defect &d = defects[i];

        std::string color = "red";
        std::string label = std::to_string(d.charge);

        out << "  " << i
            << " [label=\"" << label << "\", "
            << "pos=\"" << d.x << "," << -d.t << "!\", "
            << "style=filled, fillcolor=" << color << ", "
            << "shape=circle];\n";
    }

    for (int i = 0; i < static_cast<int>(defects.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(defects.size()); ++j) {
            int dist = spacetime_distance(defects[i], defects[j]);
            if (dist <= radius) {
                out << "  " << i << " -- " << j
                    << " [label=\"" << dist << "\"];\n";
            }
        }
    }

    out << "}\n";
}

void write_clusters_csv(
        std::ofstream &out,
        int shot,
        const HDRGResult &result,
        const std::vector<Defect> &defects) {

    for (const Cluster &c : result.clusters) {
        for (int idx : c.defect_indices) {
            const Defect &d = defects[idx];

            out << shot << ','
                << c.id << ','
                << c.radius << ','
                << c.total_charge_mod3 << ','
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
        int shot,
        const HDRGResult &result) {

    for (const CorrectionChain &chain : result.chains) {
        for (int step = 0; step < static_cast<int>(chain.path.size()); ++step) {
            out << shot << ','
                << chain.cluster_id << ','
                << chain.defect_a << ','
                << chain.defect_b << ','
                << chain.charge << ','
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

    auto shots = read_syndrome_csv(syndrome_file);

    std::ofstream summary("hdrg_summary.csv");
    std::ofstream clusters_out("clusters.csv");
    std::ofstream chains_out("correction_chains.csv");

    summary << "shot,num_defects,success,final_radius,num_clusters,num_chains\n";

    clusters_out
        << "shot,cluster_id,radius,total_charge_mod3,neutral,"
        << "defect_id,x,t,charge,final_boundary\n";

    chains_out
        << "shot,cluster_id,defect_a,defect_b,charge,path_step,x,t\n";

    int total = 0;
    int failures = 0;
    int graph_count = 0;

    for (const auto &[shot, defects] : shots) {
        HDRGResult result = hdrg_decode_one_shot(defects, distance, rounds);

        total += 1;
        if (!result.success) failures += 1;

        summary << shot << ','
                << defects.size() << ','
                << result.success << ','
                << result.final_radius << ','
                << result.clusters.size() << ','
                << result.chains.size() << '\n';

        write_clusters_csv(clusters_out, shot, result, defects);
        write_correction_chains_csv(chains_out, shot, result);

        if (graph_count < graphs_to_save) {
            std::string graph_name =
                    "graph_shot_" + std::to_string(shot) + ".dot";

            write_graph_dot(graph_name, defects, result.final_radius);
            graph_count += 1;
        }
    }

    double failure_rate =
            total > 0 ? static_cast<double>(failures) / total : 0.0;

    std::cout << "HDRG decode complete\n";
    std::cout << "decoded shots = " << total << "\n";
    std::cout << "failures = " << failures << "\n";
    std::cout << "failure rate = " << failure_rate << "\n";
    std::cout << "wrote hdrg_summary.csv\n";
    std::cout << "wrote clusters.csv\n";
    std::cout << "wrote correction_chains.csv\n";
    std::cout << "wrote first " << graph_count << " graph .dot files\n";

    return 0;
}
