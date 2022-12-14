#include <iostream>
#include <filesystem>
#include <unordered_map>
#include <mpfr.h>
#include <bvh/triangle.hpp>
#include <bvh/bvh.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/primitive_intersectors.hpp>
#include <bvh/single_ray_traverser.hpp>
#include "happly.h"
#include "mark.hpp"

const mpfr_prec_t mantissa_width = 7;
const mpfr_exp_t exponent_width = 8;

template <typename Bvh, typename Traverser, typename PrimitiveIntersector>
void traverse(Bvh &bvh, Traverser &traverser, PrimitiveIntersector &primitive_intersector,
              bool output_graph, std::unordered_map<size_t, size_t> &parent) {
    uintmax_t traversal_steps_sum = 0;
    uintmax_t intersections_sum = 0;
    uintmax_t high_precision_sum = 0;
    uintmax_t low_precision_sum = 0;
    uintmax_t count = 0;

    srand(0);
    std::ifstream ray_queries_file("ray_queries.bin", std::ios::in | std::ios::binary);
    float r[7];
    for (int i = 0; ray_queries_file.read(reinterpret_cast<char*>(&r), 7 * sizeof(float)); i++) {
        if (rand() % 1000 != 0) continue;

        bvh::Ray<float> ray(
                bvh::Vector3<float>(r[0], r[1], r[2]),
                bvh::Vector3<float>(r[3], r[4], r[5]),
                0.f,
                r[6]
        );

        typename Traverser::Statistics statistics;
        traverser.traverse(ray, primitive_intersector, statistics);

        if (output_graph) {
            std::filesystem::create_directory("graphs");
            std::string filepath = "graphs/bvh_" + std::to_string(i) + ".dot";
            std::ofstream dot_file(filepath);
            dot_file << "digraph bvh {";
            dot_file << "\n    node [shape=point]";
            dot_file << "\n    edge [arrowhead=none]";

            for (auto &n : statistics.traversed) {
                dot_file << "\n    " << parent[n] << " -> " << n;
                if (bvh.nodes[n].high_precision)  dot_file << " [color=red]";
            }

            dot_file << "\n}";
            dot_file.close();
        }

        traversal_steps_sum += statistics.traversal_steps;
        intersections_sum += statistics.intersections;
        high_precision_sum += statistics.high_precision_count;
        low_precision_sum += statistics.low_precision_count;
        count += 1;
    }

    std::cout << double(traversal_steps_sum) / double(count) << " "
              << double(intersections_sum) / double(count) << " "
              << double(high_precision_sum) / double(count) << " "
              << double(low_precision_sum) / double(count) << std::endl;
}

int main() {
    happly::PLYData ply_data("model.ply");
    std::vector<std::array<double, 3>> v_pos = ply_data.getVertexPositions();
    std::vector<std::vector<size_t>> f_idx = ply_data.getFaceIndices<size_t>();

    std::vector<bvh::Triangle<float>> triangles;
    for (auto &face : f_idx) {
        triangles.emplace_back(bvh::Vector3<float>(v_pos[face[0]][0], v_pos[face[0]][1], v_pos[face[0]][2]),
                               bvh::Vector3<float>(v_pos[face[1]][0], v_pos[face[1]][1], v_pos[face[1]][2]),
                               bvh::Vector3<float>(v_pos[face[2]][0], v_pos[face[2]][1], v_pos[face[2]][2]));
    }

    auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(triangles.data(), triangles.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), triangles.size());
    std::cout << "global bounding box: ("
              << global_bbox.min[0] << ", " << global_bbox.min[1] << ", " << global_bbox.min[2] << "), ("
              << global_bbox.max[0] << ", " << global_bbox.max[1] << ", " << global_bbox.max[2] << ")" << std::endl;

    bvh::Bvh<float> bvh;
    bvh::SweepSahBuilder<bvh::Bvh<float>> builder(bvh);
    builder.build(global_bbox, bboxes.get(), centers.get(), triangles.size());

    std::unordered_map<size_t, size_t> parent;
    std::queue<size_t> queue;
    queue.push(0);
    while (!queue.empty()) {
        size_t curr = queue.front();
        queue.pop();
        if (!bvh.nodes[curr].is_leaf()) {
            size_t left_idx = bvh.nodes[curr].first_child_or_primitive;
            size_t right_idx = left_idx + 1;
            parent[left_idx] = curr;
            parent[right_idx] = curr;
            queue.push(left_idx);
            queue.push(right_idx);
        }
    }

    std::cout << "All high: ";
    bvh::SingleRayTraverser<bvh::Bvh<float>> traverser_high(bvh);
    bvh::ClosestPrimitiveIntersector<bvh::Bvh<float>, bvh::Triangle<float>> primitive_intersector(bvh, triangles.data());
    traverse(bvh, traverser_high, primitive_intersector, false, parent);

    std::cout << "All low: ";
    for (size_t i = 0; i < bvh.node_count; i++) bvh.nodes[i].high_precision = false;
    using MPNodeIntersector = bvh::MPNodeIntersector<bvh::Bvh<float>, mantissa_width, exponent_width>;
    bvh::SingleRayTraverser<bvh::Bvh<float>, 64, MPNodeIntersector> traverser_mixed(bvh);
    traverse(bvh, traverser_mixed, primitive_intersector, false, parent);

    for (int i = -200; i <= 55; i++) {
        float t_trav_high = 0.5;
        float t_trav_low = float(i) / 100;

        HighPrecisionMarker high_precision_marker(mantissa_width, exponent_width);
        high_precision_marker.mark(bvh, t_trav_high, t_trav_low);

        std::cout << t_trav_high << " " << t_trav_low << " ";
        traverse(bvh, traverser_mixed, primitive_intersector, false, parent);
    }
}
