#include "poisson_disk_sampling.h"

#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/container/small_vector.hpp>
#include <vector>
#include <iostream>

poisson_disk_sampling::poisson_disk_sampling(int dims, int num_points, int maxAttempts, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds)
        : dimensions(dims), num_points(num_points), maxAttempts(maxAttempts), lower_bounds(lower_bounds), upper_bounds(upper_bounds), gen(std::random_device{}()) {
    double volume = calculate_volume();
    r = calculate_r_from_point_count(volume);
}

double poisson_disk_sampling::calculate_volume() {
    double volume = 1.0;
    for (int d = 0; d < dimensions; ++d) {
        volume *= (upper_bounds[d] - lower_bounds[d]);
    }
    return volume;
}

double poisson_disk_sampling::calculate_r_from_point_count(double volume) const {
    static const double unit_sphere_volume = std::pow(M_PI, dimensions / 2.0) / boost::math::tgamma(dimensions / 2.0 + 1);
    double sphere_volume_per_point = volume / num_points;
    return std::pow(sphere_volume_per_point / unit_sphere_volume, 1.0 / dimensions);
}

std::vector<std::vector<double>> poisson_disk_sampling::generatePoints() {
    std::vector<double> p0 = random_point_in_bounds();  // Initial point in bounds
    sample_points.reserve(num_points);
    active_list.reserve(num_points);

    sample_points.push_back(p0);
    active_list.push_back(p0);

    boost::random::uniform_int_distribution<> dist(0, active_list.size() - 1);

    while (!active_list.empty() && sample_points.size() < num_points) {
        int idx = dist(gen);

        // Safely reference the point
        auto& p = active_list[idx];
        bool point_found = false;

        // Try generating up to maxAttempts random points around the selected point
        for (int i = 0; i < maxAttempts; ++i) {
            auto new_point = generate_random_point_around(p);
            if (is_within_bounds(new_point) && is_valid_point(new_point)) {
                sample_points.push_back(std::move(new_point));
                active_list.push_back(sample_points.back());
                point_found = true;
                break;
            }
        }

        if (!point_found) {
            std::swap(active_list[idx], active_list.back());
            active_list.pop_back();
        }
    }

    return sample_points;
}

std::vector<double> poisson_disk_sampling::random_point_in_bounds() {
    std::vector<double> point(dimensions);  // Vector is already sized dimensions
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int d = 0; d < dimensions; ++d) {
        point[d] = lower_bounds[d] + dis(gen) * (upper_bounds[d] - lower_bounds[d]);
    }
    return point;
}

std::vector<double> poisson_disk_sampling::generate_random_point_around(const std::vector<double>& p) {
    std::vector<double> direction(dimensions);
    std::vector<double> new_point(dimensions);
    std::normal_distribution<> normal(0.0, 1.0);
    std::uniform_real_distribution<> dis_radius(r, 2 * r);
    double squared_norm = 0.0;

    for (int i = 0; i < dimensions; ++i) {
        direction[i] = normal(gen);
        squared_norm += direction[i] * direction[i];
    }

    double norm = std::sqrt(squared_norm);
    for (int i = 0; i < dimensions; ++i) {
        direction[i] /= norm;
    }

    double radius = dis_radius(gen);
    for (int i = 0; i < dimensions; ++i) {
        new_point[i] = p[i] + radius * direction[i];
    }

    return new_point;
}

bool poisson_disk_sampling::is_within_bounds(const std::vector<double>& point) {
    for (int d = 0; d < dimensions; ++d) {
        if (point[d] < lower_bounds[d] || point[d] > upper_bounds[d]) {
            return false;
        }
    }
    return true;
}

bool poisson_disk_sampling::is_valid_point(const std::vector<double>& point) {
    for (const auto& existing_point : sample_points) {
        double dist = euclidean_distance(point, existing_point);
        if (dist < r) {
            std::cout << "Points too close: " << dist << " (threshold: " << r << ")" << std::endl;
            return false;
        }
    }
    return true;
}

double poisson_disk_sampling::euclidean_distance(const std::vector<double>& p1, const std::vector<double>& p2) const {
    double sum = 0;
    for (int d = 0; d < dimensions; ++d) {
        double diff = p1[d] - p2[d];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

double poisson_disk_sampling::cosine_similarity(const std::vector<double>& p1, const std::vector<double>& p2) {
    double dot_product = 0.0, norm1 = 0.0, norm2 = 0.0;
    for (int d = 0; d < dimensions; ++d) {
        dot_product += p1[d] * p2[d];
        norm1 += p1[d] * p1[d];
        norm2 += p2[d] * p2[d];
    }
    return 1.0 - (dot_product / (std::sqrt(norm1) * std::sqrt(norm2)));
}
