#ifndef POISSON_DISK_SAMPLING_H
#define POISSON_DISK_SAMPLING_H

#include "sampling_interface.h"
#include <vector>
#include <random>
#include <boost/multi_array.hpp>

class poisson_disk_sampling  : public SamplingInterface{
public:
    poisson_disk_sampling(int dims, int num_points, int maxAttempts, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);
    std::vector<std::vector<double>> generatePoints();
    [[nodiscard]] double get_r() const { return r; }


private:
    int dimensions;
    int num_points;
    int maxAttempts;
    std::vector<double> upper_bounds;
    std::vector<double> lower_bounds;
    double r;
    boost::multi_array<std::vector<std::vector<double>>, 2> grid;
    std::vector<std::vector<double>> active_list;
    std::vector<std::vector<double>> sample_points;
    std::mt19937 gen;

    double calculate_volume();
    [[nodiscard]] double calculate_r_from_point_count(double volume) const;
    std::vector<double> random_point_in_bounds();
    std::vector<double> generate_random_point_around(const std::vector<double>& p);
    bool is_within_bounds(const std::vector<double>& point);
    bool is_valid_point(const std::vector<double>& point);
    double euclidean_distance(const std::vector<double>& p1, const std::vector<double>& p2) const;

    double cosine_similarity(const std::vector<double> &p1, const std::vector<double> &p2);
};

#endif // POISSON_DISK_SAMPLING_H