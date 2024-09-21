#ifndef DISTANCE_METRICS_H
#define DISTANCE_METRICS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

class DistanceMetrics {
public:
    // Static method to calculate Euclidean distance
    static double euclidean(const std::vector<double>& a, const std::vector<double>& b) {
        check_dimensions(a, b);
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }

    // Static method to calculate Manhattan distance
    static double manhattan(const std::vector<double>& a, const std::vector<double>& b) {
        check_dimensions(a, b);
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += std::abs(a[i] - b[i]);
        }
        return sum;
    }

    // Static method to calculate Chebyshev distance
    static double chebyshev(const std::vector<double>& a, const std::vector<double>& b) {
        check_dimensions(a, b);
        double max_diff = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            max_diff = std::max(max_diff, std::abs(a[i] - b[i]));
        }
        return max_diff;
    }

    // Static method to calculate Minkowski distance
    static double minkowski(const std::vector<double>& a, const std::vector<double>& b, double p) {
        check_dimensions(a, b);
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += std::pow(std::abs(a[i] - b[i]), p);
        }
        return std::pow(sum, 1.0 / p);
    }

private:
    // Utility method to ensure both vectors have the same dimensions
    static void check_dimensions(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same dimensions.");
        }
    }
};

#endif // DISTANCE_METRICS_H

