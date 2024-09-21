#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <boost/unordered_map.hpp>
#include "poisson_disk_sampling.h"

// Function to calculate Euclidean euclidean_distance between two points
double calculateEuclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// Function to check that no points are closer than the minimum euclidean_distance
bool validateMinDistance(const std::vector<std::vector<double>>& points, double minDistance) {
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            if (calculateEuclideanDistance(points[i], points[j]) < minDistance) {
                std::cerr << "Points too close: (" << i << ") and (" << j << ") are closer than " << minDistance << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Function to map a point to a bin in the spatial grid
std::vector<int> getBin(const std::vector<double>& point, const std::vector<double>& lowerBounds,
                        int numBins, const std::vector<double>& inverseBoundsDiff) {
    std::vector<int> bin(point.size());

    for (size_t i = 0; i < point.size(); ++i) {
        double normalized = (point[i] - lowerBounds[i]) * inverseBoundsDiff[i];
        bin[i] = static_cast<int>(std::floor(normalized * numBins));
        if (bin[i] == numBins) bin[i] = numBins - 1;
    }

    return bin;
}

// Function to validate that the points are evenly distributed in space
bool validateEvenDistribution(const std::vector<std::vector<double>>& points, const std::vector<double>& lowerBounds,
                              const std::vector<double>& upperBounds, int numBins) {
    size_t dimensions = lowerBounds.size();
    int totalBins = static_cast<int>(std::pow(numBins, dimensions));
    boost::unordered_map<std::vector<int>, int> binCount;

    std::vector<double> inverseBoundsDiff(dimensions);
    for (size_t i = 0; i < dimensions; ++i) {
        inverseBoundsDiff[i] = 1.0 / (upperBounds[i] - lowerBounds[i]);
    }

    for (const auto& point : points) {
        std::vector<int> bin = getBin(point, lowerBounds, numBins, inverseBoundsDiff);
        binCount[bin]++;
    }

    std::vector<int> counts;
    for (const auto& bin : binCount) {
        counts.push_back(bin.second);
    }

    double mean = static_cast<double>(points.size()) / totalBins;
    double variance = std::accumulate(counts.begin(), counts.end(), 0.0, [mean](double acc, int count) {
        return acc + (count - mean) * (count - mean);
    }) / totalBins;

    std::cout << "Mean: " << mean << ", Variance: " << variance << std::endl;
    double varianceRatio = variance / mean;
    std::cout << "Variance-to-Mean Ratio: " << varianceRatio << std::endl;

    const double VARIANCE_TOLERANCE = 2.0;
    return varianceRatio < VARIANCE_TOLERANCE;
}

int main() {
    int D = 4;
    int N = 100;
    int k = 30;
    std::vector<double> lowerBounds(4, 0.0);
    std::vector<double> upperBounds(4, 10.0);
    int numBins = 10;

    poisson_disk_sampling pds(D, N, k, lowerBounds, upperBounds);
    std::vector<std::vector<double>> points = pds.generatePoints();

    std::cout << "Generated " << points.size() << " points:" << std::endl;
    for (const auto& point : points) {
        std::cout << "(";
        for (size_t i = 0; i < point.size(); ++i) {
            std::cout << point[i];
            if (i < point.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ")" << std::endl;
    }

    // Extract minDistance from the poisson_disk_sampling object
    double minDistance = pds.get_r();

    // Validate the minimum euclidean_distance between points
    if (validateMinDistance(points, minDistance)) {
        std::cout << "All points are at least " << minDistance << " units apart." << std::endl;
    } else {
        std::cerr << "Validation failed: Some points are too close!" << std::endl;
    }

    // Validate even distribution across the space
    if (validateEvenDistribution(points, lowerBounds, upperBounds, numBins)) {
        std::cout << "Points are evenly distributed through the space." << std::endl;
    } else {
        std::cerr << "Validation failed: Points are not evenly distributed!" << std::endl;
    }

    return 0;
}