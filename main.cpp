#include <iostream>
#include <vector>
#include <cmath>
#include <boost/geometry/algorithms/distance.hpp>
#include "PoissonDiskSampling.h"

// Function to calculate Euclidean distance between two points
double calculateEuclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// Function to check that no points are closer than the minimum distance
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

    // Map each dimension to a bin
    for (size_t i = 0; i < point.size(); ++i) {
        double normalized = (point[i] - lowerBounds[i]) * inverseBoundsDiff[i];  // Normalize to [0, 1]
        bin[i] = static_cast<int>(std::floor(normalized * numBins));  // Map to bin index
        if (bin[i] == numBins) bin[i] = numBins - 1;  // Ensure max value doesn't exceed numBins-1
    }

    return bin;
}


// Function to validate that the points are evenly distributed in space
bool validateEvenDistribution(const std::vector<std::vector<double>>& points, const std::vector<double>& lowerBounds,
                              const std::vector<double>& upperBounds, int numBins) {
    size_t dimensions = lowerBounds.size();
    int totalBins = static_cast<int>(std::pow(numBins, dimensions));  // Total bins in n-dimensional space
    boost::unordered_map<std::vector<int>, int> binCount;

    // Precompute inverse of bounds difference for each dimension
    std::vector<double> inverseBoundsDiff(dimensions);
    for (size_t i = 0; i < dimensions; ++i) {
        inverseBoundsDiff[i] = 1.0 / (upperBounds[i] - lowerBounds[i]);
    }

    // Count how many points fall into each bin
    for (const auto& point : points) {
        std::vector<int> bin = getBin(point, lowerBounds, numBins, inverseBoundsDiff);
        binCount[bin]++;
    }

    std::vector<int> counts(totalBins, 0);
    for (const auto& bin : binCount) {
        counts.push_back(bin.second);
    }

    // Calculate mean and variance of point distribution
    double mean = static_cast<double>(points.size()) / totalBins;
    double variance = std::accumulate(counts.begin(), counts.end(), 0.0, [mean](double acc, int count) {
        return acc + (count - mean) * (count - mean);
    }) / totalBins;

    std::cout << "Mean: " << mean << ", Variance: " << variance << std::endl;
    double varianceRatio = variance / mean;
    std::cout << "Variance-to-Mean Ratio: " << varianceRatio << std::endl;

    // Acceptable variance tolerance
    const double VARIANCE_TOLERANCE = 2.0;
    return varianceRatio < VARIANCE_TOLERANCE;
}


int main() {
    // Example setup for 3D Poisson Disk Sampling
    int dimensions = 4;  // 3D space
    double minDistance = 1.0;  // Minimum distance between points
    std::vector<double> lowerBounds = {0.0, 0.0, 0.0, 0.0};  // Lower bounds for the space
    std::vector<double> upperBounds = {10.0, 10.0, 10.0, 10.0};  // Upper bounds for the space
    int numPoints = 100;  // Number of points to generate
    int numBins = 10;  // Number of bins along each axis

    PoissonDiskSampling sampler(dimensions, minDistance, lowerBounds, upperBounds);
    std::vector<std::vector<double>> points = sampler.generatePoints(numPoints);

    // Print the generated points
    std::cout << "Generated Points:" << std::endl;
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

    // Validate the minimum distance between points
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
