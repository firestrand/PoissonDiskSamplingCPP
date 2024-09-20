#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>  // For std::accumulate
#include "PoissonDiskSampling.h"

// Function to calculate Euclidean distance between two points
double calculateDistance(const std::vector<double>& a, const std::vector<double>& b) {
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
            if (calculateDistance(points[i], points[j]) < minDistance) {
                std::cerr << "Points too close: (" << i << ") and (" << j << ") are closer than " << minDistance << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Function to map a point to a bin in the spatial grid
std::vector<int> getBin(const std::vector<double>& point, const std::vector<double>& lowerBounds,
                        const std::vector<double>& upperBounds, int numBins) {
    std::vector<int> bin(point.size());
    for (size_t i = 0; i < point.size(); ++i) {
        double normalized = (point[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);  // Normalize to [0, 1]
        bin[i] = static_cast<int>(std::floor(normalized * numBins));  // Map to bin index
        if (bin[i] == numBins) bin[i] = numBins - 1;  // Ensure that the max value doesn't exceed numBins-1
    }
    return bin;
}

// Function to validate that the points are evenly distributed in space with tolerance for randomness
bool validateEvenDistribution(const std::vector<std::vector<double>>& points, const std::vector<double>& lowerBounds,
                              const std::vector<double>& upperBounds, int numBins) {
    std::unordered_map<std::vector<int>, int, VectorHash> binCount;  // To store how many points fall into each bin

    // Count how many points fall into each bin
    for (const auto& point : points) {
        std::vector<int> bin = getBin(point, lowerBounds, upperBounds, numBins);
        binCount[bin]++;
    }

    // Analyze the bin counts
    int totalBins = std::pow(numBins, points[0].size());  // Total number of bins in n-dimensional space
    std::vector<int> counts(totalBins, 0);

    // Collect the bin counts
    for (const auto& bin : binCount) {
        counts.push_back(bin.second);
    }

    // Calculate the mean and variance of the point distribution across the bins
    double mean = static_cast<double>(points.size()) / totalBins;
    double variance = 0.0;
    for (int count : counts) {
        variance += (count - mean) * (count - mean);
    }
    variance /= totalBins;

    std::cout << "Mean number of points per bin: " << mean << std::endl;
    std::cout << "Variance in bin counts: " << variance << std::endl;

    // Return true if the variance is close to the mean, indicating a Poisson-like distribution
    double expectedVariance = mean;  // For a Poisson distribution, variance should be close to mean
    double varianceRatio = variance / expectedVariance;
    std::cout << "Variance-to-Mean Ratio: " << varianceRatio << std::endl;

    // Allow some tolerance for variance ratio, e.g., accept variance within a factor of 2 of the mean
    return varianceRatio < 2.0;
}

int main() {
    // Example setup for 3D Poisson Disk Sampling
    int dimensions = 3;  // 3D space
    double minDistance = 1.0;  // Minimum distance between points
    std::vector<double> lowerBounds = {0.0, 0.0, 0.0};  // Lower bounds for the space
    std::vector<double> upperBounds = {10.0, 10.0, 10.0};  // Upper bounds for the space
    int numPoints = 100;  // Number of points to generate
    int numBins = 10;  // Number of bins along each axis

    // Create an instance of the PoissonDiskSampling class
    PoissonDiskSampling sampler(dimensions, minDistance, lowerBounds, upperBounds);

    // Generate points
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
