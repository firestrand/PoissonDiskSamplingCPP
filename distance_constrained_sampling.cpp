#include "distance_constrained_sampling.h"
#include "distance_metrics.h"
#include <boost/geometry.hpp>
#include <random>
#include <algorithm>
#include <cmath>


DistanceConstrainedSampling::DistanceConstrainedSampling(int dims, int numPoints, int maxAttempts, const std::vector<double>& lower, const std::vector<double>& upper)
        : dimensions(dims), numPoints(numPoints), maxAttempts(maxAttempts),lowerBounds(lower), upperBounds(upper),
          gen(std::random_device{}()), uniform(0.0, 1.0)
{
    cellSize = calculateCellSize();
}

double DistanceConstrainedSampling::calculateCellSize() {
    double volume = 1.0;
    for (int i = 0; i < dimensions; ++i) {
        volume *= (upperBounds[i] - lowerBounds[i]);
    }

    double volumePerPoint = volume / numPoints;
    double radius = std::pow((volumePerPoint * std::tgamma(dimensions / 2.0 + 1)) / std::pow(M_PI, dimensions / 2.0), 1.0 / dimensions);

    return 2 * radius / std::sqrt(dimensions);
}

std::vector<double> DistanceConstrainedSampling::generateRandomPoint() {
    std::vector<double> point(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        point[i] = lowerBounds[i] + uniform(gen) * (upperBounds[i] - lowerBounds[i]);
    }
    return point;
}

std::vector<int> DistanceConstrainedSampling::gridCoords(const std::vector<double>& point) {
    std::vector<int> coords(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        coords[i] = static_cast<int>((point[i] - lowerBounds[i]) / cellSize);
    }
    return coords;
}

bool DistanceConstrainedSampling::isValidPoint(const std::vector<double>& point) {
    std::vector<int> coords = gridCoords(point);
    std::vector<std::vector<int>> neighbors;

    // Recursively generatePoints neighboring coordinates
    generateNeighbors(coords, 0, neighbors);

    // Check if any of the neighboring points are too close
    for (const auto& neighbor : neighbors) {
        if (grid.find(neighbor) != grid.end()) {
            for (const auto& p : grid[neighbor]) {
                if (DistanceMetrics::euclidean(p, point) < minDistance) {
                    return false;
                }
            }
        }
    }
    return true;
}

// Helper function to recursively generatePoints neighbors in d-dimensions
void DistanceConstrainedSampling::generateNeighbors(std::vector<int>& coords, int dim, std::vector<std::vector<int>>& neighbors) {
    static const int neighborRange = 2; // +/- 2 cells in each direction

    // If we have reached the last dimension, add this neighbor configuration
    if (dim == dimensions) {
        neighbors.push_back(coords);
        return;
    }

    // Recursively visit all neighbors in this dimension
    for (int offset = -neighborRange; offset <= neighborRange; ++offset) {
        int originalCoord = coords[dim];
        coords[dim] += offset;  // Modify the coordinate in the current dimension
        generateNeighbors(coords, dim + 1, neighbors);  // Move to the next dimension
        coords[dim] = originalCoord;  // Restore the original value after recursion
    }
}

std::vector<std::vector<double>> DistanceConstrainedSampling::generatePoints() {
    while (points.size() < numPoints) {
        std::vector<double> candidate = generateRandomPoint();
        for (int attempts = 0; attempts < maxAttempts; ++attempts) {
            if (isValidPoint(candidate)) {
                points.push_back(candidate);
                grid[gridCoords(candidate)].push_back(candidate);
                break;
            }
            candidate = generateRandomPoint();
        }
    }
    return points;
}