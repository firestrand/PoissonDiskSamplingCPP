#include "PoissonDiskSampling.h"
#include <boost/geometry.hpp>
#include <random>
#include <algorithm>
#include <cmath>


PoissonDiskSampling::PoissonDiskSampling(int dims, double minDist, const std::vector<double>& lower, const std::vector<double>& upper)
        : dimensions(dims), minDistance(minDist), lowerBounds(lower), upperBounds(upper),
          gen(std::random_device{}()), uniform(0.0, 1.0)
{
    cellSize = minDistance / std::sqrt(dimensions);
}

double PoissonDiskSampling::distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = std::inner_product(a.begin(), a.end(), b.begin(), 0.0, std::plus<>(),
                                    [](double ai, double bi) { return (ai - bi) * (ai - bi); });
    return std::sqrt(sum);
}

std::vector<double> PoissonDiskSampling::generateRandomPoint() {
    std::vector<double> point(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        point[i] = lowerBounds[i] + uniform(gen) * (upperBounds[i] - lowerBounds[i]);
    }
    return point;
}

std::vector<int> PoissonDiskSampling::gridCoords(const std::vector<double>& point) {
    std::vector<int> coords(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        coords[i] = static_cast<int>((point[i] - lowerBounds[i]) / cellSize);
    }
    return coords;
}

bool PoissonDiskSampling::isValidPoint(const std::vector<double>& point) {
    std::vector<int> coords = gridCoords(point);
    std::vector<std::vector<int>> neighbors;

    // Recursively generate neighboring coordinates
    generateNeighbors(coords, 0, neighbors);

    // Check if any of the neighboring points are too close
    for (const auto& neighbor : neighbors) {
        if (grid.find(neighbor) != grid.end()) {
            for (const auto& p : grid[neighbor]) {
                if (distance(p, point) < minDistance) {
                    return false;
                }
            }
        }
    }
    return true;
}

// Helper function to recursively generate neighbors in d-dimensions
void PoissonDiskSampling::generateNeighbors(std::vector<int>& coords, int dim, std::vector<std::vector<int>>& neighbors) {
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

std::vector<std::vector<double>> PoissonDiskSampling::generatePoints(int numPoints, int maxAttempts) {
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