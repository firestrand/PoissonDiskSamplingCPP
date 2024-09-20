#include "PoissonDiskSampling.h"
#include <cmath>
#include <algorithm>

std::size_t VectorHash::operator()(const std::vector<int>& v) const {
    std::hash<int> hasher;
    std::size_t seed = 0;
    for (int i : v) {
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

PoissonDiskSampling::PoissonDiskSampling(int dims, double minDist, const std::vector<double>& lower, const std::vector<double>& upper)
        : dimensions(dims), minDistance(minDist), lowerBounds(lower), upperBounds(upper),
          gen(std::random_device{}()), dis(0.0, 1.0)
{
    cellSize = minDistance / std::sqrt(dimensions);
}

double PoissonDiskSampling::distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0;
    for (int i = 0; i < dimensions; ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

std::vector<double> PoissonDiskSampling::generateRandomPoint() {
    std::vector<double> point(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        point[i] = lowerBounds[i] + dis(gen) * (upperBounds[i] - lowerBounds[i]);
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
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            for (int k = -2; k <= 2; ++k) {
                std::vector<int> neighbor = {coords[0] + i, coords[1] + j, coords[2] + k};
                if (grid.find(neighbor) != grid.end()) {
                    for (const auto& p : grid[neighbor]) {
                        if (distance(p, point) < minDistance) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
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