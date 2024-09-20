#ifndef POISSON_DISK_SAMPLING_H
#define POISSON_DISK_SAMPLING_H

#include <vector>
#include <unordered_map>
#include <random>

struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const;
};

class PoissonDiskSampling {
private:
    int dimensions;
    double minDistance;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    std::vector<std::vector<double>> points;
    std::unordered_map<std::vector<int>, std::vector<std::vector<double>>, VectorHash> grid;
    double cellSize;

    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

    double distance(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> generateRandomPoint();
    std::vector<int> gridCoords(const std::vector<double>& point);
    bool isValidPoint(const std::vector<double>& point);

public:
    PoissonDiskSampling(int dims, double minDist, const std::vector<double>& lower, const std::vector<double>& upper);
    std::vector<std::vector<double>> generatePoints(int numPoints, int maxAttempts = 30);
};

#endif  // POISSON_DISK_SAMPLING_H