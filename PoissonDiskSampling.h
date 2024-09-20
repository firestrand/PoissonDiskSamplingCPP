#ifndef POISSON_DISK_SAMPLING_H
#define POISSON_DISK_SAMPLING_H

#include <vector>
#include <unordered_map>
#include <random>
#include <boost/unordered_map.hpp>

class PoissonDiskSampling {
private:
    int dimensions;
    double minDistance;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    std::vector<std::vector<double>> points;
    boost::unordered_map<std::vector<int>, std::vector<std::vector<double>>> grid;  // Grid for spatial hashing
    double cellSize;

    std::mt19937 gen;
    std::uniform_real_distribution<> uniform;

    static double distance(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> generateRandomPoint();
    std::vector<int> gridCoords(const std::vector<double>& point);
    bool isValidPoint(const std::vector<double>& point);

public:
    PoissonDiskSampling(int dims, double minDist, const std::vector<double>& lower, const std::vector<double>& upper);
    std::vector<std::vector<double>> generatePoints(int numPoints, int maxAttempts = 30);
    void generateNeighbors(std::vector<int> &coords, int dim, std::vector<std::vector<int>> &neighbors);
};

#endif  // POISSON_DISK_SAMPLING_H