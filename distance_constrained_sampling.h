//
// Created by Travis Silvers on 9/20/24.
//

#ifndef POISSONDISKSAMPLINGCPP_DISTANCE_CONSTRAINED_SAMPLING_H
#define POISSONDISKSAMPLINGCPP_DISTANCE_CONSTRAINED_SAMPLING_H

#include "sampling_interface.h"
#include <vector>
#include <unordered_map>
#include <random>
#include <boost/unordered_map.hpp>

class DistanceConstrainedSampling : public SamplingInterface{
private:
    int dimensions;
    int numPoints;
    int maxAttempts;
    double minDistance;
    double cellSize;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    std::vector<std::vector<double>> points;
    boost::unordered_map<std::vector<int>, std::vector<std::vector<double>>> grid;  // Grid for spatial hashing

    std::mt19937 gen;
    std::uniform_real_distribution<> uniform;

    double calculateCellSize();
    std::vector<double> generateRandomPoint();
    std::vector<int> gridCoords(const std::vector<double>& point);
    bool isValidPoint(const std::vector<double>& point);
    void generateNeighbors(std::vector<int> &coords, int dim, std::vector<std::vector<int>> &neighbors);

public:
    DistanceConstrainedSampling(int dims, int numPoints, int maxAttempts ,const std::vector<double>& lower, const std::vector<double>& upper);
    std::vector<std::vector<double>> generatePoints();
};


#endif //POISSONDISKSAMPLINGCPP_DISTANCE_CONSTRAINED_SAMPLING_H
