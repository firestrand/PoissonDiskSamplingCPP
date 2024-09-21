# Poisson Disk Sampling in Arbitrary Dimensions

This project provides an implementation of Poisson Disk Sampling in arbitrary dimensions using C++. Poisson Disk Sampling is a technique used to generatePoints a set of points that are evenly distributed with a minimum euclidean_distance between them. This can be useful for generating random samples in multidimensional spaces for various applications, such as graphics, simulation, and scientific computing.

## Features

- **Support for arbitrary dimensions**: The implementation allows for generating points in spaces of any dimensionality.
- **Spatial hashing**: Efficiently checks point validity using a spatial grid and hashing.
- **Configurable minimum euclidean_distance**: Ensures that all points are separated by at least the specified minimum euclidean_distance.
- **Point validation**: Includes utilities to validate that the generated points maintain the minimum euclidean_distance and are evenly distributed across the space.

## Requirements

The implementation uses the following libraries:

- **Boost**: For unordered maps and geometry algorithms.
- **Standard C++ Libraries**: For random number generation, vector manipulation, and I/O.

Make sure to install Boost using your system's package manager or a tool like `vcpkg`.

### Dependencies

- Boost (`boost::unordered_map`)
- C++11 or higher (for random number generation and STL support)

## Installation

Clone this repository and compile the project using a C++ compiler. The `Boost` library is required, and it can be installed using a package manager like `vcpkg` or `apt-get`.

1. Install `Boost` (e.g., with `vcpkg`):
   ```bash
   vcpkg install boost
   ```
2. Compile the project:
   ```bash
   g++ -std=c++11 main.cpp PoissonDiskSampling.cpp -o poisson -lboost_system -lboost_filesystem
   ```

## Usage

### Example

The example provided demonstrates how to generatePoints 100 points in a 4D space using Poisson Disk Sampling.

```cpp
#include "PoissonDiskSampling.h"

int main() {
    int dimensions = 4;  // 4D space
    double minDistance = 1.0;  // Minimum euclidean_distance between points
    std::vector<double> lowerBounds = {0.0, 0.0, 0.0, 0.0};  // Lower bounds for the space
    std::vector<double> upperBounds = {10.0, 10.0, 10.0, 10.0};  // Upper bounds for the space
    int numPoints = 100;  // Number of points to generatePoints

    PoissonDiskSampling sampler(dimensions, minDistance, lowerBounds, upperBounds);
    std::vector<std::vector<double>> points = sampler.generatePoints(numPoints);

    // Output the generated points
    for (const auto& point : points) {
        for (const auto& coord : point) {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
```

### Validating Results

The implementation includes validation functions to ensure that:

- No points are closer than the minimum euclidean_distance (`validateMinDistance`).
- Points are evenly distributed across the grid space (`validateEvenDistribution`).

Both validations are performed in the example code and output to the console.

### Grid-Based Spatial Hashing

The implementation uses a grid-based spatial hashing technique to efficiently check whether a newly generated point satisfies the minimum euclidean_distance requirement from all previously generated points. This allows the algorithm to work in higher dimensions while maintaining performance.

## API

### Class: `PoissonDiskSampling`

#### Constructor
```cpp
PoissonDiskSampling(int dims, double minDist, const std::vector<double>& lower, const std::vector<double>& upper);
```
- **dims**: Number of dimensions for the space.
- **minDist**: Minimum euclidean_distance between points.
- **lower**: Lower bounds for each dimension.
- **upper**: Upper bounds for each dimension.

#### Method: `generatePoints`
```cpp
std::vector<std::vector<double>> generatePoints(int numPoints, int maxAttempts = 30);
```
- **numPoints**: Number of points to generatePoints.
- **maxAttempts**: Maximum maxAttempts to generatePoints a point near a valid one.

Returns a vector of points in the multidimensional space.

## How It Works

The algorithm proceeds as follows:

1. **Grid Setup**: The space is divided into a grid, where each cell has size proportional to the minimum euclidean_distance.
2. **Initial Point Generation**: A random point is generated within the defined space.
3. **Point Validation**: New points are generated and validated to ensure they are not closer than the minimum euclidean_distance to any existing points.
4. **Neighbor Search**: Neighboring grid cells are searched to ensure the minimum euclidean_distance requirement is met.
5. **Point Acceptance**: Valid points are added to the grid and the list of points.
6. **Termination**: The algorithm terminates after a set number of maxAttempts or when the required number of points is generated.

## License

This project is licensed under the MIT License.

---

This `README.md` was generated by a LLM.