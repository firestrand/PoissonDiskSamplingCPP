#ifndef SAMPLINGINTERFACE_H
#define SAMPLINGINTERFACE_H

#include <vector>

class SamplingInterface {
public:
    virtual ~SamplingInterface() = default;
    virtual std::vector<std::vector<double>> generatePoints() = 0;
};

#endif  // SAMPLINGINTERFACE_H