#ifndef FIND_PEAK_H
#define FIND_PEAK_H


double GetMax(const std::vector<double>& data);

double GetMin(const std::vector<double>& data);

std::vector<int32_t> findPeaks(const std::vector<double>& x,
                           double height, int distance);

std::vector<int32_t> diffSignal(const std::vector<int32_t>& signal);

double meanSignal(const std::vector<double>& signal);

int32_t meanSignal(const std::vector<int32_t>& signal);

int32_t modeSignal(const std::vector<int32_t>& nums);

#endif