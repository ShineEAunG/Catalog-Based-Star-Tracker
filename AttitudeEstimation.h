#pragma once

#include <vector>
#include <opencv2/opencv.hpp>
#include "Parameters.h"

AttitudeResult EstimateAttitudeSVD(
    const std::vector<MatchedStar>& matchedStars);

void PrintAttitudeResult(
    const AttitudeResult& result);

Quaternion NormalizeQuaternion(const Quaternion& q);

double QuaternionAngleErrorDeg(
    const Quaternion& qEstimated,
    const Quaternion& qTruth);

std::tuple<double, double, double> QuaternionToEular(Quaternion& q_est);

AttitudeResult EstimateAttitudeQUEST(
    const std::vector<MatchedStar>& matchedStars);

cv::Matx33d QuaternionToRotationMatrix(const Quaternion& q_in);

cv::Matx33d CalculateAttitudeProfileMatrix(const std::vector<MatchedStar>& matchedStars);