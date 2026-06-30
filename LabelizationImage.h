#pragma once

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include "Parameters.h"

bool SaveDetectedStarLabelImage(
    const std::string& inputImagePath,
    const std::vector<cv::Point2d>& centroids,
    const std::string& outputImagePath);

bool SaveCandidateStarLabelImage(
    const std::string& inputImagePath,
    const std::vector<cv::Point2d>& centroids,
    const PyramidCandidate& candidate,
    const std::string& outputImagePath);