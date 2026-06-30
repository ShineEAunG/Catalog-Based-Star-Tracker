#pragma once

#include <vector>
#include <opencv2/opencv.hpp>
#include "Parameters.h"


ReprojectionResult VerifyCandidateByReprojection(
    const PyramidCandidate& candidate,
    const std::vector<CameraStar>& cameraStars,
    const std::vector<cv::Point2d>& detectedCentroids,
    const std::vector<CatalogStar>& catalogStars,
    const CameraModel& camera,
    double pixelThreshold,
    double magnitudeLimit,
    bool method);

ReprojectionResult SelectBestCandidateByReprojection(
    const std::vector<PyramidCandidate>& candidates,
    const std::vector<CameraStar>& cameraStars,
    const std::vector<cv::Point2d>& detectedCentroids,
    const std::vector<CatalogStar>& catalogStars,
    const CameraModel& camera,
    int topK,
    double pixelThreshold,
    double magnitudeLimit);

void PrintReprojectionResult(
    const ReprojectionResult& result);