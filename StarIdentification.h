#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <opencv2/opencv.hpp>
#include "Parameters.h"


std::vector<CameraStar> ConvertToCameraUnitVectors(
    const std::vector<cv::Point2d>& centroids,
    double focalLengthPixels,
    int imageWidth,
    int imageHeight);

double CalculateAngleDeg(
    const CameraStar& a,
    const CameraStar& b);

void SaveLabeledIdentifiedPyramidImage(
    const std::string& inputImagePath,
    const std::string& outputImagePath,
    const IdentifiedPyramid& result);

IdentifiedPyramid IdentifyPyramidFromBinaryImage(
    const std::string& binaryImagePath,
    const std::vector<StarPairAngle>& pairDatabase,
    const std::vector<CatalogStar>& catalogStars,
    double focalLengthPixels,
    double toleranceDeg,
    int maxDetectedStars = 15);

StarIdentificationDebugData IdentifyAllPyramidCandidatesForDebugging(
    const std::string& binaryImagePath,
    const std::vector<StarPairAngle>& pairDatabase,
    double focalLengthPixels,
    double toleranceDeg,
    int maxDetectedStars);
