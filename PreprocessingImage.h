#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <opencv2/opencv.hpp>


std::vector<std::string> LoadImagePaths(
    const std::string& folderPath);

bool ConvertToGrayScale(
    const std::string& imagePath,
    const std::string& grayImagePath);

bool TurnGrayToBinary(
    const std::string& grayImagePath,
    int thresholdLevel,
    const std::string& binaryImagePath);

std::tuple<std::string, std::string> PreprocessInputImage(
    const std::string& imagePath,
    const std::string& folder,
    const std::string& subFolder,
    int thresholdLevel);

std::tuple<std::vector<cv::Point2d>, int, int>
FindContoursCentroids(
    const std::string& binaryImagePath,
    int maxStars = 15,
    double minimumArea = 1.0);

//std::tuple<std::vector<cv::Point2d>, int, int>
//FindContoursCentroids(
//    const std::string& binaryImagePath);

void ShowCentroidsOnConsole(
    const std::vector<cv::Point2d>& centroids);