#pragma once

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <fstream>
#include <filesystem>
#include "Parameters.h"
#include "CatalogGenerator.h"
#include "PreprocessingImage.h"


std::vector<Star> GetPolarisCostellation();
std::vector<Star> GetCassiopeiaConstellation();
std::vector<Star> GetOrionConstellation();
std::vector<Star> GetBigDipper();
std::vector<Star> GetCygnus();

std::vector<double> RotateVector(
    const Quaternion& q,
    const std::vector<double>& v);

std::vector<CameraTruthStar> ConvertCatalogIntoCameraFrame(
    const std::vector<CatalogStar>& catalogStars,
    const Quaternion& q_truth,
    double maxBrightness,
    int width,
    int height,
    double focalLengthPixels,
    double magnitudeLimit);

std::vector<PixelCoordinate> ConvertIntoPixelCoordinates(
    const std::vector<CameraTruthStar>& cameraStars,
    int width,
    int height,
    double focalLengthPixels);


void LogDataInConsole(
    const std::vector<PixelCoordinate>& coordinates);

void WriteCoordinatesToFile(
    const std::vector<PixelCoordinate>& stars,
    const std::string& fileName,
    const std::string& folder);

void Generate3ColorsConstellation(
    const std::vector<PixelCoordinate>& pixelCoordinates,
    int width,
    int height,
    const std::string& imageName,
    const std::string& folder);

void ShowImage(const std::string& imageName,
    const std::string& windowName,
    const std::string& folder);

