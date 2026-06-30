//preprocessing pipeline for image processing, 
//including loading images, 
//converting to grayscale, 
//thresholding to binary, 
//and finding contours and centroids.

#include "PreprocessingImage.h"
#include <algorithm>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

bool CheckFolderExists(const std::string& folderName)
{
    return fs::exists(folderName);
}

std::vector<std::string> LoadImagePaths(
    const std::string& folderPath)
{
    std::vector<std::string> imagePaths;

    if (!fs::exists(folderPath))
    {
        std::cout << "Folder to load images is not found: "
            << folderPath << std::endl;
        return imagePaths;
    }

    for (const auto& image : fs::directory_iterator(folderPath))
    {
        if (!image.is_regular_file())
            continue;

        std::string filePath = image.path().string();
        std::string fileExtension = image.path().extension().string();

        if (fileExtension == ".png" ||
            fileExtension == ".jpg" ||
            fileExtension == ".jpeg" ||
            fileExtension == ".bmp")
        {
            std::cout << "Found image in " << filePath << std::endl;
            imagePaths.push_back(filePath);
        }
    }

    return imagePaths;
}

bool ConvertToGrayScale(
    const std::string& imagePath,
    const std::string& grayImagePath)
{
    cv::Mat inputImage = cv::imread(imagePath, cv::IMREAD_COLOR);

    if (inputImage.empty())
    {
        std::cerr << "ERROR: Cannot load input image: "
            << imagePath << std::endl;
        return false;
    }

    cv::Mat grayImage;
    cv::cvtColor(inputImage, grayImage, cv::COLOR_BGR2GRAY);

    if (!cv::imwrite(grayImagePath, grayImage))
    {
        std::cerr << "ERROR: Cannot save grayscale image: "
            << grayImagePath << std::endl;
        return false;
    }

    return true;
}

bool TurnGrayToBinary(
    const std::string& grayImagePath,
    int thresholdLevel,
    const std::string& binaryImagePath)
{
    cv::Mat grayImage =
        cv::imread(grayImagePath, cv::IMREAD_GRAYSCALE);

    if (grayImage.empty())
    {
        std::cerr << "ERROR: Cannot load grayscale image: "
            << grayImagePath << std::endl;
        return false;
    }

    cv::Mat binaryImage;

    cv::threshold(
        grayImage,
        binaryImage,
        thresholdLevel,
        255,
        cv::THRESH_BINARY
    );

    if (!cv::imwrite(binaryImagePath, binaryImage))
    {
        std::cerr << "ERROR: Cannot save binary image: "
            << binaryImagePath << std::endl;
        return false;
    }

    return true;
}

std::tuple<std::string, std::string> PreprocessInputImage(
    const std::string& imagePath,
    const std::string& folder,
    const std::string& subFolder,
    int thresholdLevel)
{
    std::cout << "***********Pre-Processing On Image***********" << std::endl;

    if (!CheckFolderExists(subFolder))
    {
        fs::create_directories(subFolder);
    }

    std::vector<std::string> imagePaths =
        LoadImagePaths(folder);

    std::cout << "Images found = "
        << imagePaths.size() << std::endl;

    std::string grayImagePath =
        subFolder + "gray_star_image.png";

    std::string binaryImageName =
        "binary_star_image.png";

    std::string binaryImagePath =
        subFolder + binaryImageName;

    bool grayOk =
        ConvertToGrayScale(
            imagePath,
            grayImagePath
        );

    if (!grayOk)
        return { "", "" };

    bool binaryOk =
        TurnGrayToBinary(
            grayImagePath,
            thresholdLevel,
            binaryImagePath
        );

    if (!binaryOk)
        return { "", "" };

    return {
        binaryImagePath,
        binaryImageName
    };
}

std::tuple<std::vector<cv::Point2d>, int, int>
FindContoursCentroids(
    const std::string& binaryImagePath,
    int maxStars,
    double minimumArea)
{
    cv::Mat binaryImage =
        cv::imread(binaryImagePath, cv::IMREAD_GRAYSCALE);

    if (binaryImage.empty())
    {
        std::cerr << "ERROR: Cannot load binary image: "
            << binaryImagePath << std::endl;

        return { {}, 0, 0 };
    }

    std::vector<std::vector<cv::Point>> contours;

    cv::findContours(
        binaryImage,
        contours,
        cv::RETR_EXTERNAL,
        cv::CHAIN_APPROX_SIMPLE
    );

    struct DetectedBlob
    {
        cv::Point2d centroid;
        double area;
    };

    std::vector<DetectedBlob> blobs;

    for (const auto& contour : contours)
    {
        double area =
            cv::contourArea(contour);

        if (area < minimumArea)
            continue;

        cv::Moments m =
            cv::moments(contour);

        if (m.m00 == 0.0)
            continue;

        double cx =
            m.m10 / m.m00;

        double cy =
            m.m01 / m.m00;

        blobs.push_back(
            {
                cv::Point2d(cx, cy),
                area
            });
    }

    std::sort(
        blobs.begin(),
        blobs.end(),
        [](const DetectedBlob& a, const DetectedBlob& b)
        {
            return a.area > b.area;
        }
    );

    if (maxStars > 0 &&
        static_cast<int>(blobs.size()) > maxStars)
    {
        blobs.resize(maxStars);
    }

    std::vector<cv::Point2d> centroids;

    for (const auto& blob : blobs)
    {
        centroids.push_back(blob.centroid);
    }

    std::cout << "Detected usable stars = "
        << centroids.size()
        << " / raw blobs = "
        << contours.size()
        << std::endl;

    return {
        centroids,
        binaryImage.cols,
        binaryImage.rows
    };
}

//std::tuple<std::vector<cv::Point2d>, int, int>
//FindContoursCentroids(
//    const std::string& binaryImagePath)
//{
//    cv::Mat binaryImage =
//        cv::imread(binaryImagePath, cv::IMREAD_GRAYSCALE);
//
//    if (binaryImage.empty())
//    {
//        std::cerr << "ERROR: Cannot load binary image: "
//            << binaryImagePath << std::endl;
//
//        return { {}, 0, 0 };
//    }
//
//    std::vector<std::vector<cv::Point>> contours;
//
//    cv::findContours(
//        binaryImage,
//        contours,
//        cv::RETR_EXTERNAL,
//        cv::CHAIN_APPROX_SIMPLE
//    );
//
//    std::vector<cv::Point2d> centroids;
//
//    const double minimumArea = 1.0;
//
//    for (const auto& contour : contours)
//    {
//        double area = cv::contourArea(contour);
//
//        if (area < minimumArea)
//            continue;
//
//        cv::Moments m = cv::moments(contour);
//
//        if (m.m00 == 0.0)
//            continue;
//
//        double cx = m.m10 / m.m00;
//        double cy = m.m01 / m.m00;
//
//        centroids.push_back(cv::Point2d(cx, cy));
//    }
//
//    std::cout << "Detected stars = "
//        << centroids.size() << std::endl;
//
//    return {
//        centroids,
//        binaryImage.cols,
//        binaryImage.rows
//    };
//}

void ShowCentroidsOnConsole(
    const std::vector<cv::Point2d>& centroids)
{
    std::cout << "Centroids found at:\n";

    for (int i = 0; i < static_cast<int>(centroids.size()); i++)
    {
        std::cout << "Star " << i + 1
            << " = ("
            << centroids[i].x
            << ", "
            << centroids[i].y
            << ")"
            << std::endl;
    }
}