#include "LabelizationImage.h"

#include <iostream>
#include <unordered_map>
#include <sstream>
#include <algorithm>

static cv::Mat LoadImageForDrawing(
    const std::string& inputImagePath)
{
    cv::Mat image =
        cv::imread(inputImagePath, cv::IMREAD_UNCHANGED);

    if (image.empty())
    {
        std::cerr << "ERROR: Cannot load image for labeling: "
            << inputImagePath << "\n";
        return cv::Mat();
    }

    cv::Mat displayImage;

    if (image.channels() == 1)
    {
        cv::cvtColor(image, displayImage, cv::COLOR_GRAY2BGR);
    }
    else if (image.channels() == 3)
    {
        displayImage = image.clone();
    }
    else if (image.channels() == 4)
    {
        cv::cvtColor(image, displayImage, cv::COLOR_BGRA2BGR);
    }
    else
    {
        std::cerr << "ERROR: Unsupported image channel count.\n";
        return cv::Mat();
    }

    return displayImage;
}

static cv::Point MakeSafeTextPosition(
    const cv::Point2d& centroid,
    const std::string& text,
    const cv::Size& imageSize,
    double fontScale,
    int thickness)
{
    int baseline = 0;

    cv::Size textSize =
        cv::getTextSize(
            text,
            cv::FONT_HERSHEY_SIMPLEX,
            fontScale,
            thickness,
            &baseline);

    int x = static_cast<int>(centroid.x) + 8;
    int y = static_cast<int>(centroid.y) - 8;

    if (x + textSize.width >= imageSize.width)
    {
        x = static_cast<int>(centroid.x) - textSize.width - 8;
    }

    if (y - textSize.height < 0)
    {
        y = static_cast<int>(centroid.y) + textSize.height + 8;
    }

    x = std::max(0, std::min(x, imageSize.width - 1));
    y = std::max(textSize.height, std::min(y, imageSize.height - 1));

    return cv::Point(x, y);
}

static void DrawTextWithOutline(
    cv::Mat& image,
    const std::string& text,
    const cv::Point& position,
    double fontScale)
{
    cv::putText(
        image,
        text,
        position,
        cv::FONT_HERSHEY_SIMPLEX,
        fontScale,
        cv::Scalar(0, 0, 0),
        3,
        cv::LINE_AA);

    cv::putText(
        image,
        text,
        position,
        cv::FONT_HERSHEY_SIMPLEX,
        fontScale,
        cv::Scalar(0, 255, 255),
        1,
        cv::LINE_AA);
}

static void DrawOneStarLabel(
    cv::Mat& image,
    const cv::Point2d& centroid,
    const std::string& label,
    bool isSelected)
{
    cv::Point center(
        static_cast<int>(std::round(centroid.x)),
        static_cast<int>(std::round(centroid.y)));

    cv::Scalar circleColor =
        isSelected ? cv::Scalar(0, 0, 255) : cv::Scalar(0, 255, 0);

    int circleRadius =
        isSelected ? 7 : 5;

    int circleThickness =
        isSelected ? 2 : 1;

    cv::circle(
        image,
        center,
        circleRadius,
        circleColor,
        circleThickness,
        cv::LINE_AA);

    cv::drawMarker(
        image,
        center,
        circleColor,
        cv::MARKER_CROSS,
        8,
        1,
        cv::LINE_AA);

    double fontScale = 0.45;

    cv::Point textPosition =
        MakeSafeTextPosition(
            centroid,
            label,
            image.size(),
            fontScale,
            1);

    DrawTextWithOutline(
        image,
        label,
        textPosition,
        fontScale);
}

bool SaveDetectedStarLabelImage(
    const std::string& inputImagePath,
    const std::vector<cv::Point2d>& centroids,
    const std::string& outputImagePath)
{
    cv::Mat image =
        LoadImageForDrawing(inputImagePath);

    if (image.empty())
        return false;

    for (int i = 0; i < static_cast<int>(centroids.size()); i++)
    {
        std::ostringstream label;
        label << "S" << i;

        DrawOneStarLabel(
            image,
            centroids[i],
            label.str(),
            false);
    }

    bool saved =
        cv::imwrite(outputImagePath, image);

    if (!saved)
    {
        std::cerr << "ERROR: Cannot save labeled image: "
            << outputImagePath << "\n";
        return false;
    }

    std::cout << "Labeled detected-star image saved: "
        << outputImagePath << "\n";

    return true;
}

static std::unordered_map<int, int> BuildCandidateHipMap(
    const PyramidCandidate& candidate)
{
    std::unordered_map<int, int> hipByImageIndex;

    hipByImageIndex[candidate.imageA] = candidate.hipA;
    hipByImageIndex[candidate.imageB] = candidate.hipB;
    hipByImageIndex[candidate.imageC] = candidate.hipC;
    hipByImageIndex[candidate.imageD] = candidate.hipD;

    return hipByImageIndex;
}

bool SaveCandidateStarLabelImage(
    const std::string& inputImagePath,
    const std::vector<cv::Point2d>& centroids,
    const PyramidCandidate& candidate,
    const std::string& outputImagePath)
{
    cv::Mat image =
        LoadImageForDrawing(inputImagePath);

    if (image.empty())
        return false;

    std::unordered_map<int, int> hipByImageIndex =
        BuildCandidateHipMap(candidate);

    for (int i = 0; i < static_cast<int>(centroids.size()); i++)
    {
        std::ostringstream label;

        auto it =
            hipByImageIndex.find(i);

        bool isSelected =
            it != hipByImageIndex.end();

        if (isSelected)
        {
            label << "S" << i << " / HIP " << it->second;
        }
        else
        {
            label << "S" << i;
        }

        DrawOneStarLabel(
            image,
            centroids[i],
            label.str(),
            isSelected);
    }

    bool saved =
        cv::imwrite(outputImagePath, image);

    if (!saved)
    {
        std::cerr << "ERROR: Cannot save candidate-labeled image: "
            << outputImagePath << "\n";
        return false;
    }

    std::cout << "Candidate-labeled image saved: "
        << outputImagePath << "\n";

    return true;
}