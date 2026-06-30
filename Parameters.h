#pragma once

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>

struct CatalogStar
{
    int hip = 0;

    double ra_deg = 0.0;
    double dec_deg = 0.0;
    double magnitude = 0.0;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct StarPairAngle
{
    int hipA = 0;
    int hipB = 0;
    double angle_deg = 0.0;
};

struct CameraStar
{
    int index = -1;
    std::string name;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct PyramidCandidate
{
    bool found = false;

    int imageA = -1;
    int imageB = -1;
    int imageC = -1;
    int imageD = -1;

    int hipA = -1;
    int hipB = -1;
    int hipC = -1;
    int hipD = -1;

    double totalAngleErrorDeg = 0.0;
    double averageAngleErrorDeg = 0.0;
    int voteScore = 0;

    std::vector<cv::Point2d> centroids;
};

struct MatchedStar
{
    int imageIndex = -1;
    int hip = -1;

    double bx = 0.0;
    double by = 0.0;
    double bz = 0.0;

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;
};

struct IdentifiedPyramid
{
    bool found = false;

    int imageA = -1;
    int imageB = -1;
    int imageC = -1;
    int imageD = -1;

    int hipA = -1;
    int hipB = -1;
    int hipC = -1;
    int hipD = -1;

    int voteScore = 0;
    double averageAngleErrorDeg = 0.0;

    std::vector<cv::Point2d> centroids;
    std::vector<CameraStar> cameraStars;

    std::vector<MatchedStar> matchedStars;

    int imageWidth = 0;
    int imageHeight = 0;
};

struct StarIdentificationDebugData
{
    bool success = false;

    std::vector<cv::Point2d> centroids;
    std::vector<CameraStar> cameraStars;
    std::vector<PyramidCandidate> candidates;

    int imageWidth = 0;
    int imageHeight = 0;
};

struct Quaternion
{
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct AttitudeResult
{
    bool success = false;

    cv::Matx33d R = cv::Matx33d::eye();
    Quaternion q;

    int usedStars = 0;
    double rmsResidualDeg = 0.0;
};

struct CameraModel
{
    double fx = 0.0;
    double fy = 0.0;
    double cx = 0.0;
    double cy = 0.0;

    int width = 0;
    int height = 0;
};

struct ReprojectionResult
{
    bool success = false;

    PyramidCandidate candidate;
    AttitudeResult attitude;

    int inlierCount = 0;
    double rmsPixelError = 1.0e9;
    double rmsAngularErrorDeg = 1.0e9;

    std::vector<MatchedStar> matchedStars;
};

struct StarIdentificationResult
{
    bool success = false;

    PyramidCandidate bestCandidate;

    std::vector<cv::Point2d> centroids;
    std::vector<CameraStar> cameraStars;
    std::vector<MatchedStar> matchedStars;

    int imageWidth = 0;
    int imageHeight = 0;
};

struct TruthStar
{
    int hip = 0;
    double u = 0.0;
    double v = 0.0;
    double brightness = 0.0;
};

struct ValidationSummary
{
    int truthStars = 0;
    int selectedPyramidStars = 0;
    int correctPyramidStars = 0;

    double svdResidualArcsec = 0.0;
    double questResidualArcsec = 0.0;

    double svdAttitudeErrorDeg = -1.0;
    double questAttitudeErrorDeg = -1.0;

    bool pass = false;
};

struct ProjectedCatalogStar
{
    int hip = 0;

    double u = 0.0;
    double v = 0.0;

    cv::Vec3d cameraVector;
};

bool CheckFolderExists(const std::string& folderName);

struct Star
{
    std::string name;
    double rann, dec, magnitude;
};

struct CameraTruthStar
{
    int hip;
    double x;
    double y;
    double z;
    double brightness;
};

struct PixelCoordinate
{
    int hip;
    double u;
    double v;
    double brightness;
};