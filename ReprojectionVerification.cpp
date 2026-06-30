#include "ReprojectionVerification.h"
#include "AttitudeEstimation.h"
#include <iostream>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>

static cv::Vec3d NormalizeVectorLocal(const cv::Vec3d& v)
{
    double n = cv::norm(v);

    if (n <= 0.0)
        return cv::Vec3d(0.0, 0.0, 0.0);

    return v / n;
}

static double AngularErrorDegLocal(
    const cv::Vec3d& a,
    const cv::Vec3d& b)
{
    cv::Vec3d an = NormalizeVectorLocal(a);
    cv::Vec3d bn = NormalizeVectorLocal(b);

    double dot =
        an[0] * bn[0] +
        an[1] * bn[1] +
        an[2] * bn[2];

    dot = std::clamp(dot, -1.0, 1.0);

    return std::acos(dot) * 180.0 / CV_PI;
}

static std::unordered_map<int, CatalogStar> BuildCatalogMap(
    const std::vector<CatalogStar>& catalogStars)
{
    std::unordered_map<int, CatalogStar> catalogMap;

    for (const auto& star : catalogStars)
    {
        catalogMap[star.hip] = star;
    }

    return catalogMap;
}

static bool FindCameraStarByIndex(
    const std::vector<CameraStar>& cameraStars,
    int imageIndex,
    CameraStar& output)
{
    for (const auto& star : cameraStars)
    {
        if (star.index == imageIndex)
        {
            output = star;
            return true;
        }
    }

    return false;
}

static std::vector<std::pair<int, int>> GetCandidatePairs(
    const PyramidCandidate& candidate)
{
    std::vector<std::pair<int, int>> pairs;

    pairs.push_back({ candidate.imageA, candidate.hipA });
    pairs.push_back({ candidate.imageB, candidate.hipB });
    pairs.push_back({ candidate.imageC, candidate.hipC });
    pairs.push_back({ candidate.imageD, candidate.hipD });

    return pairs;
}

static std::vector<MatchedStar> BuildMatchedStarsFromCandidate(
    const PyramidCandidate& candidate,
    const std::vector<CameraStar>& cameraStars,
    const std::unordered_map<int, CatalogStar>& catalogMap)
{
    std::vector<MatchedStar> matchedStars;

    std::vector<std::pair<int, int>> pairs =
        GetCandidatePairs(candidate);

    for (const auto& pair : pairs)
    {
        int imageIndex = pair.first;
        int hip = pair.second;

        CameraStar cameraStar;

        if (!FindCameraStarByIndex(
            cameraStars,
            imageIndex,
            cameraStar))
        {
            continue;
        }

        auto it = catalogMap.find(hip);

        if (it == catalogMap.end())
        {
            continue;
        }

        const CatalogStar& catalogStar = it->second;

        MatchedStar matched;

        matched.imageIndex = imageIndex;
        matched.hip = hip;

        matched.bx = cameraStar.x;
        matched.by = cameraStar.y;
        matched.bz = cameraStar.z;

        matched.rx = catalogStar.x;
        matched.ry = catalogStar.y;
        matched.rz = catalogStar.z;

        matchedStars.push_back(matched);
    }

    return matchedStars;
}

static bool ProjectCatalogStarToImage(
    const CatalogStar& catalogStar,
    const cv::Matx33d& R,
    const CameraModel& camera,
    ProjectedCatalogStar& projected)
{
    cv::Vec3d r(
        catalogStar.x,
        catalogStar.y,
        catalogStar.z);

    r = NormalizeVectorLocal(r);

    // R maps catalog/inertial vector to camera/body vector.
    // b = R * r
    cv::Vec3d b = R * r;
    b = NormalizeVectorLocal(b);

    if (b[2] <= 0.0)
        return false;

    double u =
        camera.cx +
        camera.fx * (b[0] / b[2]);

    double v =
        camera.cy -
        camera.fy * (b[1] / b[2]);

    if (u < 0.0 || u >= camera.width)
        return false;

    if (v < 0.0 || v >= camera.height)
        return false;

    projected.hip = catalogStar.hip;
    projected.u = u;
    projected.v = v;
    projected.cameraVector = b;

    return true;
}

static cv::Vec3d CentroidToCameraVector(
    const cv::Point2d& centroid,
    const CameraModel& camera)
{
    double x = (centroid.x - camera.cx) / camera.fx;
    double y = (camera.cy - centroid.y) / camera.fy;
    double z = 1.0;

    return NormalizeVectorLocal(cv::Vec3d(x, y, z));
}

ReprojectionResult VerifyCandidateByReprojection(
    const PyramidCandidate& candidate,
    const std::vector<CameraStar>& cameraStars,
    const std::vector<cv::Point2d>& detectedCentroids,
    const std::vector<CatalogStar>& catalogStars,
    const CameraModel& camera,
    double pixelThreshold,
    double magnitudeLimit)
{
    ReprojectionResult result;
    result.candidate = candidate;

    std::unordered_map<int, CatalogStar> catalogMap =
        BuildCatalogMap(catalogStars);

    std::vector<MatchedStar> matchedStars =
        BuildMatchedStarsFromCandidate(
            candidate,
            cameraStars,
            catalogMap);

    result.matchedStars = matchedStars;

    if (matchedStars.size() < 3)
    {
        return result;
    }
    
	// change to EstimateAttitudeQUEST(matchedStars) back and forth according to user input method
    AttitudeResult attitude = EstimateAttitudeSVD(matchedStars);
    //AttitudeResult attitude = EstimateAttitudeQUEST(matchedStars);

    if (!attitude.success)
    {
        return result;
    }

    result.attitude = attitude;

    std::vector<bool> usedDetected(
        detectedCentroids.size(),
        false);

    double sumPixelSquared = 0.0;
    double sumAngularSquared = 0.0;

    int inliers = 0;

    for (const auto& catalogStar : catalogStars)
    {
        if (catalogStar.magnitude > magnitudeLimit)
            continue;

        ProjectedCatalogStar projected;

        bool ok =
            ProjectCatalogStarToImage(
                catalogStar,
                attitude.R,
                camera,
                projected);

        if (!ok)
            continue;

        int bestIndex = -1;
        double bestDistance = std::numeric_limits<double>::max();

        for (int i = 0; i < static_cast<int>(detectedCentroids.size()); i++)
        {
            if (usedDetected[i])
                continue;

            double dx =
                detectedCentroids[i].x -
                projected.u;

            double dy =
                detectedCentroids[i].y -
                projected.v;

            double distance =
                std::sqrt(dx * dx + dy * dy);

            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestIndex = i;
            }
        }

        if (bestIndex >= 0 && bestDistance <= pixelThreshold)
        {
            usedDetected[bestIndex] = true;
            inliers++;

            sumPixelSquared +=
                bestDistance * bestDistance;

            cv::Vec3d measuredVector =
                CentroidToCameraVector(
                    detectedCentroids[bestIndex],
                    camera);

            double angularErrorDeg =
                AngularErrorDegLocal(
                    measuredVector,
                    projected.cameraVector);

            sumAngularSquared +=
                angularErrorDeg * angularErrorDeg;
        }
    }

    result.inlierCount = inliers;

    if (inliers > 0)
    {
        result.rmsPixelError =
            std::sqrt(sumPixelSquared / inliers);

        result.rmsAngularErrorDeg =
            std::sqrt(sumAngularSquared / inliers);
    }

    result.success = true;

    return result;
}

ReprojectionResult SelectBestCandidateByReprojection(
    const std::vector<PyramidCandidate>& candidates,
    const std::vector<CameraStar>& cameraStars,
    const std::vector<cv::Point2d>& detectedCentroids,
    const std::vector<CatalogStar>& catalogStars,
    const CameraModel& camera,
    int topK,
    double pixelThreshold,
    double magnitudeLimit)
{
    ReprojectionResult bestResult;
    const int minimumUsefulInliers = 5;
    int count =
        std::min(
            topK,
            static_cast<int>(candidates.size()));

    for (int i = 0; i < count; i++)
    {
        ReprojectionResult current =
            VerifyCandidateByReprojection(
                candidates[i],
                cameraStars,
                detectedCentroids,
                catalogStars,
                camera,
                pixelThreshold,
                magnitudeLimit);

        if (!current.success)
            continue;
        if (current.inlierCount < minimumUsefulInliers)
            continue;

        bool replace = false;

        if (!bestResult.success)
        {
            replace = true;
        }
        else if (current.inlierCount > bestResult.inlierCount)
        {
            replace = true;
        }
        else if (current.inlierCount == bestResult.inlierCount &&
            current.rmsPixelError < bestResult.rmsPixelError)
        {
            replace = true;
        }
        else if (current.inlierCount == bestResult.inlierCount &&
            std::abs(current.rmsPixelError - bestResult.rmsPixelError) < 1.0e-9 &&
            current.attitude.rmsResidualDeg < bestResult.attitude.rmsResidualDeg)
        {
            replace = true;
        }
        else if (current.inlierCount == bestResult.inlierCount &&
            std::abs(current.rmsPixelError - bestResult.rmsPixelError) < 1.0e-9 &&
            std::abs(current.attitude.rmsResidualDeg - bestResult.attitude.rmsResidualDeg) < 1.0e-12 &&
            current.candidate.voteScore > bestResult.candidate.voteScore)
        {
            replace = true;
        }

        if (replace)
        {
            bestResult = current;
        }
    }

    return bestResult;
}

void PrintReprojectionResult(
    const ReprojectionResult& result)
{
    std::string choosedMethod;

    if (!result.success)
    {
        std::cout << "\nReprojection verification failed.\n";
        return;
    }

    const PyramidCandidate& c = result.candidate;

    std::cout << "\n***********Reprojection Verification Result***********\n";

    std::cout << "Selected candidate:\n";

    std::cout << "Image Stars ["
        << c.imageA << ", "
        << c.imageB << ", "
        << c.imageC << ", "
        << c.imageD << "] -> HIP ["
        << c.hipA << ", "
        << c.hipB << ", "
        << c.hipC << ", "
        << c.hipD << "]\n";

    std::cout << "Vote score = "
        << c.voteScore << "\n";

    std::cout << "Average pyramid angle error = "
        << c.averageAngleErrorDeg
        << " deg\n";

    std::cout << "Reprojection inliers = "
        << result.inlierCount << "\n";

    std::cout << "RMS reprojection pixel error = "
        << result.rmsPixelError
        << " px\n";

    std::cout << "RMS reprojection angular error = "
        << result.rmsAngularErrorDeg
        << " deg = "
        << result.rmsAngularErrorDeg * 3600.0
        << " arcsec\n";

    std::cout << "Attitude solver residual = "
        << result.attitude.rmsResidualDeg
        << " deg = "
        << result.attitude.rmsResidualDeg * 3600.0
        << " arcsec\n";
}