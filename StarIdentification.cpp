//this session includes the algorithm for star identification 
//from a binary image, which involves detecting stars, 
//converting their positions to unit vectors, 
//searching for matching star pairs in a catalog, 
//and selecting the best matching pyramid of stars based on voting scores and angle errors. 
//The code also includes functions for calculating angles between stars, 
//building matched stars for attitude estimation, 
//and saving labeled images with identified stars for visualization

#include "StarIdentification.h"
#include "PreprocessingImage.h"
#include "Pattern_Generator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "CatalogGenerator.h"

constexpr double PI = 3.14159265358979323846;

//stores the 6 angular distances between 4 detected image stars
struct ObservedPyramidAngles
{
    double AB = 0.0;
    double AC = 0.0;
    double AD = 0.0;
    double BC = 0.0;
    double BD = 0.0;
    double CD = 0.0;
};

//a unique key for a HIP pair for fast lookup (a+b=b+a)
static std::uint64_t MakePairKey(int hip1, int hip2)
{
    int a = std::min(hip1, hip2);
    int b = std::max(hip1, hip2);

    return (static_cast<std::uint64_t>(a) << 32) |
        static_cast<std::uint32_t>(b);
}

//counter checker if there is known HIP in the pair and returns the other HIP if found
//AB => AC => AD then this function is used to find the HIP in BC, BD, CD pairs for consistency check
static bool GetOtherStar(
    const StarPairAngle& pair,
    int knownHip,
    int& otherHip)
{
    if (pair.hipA == knownHip)
    {
        otherHip = pair.hipB;
        return true;
    }

    if (pair.hipB == knownHip)
    {
        otherHip = pair.hipA;
        return true;
    }

    return false;
}

//This function converts a list of candidate catalog pairs into a fast lookup table
static std::unordered_map<std::uint64_t, double>
BuildPairAngleLookup(const std::vector<StarPairAngle>& pairs)
{
    std::unordered_map<std::uint64_t, double> lookup;

    lookup.reserve(pairs.size());

    for (const auto& p : pairs)
    {
        std::uint64_t key = MakePairKey(p.hipA, p.hipB);
        lookup[key] = p.angle_deg;
    }

    return lookup;
}


//trying to get the angle between two HIPs from the lookup table found or not
static bool TryGetPairAngle(
    const std::unordered_map<std::uint64_t, double>& lookup,
    int hip1,
    int hip2,
    double& angleDeg)
{
    std::uint64_t key = MakePairKey(hip1, hip2);

    auto it = lookup.find(key);

    if (it == lookup.end())
        return false;

    angleDeg = it->second;
    return true;
}


//create unique string for one candidate pyramid
static std::string MakeCandidateKey(
    int imageA,
    int imageB,
    int imageC,
    int imageD,
    int hipA,
    int hipB,
    int hipC,
    int hipD)
{
    std::ostringstream oss;

    oss << imageA << "_"
        << imageB << "_"
        << imageC << "_"
        << imageD << "_"
        << hipA << "_"
        << hipB << "_"
        << hipC << "_"
        << hipD;

    return oss.str();
}

//converts each star to 3d vector of each star in image sensor
std::vector<CameraStar> ConvertToCameraUnitVectors(
    const std::vector<cv::Point2d>& centroids,
    double focalLengthPixels,
    int imageWidth,
    int imageHeight)
{
    std::vector<CameraStar> cameraStars;

    if (imageWidth <= 0 || imageHeight <= 0)
    {
        std::cerr << "ERROR: Invalid image size.\n";
        return cameraStars;
    }

    double cx = static_cast<double>(imageWidth) / 2.0;
    double cy = static_cast<double>(imageHeight) / 2.0;

    for (int i = 0; i < static_cast<int>(centroids.size()); i++)
    {
        double x = centroids[i].x - cx;
        double y = cy - centroids[i].y;
        double z = focalLengthPixels;

        double norm = std::sqrt(x * x + y * y + z * z);

        if (norm == 0.0)
            continue;

        CameraStar star{};
        star.index = i;
        star.name = "Star_" + std::to_string(i + 1);
        star.x = x / norm;
        star.y = y / norm;
        star.z = z / norm;

        cameraStars.push_back(star);
    }

    return cameraStars;
}


//angular distance calculation
double CalculateAngleDeg(
    const CameraStar& a,
    const CameraStar& b)
{
    double dot =
        a.x * b.x +
        a.y * b.y +
        a.z * b.z;

    dot = std::clamp(dot, -1.0, 1.0);

    double angleRad = std::acos(dot);
    double angleDeg = angleRad * 180.0 / PI;

    return angleDeg;
}

//calculates the 6 observed pyramid angles of 4 vectors
static ObservedPyramidAngles CalculateObservedPyramidAngles(
    const std::vector<CameraStar>& cameraStars,
    int a,
    int b,
    int c,
    int d)
{
    ObservedPyramidAngles angles{};

    angles.AB = CalculateAngleDeg(cameraStars[a], cameraStars[b]);
    angles.AC = CalculateAngleDeg(cameraStars[a], cameraStars[c]);
    angles.AD = CalculateAngleDeg(cameraStars[a], cameraStars[d]);
    angles.BC = CalculateAngleDeg(cameraStars[b], cameraStars[c]);
    angles.BD = CalculateAngleDeg(cameraStars[b], cameraStars[d]);
    angles.CD = CalculateAngleDeg(cameraStars[c], cameraStars[d]);

    return angles;
}

static bool IsObservedPyramidUsable(
    const ObservedPyramidAngles& observed,
    double minimumPairAngleDeg = 0.10,
    double minimumLargestAngleDeg = 1.00)
{
    std::vector<double> angles =
    {
        observed.AB,
        observed.AC,
        observed.AD,
        observed.BC,
        observed.BD,
        observed.CD
    };

    double largestAngle =
        0.0;

    for (double angle : angles)
    {
        if (angle < minimumPairAngleDeg)
            return false;

        if (angle > largestAngle)
            largestAngle = angle;
    }

    if (largestAngle < minimumLargestAngleDeg)
        return false;

    return true;
}

//this function is used to find the other star in BC, BD, CD pairs for consistency check
static std::vector<PyramidCandidate> GenerateConsistentPyramidCandidates(
    const std::vector<StarPairAngle>& AB,
    const std::vector<StarPairAngle>& AC,
    const std::vector<StarPairAngle>& AD,
    const std::vector<StarPairAngle>& BC,
    const std::vector<StarPairAngle>& BD,
    const std::vector<StarPairAngle>& CD,
    const ObservedPyramidAngles& observed,
    int imageA,
    int imageB,
    int imageC,
    int imageD)
{
    std::vector<PyramidCandidate> candidates;
    std::unordered_set<std::string> seen;

    auto bcLookup = BuildPairAngleLookup(BC);
    auto bdLookup = BuildPairAngleLookup(BD);
    auto cdLookup = BuildPairAngleLookup(CD);

    for (const auto& ab : AB)
    {
        int firstHip = ab.hipA;
        int secondHip = ab.hipB;

        for (int orientation = 0; orientation < 2; orientation++)
        {
            int hipA = (orientation == 0) ? firstHip : secondHip;
            int hipB = (orientation == 0) ? secondHip : firstHip;

            for (const auto& ac : AC)
            {
                int hipC = 0;

                if (!GetOtherStar(ac, hipA, hipC))
                    continue;

                if (hipC == hipB)
                    continue;

                for (const auto& ad : AD)
                {
                    int hipD = 0;

                    if (!GetOtherStar(ad, hipA, hipD))
                        continue;

                    if (hipD == hipB || hipD == hipC)
                        continue;

                    double catalogBC = 0.0;
                    double catalogBD = 0.0;
                    double catalogCD = 0.0;

                    if (!TryGetPairAngle(bcLookup, hipB, hipC, catalogBC))
                        continue;

                    if (!TryGetPairAngle(bdLookup, hipB, hipD, catalogBD))
                        continue;

                    if (!TryGetPairAngle(cdLookup, hipC, hipD, catalogCD))
                        continue;

                    std::string key =
                        MakeCandidateKey(
                            imageA, imageB, imageC, imageD,
                            hipA, hipB, hipC, hipD
                        );

                    if (seen.find(key) != seen.end())
                        continue;

                    seen.insert(key);

                    PyramidCandidate candidate{};
                    candidate.found = true;

                    candidate.imageA = imageA;
                    candidate.imageB = imageB;
                    candidate.imageC = imageC;
                    candidate.imageD = imageD;

                    candidate.hipA = hipA;
                    candidate.hipB = hipB;
                    candidate.hipC = hipC;
                    candidate.hipD = hipD;

                    double errorAB = std::abs(observed.AB - ab.angle_deg);
                    double errorAC = std::abs(observed.AC - ac.angle_deg);
                    double errorAD = std::abs(observed.AD - ad.angle_deg);
                    double errorBC = std::abs(observed.BC - catalogBC);
                    double errorBD = std::abs(observed.BD - catalogBD);
                    double errorCD = std::abs(observed.CD - catalogCD);

                    candidate.totalAngleErrorDeg =
                        errorAB + errorAC + errorAD +
                        errorBC + errorBD + errorCD;

                    candidate.averageAngleErrorDeg =
                        candidate.totalAngleErrorDeg / 6.0;

                    candidates.push_back(candidate);
                }
            }
        }
    }

    return candidates;
}

//*************OPTIMIZATION NEEDED HERE*************
//function to searche all possible 4-star image groups 
//this is brute-force search and can be optimized with more advanced techniques
static std::vector<PyramidCandidate> FindAllPyramidCandidates(
    const std::vector<CameraStar>& cameraStars,
    const std::vector<StarPairAngle>& pairDatabase,
    double toleranceDeg)
{
    std::vector<PyramidCandidate> allCandidates;

    int n = static_cast<int>(cameraStars.size());

    if (n < 4)
    {
        std::cerr << "ERROR: At least 4 detected stars are required.\n";
        return allCandidates;
    }

    //this here nested loop is brute-force search for all possible 4-star combinations in the detected image stars
    for (int a = 0; a < n - 3; a++)
    {
        for (int b = a + 1; b < n - 2; b++)
        {
            for (int c = b + 1; c < n - 1; c++)
            {
                for (int d = c + 1; d < n; d++)
                {
                    ObservedPyramidAngles observed =
                        CalculateObservedPyramidAngles(
                            cameraStars,
                            a, b, c, d
                        );

                    if (!IsObservedPyramidUsable(observed))
                        continue;

                    std::vector<StarPairAngle> AB =
                        SearchPairDatabase(pairDatabase, observed.AB, toleranceDeg);

                    std::vector<StarPairAngle> AC =
                        SearchPairDatabase(pairDatabase, observed.AC, toleranceDeg);

                    std::vector<StarPairAngle> AD =
                        SearchPairDatabase(pairDatabase, observed.AD, toleranceDeg);

                    std::vector<StarPairAngle> BC =
                        SearchPairDatabase(pairDatabase, observed.BC, toleranceDeg);

                    std::vector<StarPairAngle> BD =
                        SearchPairDatabase(pairDatabase, observed.BD, toleranceDeg);

                    std::vector<StarPairAngle> CD =
                        SearchPairDatabase(pairDatabase, observed.CD, toleranceDeg);

                    if (AB.empty() || AC.empty() || AD.empty() ||
                        BC.empty() || BD.empty() || CD.empty())
                    {
                        continue;
                    }

                    std::vector<PyramidCandidate> candidates =
                        GenerateConsistentPyramidCandidates(
                            AB, AC, AD,
                            BC, BD, CD,
                            observed,
                            a, b, c, d
                        );

                    allCandidates.insert(
                        allCandidates.end(),
                        candidates.begin(),
                        candidates.end()
                    );
                }
            }
        }
    }

    return allCandidates;
}

//voting function
static void ApplyVotingScores(
    std::vector<PyramidCandidate>& candidates)
{
    std::unordered_map<int, std::unordered_map<int, int>> votes;

    for (const auto& c : candidates)
    {
        votes[c.imageA][c.hipA]++;
        votes[c.imageB][c.hipB]++;
        votes[c.imageC][c.hipC]++;
        votes[c.imageD][c.hipD]++;
    }

    for (auto& c : candidates)
    {
        c.voteScore =
            votes[c.imageA][c.hipA] +
            votes[c.imageB][c.hipB] +
            votes[c.imageC][c.hipC] +
            votes[c.imageD][c.hipD];
    }
}

//if there is 2 pyramids output in voting 
//this function will select the one with the highest vote score and lowest average angle error as the best candidate.
static IdentifiedPyramid SelectBestPyramidCandidate(
    std::vector<PyramidCandidate>& candidates)
{
    IdentifiedPyramid result{};

    if (candidates.empty())
        return result;

    ApplyVotingScores(candidates);

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const PyramidCandidate& a, const PyramidCandidate& b)
        {
            if (a.voteScore != b.voteScore)
                return a.voteScore > b.voteScore;

            return a.averageAngleErrorDeg < b.averageAngleErrorDeg;
        }
    );

    std::cout << "\nTotal pyramid candidates = "
        << candidates.size() << std::endl;

    std::cout << "\nTop pyramid candidates:\n";

    int printCount =
        std::min(10, static_cast<int>(candidates.size()));

    for (int i = 0; i < printCount; i++)
    {
        const auto& c = candidates[i];

        std::cout << i + 1 << ") "
            << "Image Stars ["
            << c.imageA << ", "
            << c.imageB << ", "
            << c.imageC << ", "
            << c.imageD << "] -> HIP ["
            << c.hipA << ", "
            << c.hipB << ", "
            << c.hipC << ", "
            << c.hipD << "]"
            << " | voteScore = " << c.voteScore
            << " | avgAngleError = "
            << c.averageAngleErrorDeg
            << " deg\n";
    }

    const PyramidCandidate& best = candidates.front();

    result.found = true;

    result.imageA = best.imageA;
    result.imageB = best.imageB;
    result.imageC = best.imageC;
    result.imageD = best.imageD;

    result.hipA = best.hipA;
    result.hipB = best.hipB;
    result.hipC = best.hipC;
    result.hipD = best.hipD;

    result.voteScore = best.voteScore;
    result.averageAngleErrorDeg = best.averageAngleErrorDeg;

    return result;
}

//in real star tracker pipeline this function is not appropriate 
// the input to attitude estimation should be not only one pair of matched pyramid but all the matched pyramid in the image.
static std::vector<MatchedStar> BuildMatchedStarsForAttitude(
    const IdentifiedPyramid& result,
    const std::vector<CatalogStar>& catalogStars)
{
    std::vector<MatchedStar> matchedStars;

    if (!result.found)
    {
        std::cerr << "ERROR: Cannot build matched stars. No pyramid found.\n";
        return matchedStars;
    }

    std::unordered_map<int, CatalogStar> catalogMap;
    catalogMap.reserve(catalogStars.size());

    for (const auto& star : catalogStars)
    {
        catalogMap[star.hip] = star;
    }

    const int imageIndices[4] =
    {
        result.imageA,
        result.imageB,
        result.imageC,
        result.imageD
    };

    const int hipIds[4] =
    {
        result.hipA,
        result.hipB,
        result.hipC,
        result.hipD
    };

    for (int i = 0; i < 4; i++)
    {
        int imageIndex = imageIndices[i];
        int hip = hipIds[i];

        if (imageIndex < 0 ||
            imageIndex >= static_cast<int>(result.cameraStars.size()))
        {
            std::cerr << "ERROR: Invalid camera star index: "
                << imageIndex << "\n";
            continue;
        }

        auto catalogIt = catalogMap.find(hip);

        if (catalogIt == catalogMap.end())
        {
            std::cerr << "ERROR: HIP "
                << hip
                << " not found in catalogStars.\n";
            continue;
        }

        const CameraStar& cameraStar =
            result.cameraStars[imageIndex];

        const CatalogStar& catalogStar =
            catalogIt->second;

        MatchedStar matched{};

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

void SaveLabeledIdentifiedPyramidImage(
    const std::string& inputImagePath,
    const std::string& outputImagePath,
    const IdentifiedPyramid& result)
{
    if (!result.found)
    {
        std::cout << "No identified pyramid to label.\n";
        return;
    }

    if (result.centroids.empty())
    {
        std::cerr << "ERROR: No centroid data available.\n";
        return;
    }

    cv::Mat image = cv::imread(inputImagePath, cv::IMREAD_COLOR);

    if (image.empty())
    {
        std::cerr << "ERROR: Cannot load image: "
            << inputImagePath << std::endl;
        return;
    }

    struct LabelItem
    {
        int imageIndex;
        int hipId;
    };

    const std::array<LabelItem, 4> labels =
    { {
        { result.imageA, result.hipA },
        { result.imageB, result.hipB },
        { result.imageC, result.hipC },
        { result.imageD, result.hipD }
    } };

    constexpr int Radius = 6;
    constexpr int Thickness = 1;
    constexpr double FontScale = 0.6;
    constexpr int FontFace = cv::FONT_HERSHEY_SIMPLEX;

    const cv::Scalar CircleColor(0, 255, 255);   // Yellow
    const cv::Scalar TextColor(0, 255, 255);

    const int centroidCount = static_cast<int>(result.centroids.size());

    for (const auto& label : labels)
    {
        if (label.imageIndex < 0 || label.imageIndex >= centroidCount)
            continue;

        const cv::Point2d& centroid = result.centroids[label.imageIndex];

        const cv::Point center(
            cvRound(centroid.x),
            cvRound(centroid.y));

        // Draw star marker
        cv::circle(
            image,
            center,
            Radius,
            CircleColor,
            cv::FILLED);

        // Draw HIP label
        cv::putText(
            image,
            "HIP " + std::to_string(label.hipId),
            cv::Point(center.x + 10, center.y - 10),
            FontFace,
            FontScale,
            TextColor,
            Thickness,
            cv::LINE_AA);
    }

    if (!cv::imwrite(outputImagePath, image))
    {
        std::cerr << "ERROR: Cannot save labeled image: "
            << outputImagePath << std::endl;
        return;
    }

    std::cout << "Labeled identified pyramid image saved to: "
        << outputImagePath << std::endl;
}

IdentifiedPyramid IdentifyPyramidFromBinaryImage(
    const std::string& binaryImagePath,
    const std::vector<StarPairAngle>& pairDatabase,
    const std::vector<CatalogStar>& catalogStars,
    double focalLengthPixels,
    double toleranceDeg,
    int maxDetectedStars)
{
    auto [centroids, imageWidth, imageHeight] =
        FindContoursCentroids(
            binaryImagePath,
            maxDetectedStars,
            1.0
        );

    if (centroids.empty())
    {
        std::cerr << "ERROR: No stars detected in image.\n";
        return {};
    }

    std::vector<CameraStar> cameraStars =
        ConvertToCameraUnitVectors(
            centroids,
            focalLengthPixels,
            imageWidth,
            imageHeight
        );

    if (cameraStars.size() < 4)
    {
        std::cerr << "ERROR: Fewer than 4 valid stars in captured image.\n";
        return {};
    }

    if (pairDatabase.empty())
    {
        std::cerr << "ERROR: Pair database is empty.\n";
        return {};
    }

    std::vector<PyramidCandidate> candidates =
        FindAllPyramidCandidates(
            cameraStars,
            pairDatabase,
            toleranceDeg
        );

    IdentifiedPyramid result =
        SelectBestPyramidCandidate(candidates);

    result.centroids = std::move(centroids);
    result.cameraStars = std::move(cameraStars);
    result.imageWidth = imageWidth;
    result.imageHeight = imageHeight;

    result.matchedStars =
        BuildMatchedStarsForAttitude(result, catalogStars);

    return result;
}

StarIdentificationDebugData IdentifyAllPyramidCandidatesForDebugging(
    const std::string& binaryImagePath,
    const std::vector<StarPairAngle>& pairDatabase,
    double focalLengthPixels,
    double toleranceDeg,
    int maxDetectedStars)
{

    StarIdentificationDebugData output;

    auto [centroids, imageWidth, imageHeight] =
        FindContoursCentroids(binaryImagePath, maxDetectedStars, 1.0);

    if (centroids.empty())
    {
        std::cerr << "ERROR: No stars detected.\n";
        return output;
    }

    std::vector<CameraStar> cameraStars =
        ConvertToCameraUnitVectors(
            centroids,
            focalLengthPixels,
            imageWidth,
            imageHeight);

    if (pairDatabase.empty())
    {
        std::cerr << "ERROR: Pair database is empty.\n";
        return output;
    }

    std::vector<PyramidCandidate> candidates =
        FindAllPyramidCandidates(
            cameraStars,
            pairDatabase,
            toleranceDeg);

    if (candidates.empty())
    {
        std::cerr << "ERROR: No pyramid candidates found.\n";
        return output;
    }

    ApplyVotingScores(candidates);

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const PyramidCandidate& a, const PyramidCandidate& b)
        {
            if (a.voteScore != b.voteScore)
                return a.voteScore > b.voteScore;

            return a.averageAngleErrorDeg < b.averageAngleErrorDeg;
        });

    output.success = true;
    output.centroids = centroids;
    output.cameraStars = cameraStars;
    output.candidates = candidates;
    output.imageWidth = imageWidth;
    output.imageHeight = imageHeight;

    return output;
}

//void SaveLabeledIdentifiedPyramidImage(
//    const std::string& inputImagePath,
//    const std::string& outputImagePath,
//    const IdentifiedPyramid& result)
//{
//    if (!result.found)
//    {
//        std::cout << "No identified pyramid to label.\n";
//        return;
//    }
//
//    cv::Mat image =
//        cv::imread(inputImagePath, cv::IMREAD_COLOR);
//
//    if (image.empty())
//    {
//        std::cerr << "ERROR: Cannot load image for labeling: "
//            << inputImagePath << std::endl;
//        return;
//    }
//
//    if (result.centroids.empty())
//    {
//        std::cerr << "ERROR: No centroid data available for labeling.\n";
//        return;
//    }
//
//    struct LabelItem
//    {
//        int imageIndex;
//        int hipId;
//    };
//
//    std::vector<LabelItem> labels =
//    {
//        { result.imageA, result.hipA },
//        { result.imageB, result.hipB },
//        { result.imageC, result.hipC },
//        { result.imageD, result.hipD }
//    };
//
//    int fontFace = cv::FONT_HERSHEY_SIMPLEX;
//    double fontScale = 0.6;
//    int thickness = 1;
//
//    for (const auto& label : labels)
//    {
//        if (label.imageIndex < 0 ||
//            label.imageIndex >= static_cast<int>(result.centroids.size()))
//        {
//            continue;
//        }
//
//        cv::Point center(
//            static_cast<int>(std::round(result.centroids[label.imageIndex].x)),
//            static_cast<int>(std::round(result.centroids[label.imageIndex].y))
//        );
//
//        cv::circle(
//            image,
//            center,
//            6,
//            cv::Scalar(0, 255, 255),
//            -1
//        );
//
//        std::string text =
//            "HIP " + std::to_string(label.hipId);
//
//        cv::Point textPosition(
//            center.x + 10,
//            center.y - 10
//        );
//
//        cv::putText(
//            image,
//            text,
//            textPosition,
//            fontFace,
//            fontScale,
//            cv::Scalar(0, 255, 255),
//            thickness
//        );
//    }
//
//    bool ok =
//        cv::imwrite(outputImagePath, image);
//
//    if (!ok)
//    {
//        std::cerr << "ERROR: Cannot save labeled image: "
//            << outputImagePath << std::endl;
//        return;
//    }
//
//    std::cout << "Labeled identified pyramid image saved to: "
//        << outputImagePath << std::endl;
//}

//print the selected best matching pyramid with the corresponding image stars, HIPs, vote score, and average angle error.
//void PrintIdentificationResult(
//    const IdentifiedPyramid& result)
//{
//    if (!result.found)
//    {
//        std::cout << "No consistent catalog pyramid found.\n";
//        return;
//    }
//
//    std::cout << "\nBest pyramid selected:\n";
//
//    std::cout << "Image Star " << result.imageA + 1
//        << "  ->  HIP " << result.hipA << std::endl;
//
//    std::cout << "Image Star " << result.imageB + 1
//        << "  ->  HIP " << result.hipB << std::endl;
//
//    std::cout << "Image Star " << result.imageC + 1
//        << "  ->  HIP " << result.hipC << std::endl;
//
//    std::cout << "Image Star " << result.imageD + 1
//        << "  ->  HIP " << result.hipD << std::endl;
//
//    std::cout << "Vote score          = "
//        << result.voteScore << std::endl;
//
//    std::cout << "Average angle error = "
//        << result.averageAngleErrorDeg
//        << " deg\n";
//}