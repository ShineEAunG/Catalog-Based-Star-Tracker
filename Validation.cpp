#include "Validation.h"
#include "AttitudeEstimation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

static std::string TrimLocal(const std::string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");

    if (start == std::string::npos)
        return "";

    size_t end = s.find_last_not_of(" \t\r\n");

    return s.substr(start, end - start + 1);
}

std::vector<TruthStar> LoadTruthMetadataCSV(
    const std::string& csvPath)
{
    std::vector<TruthStar> truthStars;

    std::ifstream file(csvPath);

    if (!file.is_open())
    {
        std::cerr << "WARNING: Cannot open truth metadata: "
            << csvPath << "\n";

        return truthStars;
    }

    std::string line;

    // Skip header 
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (TrimLocal(line).empty())
            continue;

        std::stringstream ss(line);
        std::string token;

        TruthStar star;

        try
        {
            std::getline(ss, token, ',');
            star.hip = std::stoi(TrimLocal(token));

            std::getline(ss, token, ',');
            star.u = std::stod(TrimLocal(token));

            std::getline(ss, token, ',');
            star.v = std::stod(TrimLocal(token));

            if (std::getline(ss, token, ','))
            {
                star.brightness = std::stod(TrimLocal(token));
            }

            truthStars.push_back(star);
        }
        catch (...)
        {
            continue;
        }
    }

    return truthStars;
}

static std::set<int> BuildTruthHipSet(
    const std::vector<TruthStar>& truthStars)
{
    std::set<int> hipSet;

    for (const auto& star : truthStars)
    {
        hipSet.insert(star.hip);
    }

    return hipSet;
}

static std::vector<int> GetCandidateHips(
    const PyramidCandidate& candidate)
{
    std::vector<int> hips;

    hips.push_back(candidate.hipA);
    hips.push_back(candidate.hipB);
    hips.push_back(candidate.hipC);
    hips.push_back(candidate.hipD);

    return hips;
}

ValidationSummary ValidateSyntheticResult(
    const ReprojectionResult& selectedResult,
    const std::vector<TruthStar>& truthStars,
    const AttitudeResult& svdResult,
    const AttitudeResult& questResult,
    bool hasTruthQuaternion,
    Quaternion qTruth)
{
    ValidationSummary summary;

    summary.truthStars =
        static_cast<int>(truthStars.size());

    if (!selectedResult.success)
    {
        summary.pass = false;
        return summary;
    }

    std::set<int> truthHipSet =
        BuildTruthHipSet(truthStars);

    std::vector<int> candidateHips =
        GetCandidateHips(selectedResult.candidate);

    summary.selectedPyramidStars =
        static_cast<int>(candidateHips.size());

    int correct = 0;

    for (int hip : candidateHips)
    {
        if (truthHipSet.find(hip) != truthHipSet.end())
        {
            correct++;
        }
    }

    summary.correctPyramidStars = correct;

    if (svdResult.success)
    {
        summary.svdResidualArcsec =
            svdResult.rmsResidualDeg * 3600.0;
    }

    if (questResult.success)
    {
        summary.questResidualArcsec =
            questResult.rmsResidualDeg * 3600.0;
    }

    if (hasTruthQuaternion)
    {
        Quaternion svdQ = svdResult.q;
        Quaternion questQ = questResult.q;

        summary.svdAttitudeErrorDeg =
            QuaternionAngleErrorDeg(
                svdQ,
                qTruth);

        summary.questAttitudeErrorDeg =
            QuaternionAngleErrorDeg(
                questQ,
                qTruth);
    }

    bool pyramidPass =
        summary.correctPyramidStars ==
        summary.selectedPyramidStars;

    bool reprojectionPass =
        selectedResult.inlierCount >=
        summary.selectedPyramidStars;

    bool residualPass =
        selectedResult.attitude.rmsResidualDeg < 0.05;

    summary.pass =
        pyramidPass &&
        reprojectionPass &&
        residualPass;

    return summary;
}

void PrintValidationSummary(
    const ValidationSummary& summary)
{
    std::cout << "\n*********** VALIDATION SUMMARY ***********\n";

    std::cout << "Truth stars in metadata       = "
        << summary.truthStars << "\n";

    std::cout << "Selected pyramid stars        = "
        << summary.selectedPyramidStars << "\n";

    std::cout << "Correct selected HIP stars    = "
        << summary.correctPyramidStars
        << " / "
        << summary.selectedPyramidStars
        << "\n";

    std::cout << "SVD residual                  = "
        << summary.svdResidualArcsec
        << " arcsec\n";

    std::cout << "QUEST residual                = "
        << summary.questResidualArcsec
        << " arcsec\n";

    if (summary.svdAttitudeErrorDeg >= 0.0)
    {
        std::cout << "SVD attitude error            = "
            << summary.svdAttitudeErrorDeg
            << " deg\n";
    }

    if (summary.questAttitudeErrorDeg >= 0.0)
    {
        std::cout << "QUEST attitude error          = "
            << summary.questAttitudeErrorDeg
            << " deg\n";
    }

    std::cout << "Result                        = "
        << (summary.pass ? "PASS" : "FAIL")
        << "\n";
}