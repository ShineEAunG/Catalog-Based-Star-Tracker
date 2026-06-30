#pragma once

#include <string>
#include <vector>
#include <set>
#include "Parameters.h"

std::vector<TruthStar> LoadTruthMetadataCSV(
    const std::string& csvPath);

ValidationSummary ValidateSyntheticResult(
    const ReprojectionResult& selectedResult,
    const std::vector<TruthStar>& truthStars,
    const AttitudeResult& svdResult,
    const AttitudeResult& questResult,
    bool hasTruthQuaternion,
    Quaternion qTruth);

void PrintValidationSummary(
    const ValidationSummary& summary);
