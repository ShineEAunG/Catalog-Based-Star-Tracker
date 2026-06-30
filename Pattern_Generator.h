#pragma once

#include <vector>
#include <string>
#include "Parameters.h"

std::vector<StarPairAngle> GeneratePairDatabase(
    const std::vector<CatalogStar>& stars
);

void SavePairDatabaseBinary(
    const std::vector<StarPairAngle>& pairs,
    const std::string& outputFile
);

std::vector<StarPairAngle> LoadPairDatabaseBinary(
    const std::string& inputFile
);

std::vector<StarPairAngle> SearchPairDatabase(
    const std::vector<StarPairAngle>& pairs,
    double targetAngle,
    double tolerance
);