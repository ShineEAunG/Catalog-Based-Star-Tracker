#pragma once

#include <string>
#include <vector>
#include "Parameters.h"

std::vector<CatalogStar> LoadHipparcosCatalog(
    const std::string& fileName,
    double magnitudeLimit
);

std::string SaveCatalogBinary(
    const std::vector<CatalogStar>& stars,
    const std::string& outputFile
);

std::vector<CatalogStar> LoadCatalogBinary(
    const std::string& inputFile
);
