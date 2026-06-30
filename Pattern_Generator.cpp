//generate a database of star pairs based on their angular separation. 
//includes functions to calculate the angle between two stars, 
//generate a pair database, save and load the database in binary format, 
//and search for star pairs within a specified angular tolerance.

#include "Pattern_Generator.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

static double CalculateCatalogAngle(
    const CatalogStar& a,
    const CatalogStar& b)
{
    double dot =
        a.x * b.x +
        a.y * b.y +
        a.z * b.z;

    dot = std::max(-1.0, std::min(1.0, dot));

    double angle_rad = std::acos(dot);
    return angle_rad * 180.0 / 3.14159265358979323846;
}

std::vector<StarPairAngle> GeneratePairDatabase(
    const std::vector<CatalogStar>& stars)
{
    std::vector<StarPairAngle> pairs;

    int n = static_cast<int>(stars.size());

    long long totalPairs =
        static_cast<long long>(n) * (n - 1) / 2;

    std::cout << "Generating pair database...\n";
    std::cout << "Stars = " << n << std::endl;
    std::cout << "Expected pairs = " << totalPairs << std::endl;

    pairs.reserve(totalPairs);

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double angle = CalculateCatalogAngle(stars[i], stars[j]);

            StarPairAngle pair{};
            pair.hipA = stars[i].hip;
            pair.hipB = stars[j].hip;
            pair.angle_deg = angle;

            pairs.push_back(pair);
        }
    }

    std::sort(
        pairs.begin(),
        pairs.end(),
        [](const StarPairAngle& a, const StarPairAngle& b)
        {
            return a.angle_deg < b.angle_deg;
        }
    );

    std::cout << "Generated pairs = "
        << pairs.size()
        << std::endl;

    return pairs;
}

void SavePairDatabaseBinary(
    const std::vector<StarPairAngle>& pairs,
    const std::string& outputFile)
{
    std::ofstream out(outputFile, std::ios::binary);

    if (!out)
    {
        std::cout << "ERROR: Cannot create pair database file\n";
        return;
    }

    size_t count = pairs.size();

    out.write(
        reinterpret_cast<const char*>(&count),
        sizeof(size_t)
    );

    if (!pairs.empty())
    {
        out.write(
            reinterpret_cast<const char*>(pairs.data()),
            pairs.size() * sizeof(StarPairAngle)
        );
    }

    if (!out.good())
    {
        std::cerr << "ERROR: Failed while writing pair database.\n";
        return;
    }

    out.close();

    std::cout << "Pair database saved: "
        << outputFile
        << std::endl;
}

std::vector<StarPairAngle> LoadPairDatabaseBinary(
    const std::string& inputFile)
{
    std::vector<StarPairAngle> pairs;

    std::ifstream in(inputFile, std::ios::binary);

    if (!in)
    {
        std::cout << "ERROR: Cannot open pair database binary file\n";
        return pairs;
    }

    size_t count = 0;
    in.read(reinterpret_cast<char*>(&count), sizeof(size_t));

    if (!in.good())
    {
        std::cerr << "ERROR: Cannot read pair database count.\n";
        return {};
    }
    pairs.resize(count);

    in.read(reinterpret_cast<char*>(pairs.data()),
        count * sizeof(StarPairAngle));

    in.close();

    std::cout << "Pair database loaded: " << inputFile << std::endl;
    std::cout << "Loaded pairs = " << pairs.size() << std::endl;

    return pairs;
}

std::vector<StarPairAngle> SearchPairDatabase(
    const std::vector<StarPairAngle>& pairs,
    double targetAngle,
    double tolerance)
{
    std::vector<StarPairAngle> candidates;

    double lowerAngle = targetAngle - tolerance;
    double upperAngle = targetAngle + tolerance;

    auto lower = std::lower_bound(
        pairs.begin(),
        pairs.end(),
        lowerAngle,
        [](const StarPairAngle& pair, double value)
        {
            return pair.angle_deg < value;
        }
    );

    auto upper = std::upper_bound(
        pairs.begin(),
        pairs.end(),
        upperAngle,
        [](double value, const StarPairAngle& pair)
        {
            return value < pair.angle_deg;
        }
    );

    candidates.assign(lower, upper);
    return candidates;
}