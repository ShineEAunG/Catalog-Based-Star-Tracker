//session to convert HIPPARCOS catalog(.bat) to .bin 

#include "CatalogGenerator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <vector>
#include <string>

#ifndef CV_PI
#define CV_PI 3.14159265358979323846
#endif

static std::string Trim(const std::string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");

    if (start == std::string::npos)
        return "";

    return s.substr(start, end - start + 1);
}

static bool TryParseDouble(const std::string& text, double& value)
{
    std::string t = Trim(text);

    if (t.empty())
        return false;

    try
    {
        value = std::stod(t);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

static bool TryParseInt(const std::string& text, int& value)
{
    std::string t = Trim(text);

    if (t.empty())
        return false;

    try
    {
        value = std::stoi(t);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

std::vector<CatalogStar> LoadHipparcosCatalog(
    const std::string& fileName,
    double magnitudeLimit)
{
    std::vector<CatalogStar> stars;

    std::ifstream file(fileName);

    if (!file.is_open())
    {
        std::cout << "ERROR: Cannot open Hipparcos file: "
            << fileName << std::endl;
        return stars;
    }

    std::string line;
    int totalLines = 0;
    int skippedLines = 0;

    while (std::getline(file, line))
    {
        totalLines++;

        if (line.length() < 76)
        {
            skippedLines++;
            continue;
        }

        int hip;
        double vmag;
        double ra_deg;
        double dec_deg;

        // ReadMe byte positions:
        // HIP   9-14   -> substr(8, 6)
        // Vmag  42-46  -> substr(41, 5)
        // RAdeg 52-63  -> substr(51, 12)
        // DEdeg 65-76  -> substr(64, 12)

        bool okHip = TryParseInt(line.substr(8, 6), hip);
        bool okMag = TryParseDouble(line.substr(41, 5), vmag);
        bool okRa = TryParseDouble(line.substr(51, 12), ra_deg);
        bool okDec = TryParseDouble(line.substr(64, 12), dec_deg);

        if (!okHip || !okMag || !okRa || !okDec)
        {
            skippedLines++;
            continue;
        }

        if (vmag > magnitudeLimit)
            continue;

        double ra_rad = ra_deg * CV_PI / 180.0;
        double dec_rad = dec_deg * CV_PI / 180.0;

        CatalogStar star{};
        star.hip = hip;
        star.ra_deg = ra_deg;
        star.dec_deg = dec_deg;
        star.magnitude = vmag;

        star.x = std::cos(dec_rad) * std::cos(ra_rad);
        star.y = std::cos(dec_rad) * std::sin(ra_rad);
        star.z = std::sin(dec_rad);

        stars.push_back(star);
    }

    file.close();

    std::cout << "Hipparcos total lines  = " << totalLines << std::endl;
    std::cout << "Skipped invalid lines  = " << skippedLines << std::endl;
    std::cout << "Loaded filtered stars  = " << stars.size() << std::endl;

    return stars;
}


std::string SaveCatalogBinary(
    const std::vector<CatalogStar>& stars,
    const std::string& outputFolder)
{
    std::string outputFile = "hip_catalog_mag6.bin";

    if (!std::filesystem::exists(outputFolder))
    {
        std::filesystem::create_directories(outputFolder);
    }

    // combine paths using C++17 filesystem
    outputFile = (std::filesystem::path(outputFolder) / outputFile).string();

    std::ofstream out(outputFile, std::ios::binary);
    if (!out)
    {
        std::cerr << "ERROR: Cannot create catalog binary file: " << outputFile << "\n";
        return "";
    }

    size_t count = stars.size();

    out.write(reinterpret_cast<const char*>(&count), sizeof(size_t));

    if (!stars.empty())
    {
        out.write(reinterpret_cast<const char*>(stars.data()),
            stars.size() * sizeof(CatalogStar));
    }

    out.close();
    return outputFile;
}


std::vector<CatalogStar> LoadCatalogBinary(
    const std::string& inputFile)
{
    std::vector<CatalogStar> stars;

    std::ifstream in(inputFile, std::ios::binary);

    if (!in)
    {
        std::cout << "ERROR: Cannot open catalog binary file\n";
        return stars;
    }

    size_t count = 0;
    in.read(reinterpret_cast<char*>(&count), sizeof(size_t));

    stars.resize(count);

    in.read(reinterpret_cast<char*>(stars.data()),
        count * sizeof(CatalogStar));

    in.close();

    std::cout << "Catalog binary loaded: " << inputFile << std::endl;
    std::cout << "Loaded stars = " << stars.size() << std::endl;

    return stars;
}