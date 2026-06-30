#include "TestImageGeneration.h"

//****************Star Field Data Section****************
std::vector<Star> GetPolarisCostellation()
{
    std::vector<Star> PolarisConstellation;
    Star polaris = { "Polaris", 46.8004, 89.2642, 2.09 };
    Star yildun = { "Yildun", 260.93125, 86.56760, 4.34 };
    Star epsilonUrsaMinor = { "Epsilon Ursa Minor", 250.8225, 81.99090, 4.34 };

    Star zetaUrsaMinor = { "Zeta Ursa Minor", 235.7775, 77.71290, 4.28 };
    Star konchab = { "Konchab", 222.6504, 74.04700, 2.20 };
    Star pherkad = { "Pherkad", 230.16375, 71.74080, 3.03 };
    Star etaUrsaMinor = { "Eta Ursa Minor", 244.1775, 75.69450, 5.03 };//286.0375

    PolarisConstellation.push_back(polaris);
    PolarisConstellation.push_back(yildun);
    PolarisConstellation.push_back(epsilonUrsaMinor);

    PolarisConstellation.push_back(zetaUrsaMinor);
    PolarisConstellation.push_back(konchab);
    PolarisConstellation.push_back(pherkad);
    PolarisConstellation.push_back(etaUrsaMinor);

    std::cout << "Returning Polaris MetaData\n";
    return PolarisConstellation;
}
std::vector<Star> GetCassiopeiaConstellation()
{
    return
    {
        {"Schedar",   10.1260, 56.5373, 2.24},
        {"Caph",       2.2937, 59.1498, 2.28},
        {"Gamma Cas", 14.1771, 60.7167, 2.15},
        {"Ruchbah",   21.4539, 60.2353, 2.68},
        {"Segin",     28.5988, 63.6701, 3.35}
    };
}
std::vector<Star> GetOrionConstellation()
{
    return
    {
        {"Betelgeuse",  88.7929,  7.4071, 0.42},
        {"Bellatrix",   81.2828,  6.3497, 1.64},
        {"Alnitak",     85.1897, -1.9426, 1.74},
        {"Alnilam",     84.0533, -1.2019, 1.69},
        {"Mintaka",     83.0017, -0.2991, 2.23},
        {"Saiph",       86.9391, -9.6696, 2.07},
        {"Rigel",       78.6345, -8.2016, 0.13}
    };
}
std::vector<Star> GetBigDipper()
{
    return
    {
        {"Dubhe",    165.4603, 61.7510, 1.79},
        {"Merak",    165.9320, 56.3824, 2.37},
        {"Phecda",   177.2649, 53.6948, 2.41},
        {"Megrez",   183.8563, 57.0326, 3.32},
        {"Alioth",   193.5073, 55.9598, 1.76},
        {"Mizar",    200.9814, 54.9254, 2.23},
        {"Alkaid",   206.8852, 49.3133, 1.85}
    };
}
std::vector<Star> GetCygnus()
{
    return
    {
        {"Deneb",   310.3579, 45.2803, 1.25},
        {"Sadr",    305.5571, 40.2567, 2.23},
        {"Albireo", 292.6803, 27.9597, 3.18},
        {"Gienah",  311.5528, 33.9703, 2.46},
        {"Delta Cyg",296.2437,45.1308, 2.87}
    };
}

//****************Vector Conversion Section****************
//known attitude quaternion q and a vector v,
//this function rotates the vector v by the quaternion q 
std::vector<double> RotateVector(const Quaternion& q, 
    const std::vector<double>& v)
{
    double qw = q.w, qx = q.x, qy = q.y, qz = q.z;

    // quaternion-vector multiplication (optimized form)
    double ix = qw * v[0] + qy * v[2] - qz * v[1];
    double iy = qw * v[1] + qz * v[0] - qx * v[2];
    double iz = qw * v[2] + qx * v[1] - qy * v[0];
    double iw = -qx * v[0] - qy * v[1] - qz * v[2];

    std::vector<double> res(3);

    res[0] = ix * qw + iw * -qx + iy * -qz - iz * -qy;
    res[1] = iy * qw + iw * -qy + iz * -qx - ix * -qz;
    res[2] = iz * qw + iw * -qz + ix * -qy - iy * -qx;

    return res;
}

std::vector<CameraTruthStar> ConvertCatalogIntoCameraFrame(
    const std::vector<CatalogStar>& catalogStars,
    const Quaternion& q_truth,
    double maxBrightness,
    int width,
    int height,
    double focalLengthPixels,
    double magnitudeLimit)
{
    std::vector<CameraTruthStar> cameraStars;

    double limitX =
        (static_cast<double>(width) / 2.0) / focalLengthPixels;

    double limitY =
        (static_cast<double>(height) / 2.0) / focalLengthPixels;

    for (const auto& star : catalogStars)
    {
        if (star.magnitude > magnitudeLimit)
            continue;
        std::vector<double> inertialVector =
        {
            star.x,
            star.y,
            star.z
        };

        std::vector<double> cameraVector =
            RotateVector(q_truth, inertialVector);

        double x = cameraVector[0];
        double y = cameraVector[1];
        double z = cameraVector[2];

        // remove stars behind the camera (z <= 0)
        if (z <= 0.0)
            continue;

        // Rectangular camera FOV test
        if (std::fabs(x / z) > limitX)
            continue;

        if (std::fabs(y / z) > limitY)
            continue;

        double brightness =
            maxBrightness / std::pow(2.512, star.magnitude);

        cameraStars.push_back(
            {
                star.hip,
                x,
                y,
                z,
                brightness
            });
    }

    return cameraStars;
}

std::vector<PixelCoordinate> ConvertIntoPixelCoordinates(
    const std::vector<CameraTruthStar>& cameraStars,
    int width,
    int height,
    double focalLengthPixels)
{
    std::vector<PixelCoordinate> pixelCoordinates;

    double centerX =
        static_cast<double>(width) / 2.0;

    double centerY =
        static_cast<double>(height) / 2.0;

    for (const auto& star : cameraStars)
    {
        if (star.z <= 0.0)
            continue;

        double imageX =
            focalLengthPixels * (star.x / star.z);

        double imageY =
            focalLengthPixels * (star.y / star.z);

        // x = centroid.x - cx
        // y = cy - centroid.y
        double pixelU =
            centerX + imageX;

        double pixelV =
            centerY - imageY;

        if (pixelU < 0.0 || pixelU >= width)
            continue;

        if (pixelV < 0.0 || pixelV >= height)
            continue;

        pixelCoordinates.push_back(
            {
                star.hip,
                pixelU,
                pixelV,
                star.brightness
            });
    }

    return pixelCoordinates;
}


//****************File Explorer Section****************
//template<typename T>//for dynamic type we use template
void LogDataInConsole(
    const std::vector<PixelCoordinate>& coordinates)
{
    for (const auto& coordinate : coordinates)
    {
        std::cout << "HIP: "
            << coordinate.hip << std::endl;

        std::cout << "u: "
            << coordinate.u << std::endl;

        std::cout << "v: "
            << coordinate.v << std::endl;

        std::cout << "Brightness: "
            << coordinate.brightness << std::endl;

        std::cout << "------------------------------------------------------"
            << std::endl;
    }
}

void WriteCoordinatesToFile(
    const std::vector<PixelCoordinate>& stars,
    const std::string& fileName,
    const std::string& folder)
{
    std::string filePath =
        folder + fileName;

    std::ofstream file(filePath);

    file << "hip,u,v,brightness\n";

    for (const auto& star : stars)
    {
        file << star.hip << ","
            << star.u << ","
            << star.v << ","
            << static_cast<int>(star.brightness)
            << std::endl;
    }

    file.close();
}

//****************Star Field Image Section****************
void Generate3ColorsConstellation(
    const std::vector<PixelCoordinate>& pixelCoordinates,
    int width,
    int height,
    const std::string& imageName,
    const std::string& folder)
{
    if (!std::filesystem::exists(folder))
    {
        std::filesystem::create_directories(folder);
    }

    std::string imagePath =
        folder + imageName;

    cv::Mat newImage =
        cv::Mat::zeros(height, width, CV_8UC3);

    for (const auto& star : pixelCoordinates)
    {
        int radiusOfStar =
            2 + static_cast<int>(star.brightness / 80.0);

        radiusOfStar =
            std::clamp(radiusOfStar, 1, 5);

        cv::circle(
            newImage,
            cv::Point(
                static_cast<int>(std::round(star.u)),
                static_cast<int>(std::round(star.v))
            ),
            radiusOfStar,
            cv::Scalar(255, 255, 255),
            -1
        );
    }

    bool ok =
        cv::imwrite(imagePath, newImage);

    if (!ok)
    {
        std::cerr << "ERROR: Cannot write image: "
            << imagePath << std::endl;
    }
    else
    {
        std::cout << "Generated image: "
            << imagePath << std::endl;
    }
}

void ShowImage(const std::string& imageName, const std::string& windowName, const std::string& folder)
{
    std::string imagePath = folder + imageName;
    cv::Mat imageToShow = cv::imread(imagePath);
    if (imageToShow.empty())
    {
        std::cout << "Error loading Image" << std::endl;
        return;
    }
    cv::imshow(windowName, imageToShow);
}


