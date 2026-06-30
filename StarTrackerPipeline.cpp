#include "StarTrackerPipeline.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>

#include "CatalogGenerator.h"
#include "PreprocessingImage.h"
#include "StarIdentification.h"
#include "AttitudeEstimation.h"
#include "ReprojectionVerification.h"
#include "Validation.h"
#include "Pattern_Generator.h"
#include "TestImageGeneration.h"
#include "LabelizationImage.h"

struct CameraSettings
{
    int width = 4056;
    int height = 3040;

    double focalLengthMm = 4.0;
    double pixelSizeMm = 0.0015;

    double FocalLengthPixels() const
    {
        return focalLengthMm / pixelSizeMm;
    }
};

struct AppPaths
{
    std::string catalogBinaryPath =
        "../Hipparcos_Catalog_Binaries/hip_catalog_mag6.bin";

    std::string pairDatabasePath =
        "../Hipparcos_Catalog_Binaries/hip_pair_mag6.bin";

    std::string imageFolder =
        "../ImgFolder/";

    std::string subFolder =
        "../ImgFolder/subFolder/";

    std::string outputFolder =
        "../ImgFolder/Output/";

    std::string inputImageName =
        "Synthetic_Star_Field.png";

    std::string truthMetadataName =
        "Synthetic_Star_Field_Metadata.txt";

    std::string InputImagePath() const
    {
        return imageFolder + inputImageName;
    }

    std::string TruthMetadataPath() const
    {
        return imageFolder + truthMetadataName;
    }

    std::string LabeledDetectedImagePath() const
    {
        return outputFolder + "labeled_detected_stars.png";
    }

    std::string LabeledCandidateImagePath() const
    {
        return outputFolder + "labeled_candidate_stars.png";
    }
};

struct IdentificationSettings
{
    int thresholdLevel = 39;
    double toleranceDeg = 0.05;
    int maxDetectedStars = 15;

    int topK = 30;
    double pixelThreshold = 5.0;
    double reprojectionMagnitudeLimit = 6.0;
};

struct SyntheticImageSettings
{
    double roll = -0.0689274;
    double pitch = 0.5808292;
    double yaw = 1.4456391;

    double maxBrightness = 255.0;
    double generationMagnitudeLimit = 2.5;

    std::string outputImageName =
        "Synthetic_Star_Field.png";

    std::string outputMetadataName =
        "Synthetic_Star_Field_Metadata.txt";
};

// ============================================================
// Quaternion / attitude helper
// ============================================================

static Quaternion EulerToQuaternion(
    double roll,
    double pitch,
    double yaw)
{
    Quaternion q;

    double cy = std::cos(yaw * 0.5);
    double sy = std::sin(yaw * 0.5);

    double cr = std::cos(roll * 0.5);
    double sr = std::sin(roll * 0.5);

    double cp = std::cos(pitch * 0.5);
    double sp = std::sin(pitch * 0.5);

    q.w = cy * cr * cp + sy * sr * sp;
    q.x = cy * sr * cp - sy * cr * sp;
    q.y = cy * cr * sp + sy * sr * cp;
    q.z = sy * cr * cp - cy * sr * sp;

    return NormalizeQuaternion(q);
}

static Quaternion BuildTruthQuaternion(
    const SyntheticImageSettings& settings)
{
    return EulerToQuaternion(
        settings.roll,
        settings.pitch,
        settings.yaw);
}

static void PrintTruthAttitude(
    const SyntheticImageSettings& settings,
    const Quaternion& qTruth)
{
    std::cout << "Euler Angle"
        << "\nRoll  = " << settings.roll
        << "\nPitch = " << settings.pitch
        << "\nYaw   = " << settings.yaw
        << "\n";

    std::cout << "Quaternion"
        << "\nw = " << qTruth.w
        << "\nx = " << qTruth.x
        << "\ny = " << qTruth.y
        << "\nz = " << qTruth.z
        << "\n";
}

// ============================================================
// General helpers
// ============================================================

static bool EnsureFolderExists(
    const std::string& folderPath)
{
    if (CheckFolderExists(folderPath))
        return true;

    try
    {
        return std::filesystem::create_directory(folderPath);
    }
    catch (const std::exception& e)
    {
        std::cerr << "ERROR: Cannot create folder: "
            << folderPath << "\n";

        std::cerr << e.what() << "\n";

        return false;
    }
}

static char AskUserOption()
{
    char userInput;

    std::cout << "Select an option:\n";
    std::cout << "1. Run Star Identification and Attitude Estimation\n";
    std::cout << "2. Test Image Generation\n";
    std::cin >> userInput;

    return userInput;
}

static CameraModel BuildCameraModel(
    const CameraSettings& cameraSettings,
    const StarIdentificationDebugData& identification)
{
    CameraModel camera;

    camera.fx = cameraSettings.FocalLengthPixels();
    camera.fy = cameraSettings.FocalLengthPixels();

    camera.width = identification.imageWidth;
    camera.height = identification.imageHeight;

    camera.cx = identification.imageWidth / 2.0;
    camera.cy = identification.imageHeight / 2.0;

    return camera;
}

// ============================================================
// Pipeline stage: preprocessing
// ============================================================

static bool RunPreprocessing(
    const AppPaths& paths,
    const IdentificationSettings& settings,
    std::string& binaryImagePath,
    std::string& binaryImageName)
{
    if (!CheckFolderExists(paths.imageFolder))
    {
        std::cerr << "ERROR: Image folder is not found.\n";
        return false;
    }

    auto result =
        PreprocessInputImage(
            paths.InputImagePath(),
            paths.imageFolder,
            paths.subFolder,
            settings.thresholdLevel);

    std::tie(binaryImagePath, binaryImageName) = result;

    if (binaryImagePath.empty())
    {
        std::cerr << "ERROR: Image preprocessing failed.\n";
        return false;
    }

    return true;
}

// ============================================================
// Pipeline stage: star identification
// ============================================================

static bool RunStarIdentification(
    const std::string& binaryImagePath,
    const std::vector<StarPairAngle>& pairDatabase,
    const CameraSettings& cameraSettings,
    const IdentificationSettings& settings,
    StarIdentificationDebugData& identification)
{
    identification =
        IdentifyAllPyramidCandidatesForDebugging(
            binaryImagePath,
            pairDatabase,
            cameraSettings.FocalLengthPixels(),
            settings.toleranceDeg,
            settings.maxDetectedStars);

    if (!identification.success)
    {
        std::cerr << "Star identification failed.\n";
        return false;
    }

    std::cout << "\nTotal pyramid candidates = "
        << identification.candidates.size()
        << "\n";

    return true;
}

// ============================================================
// Pipeline stage: reprojection verification
// ============================================================

static bool RunReprojectionVerification(
    const std::vector<CatalogStar>& catalogStars,
    const StarIdentificationDebugData& identification,
    const CameraSettings& cameraSettings,
    const IdentificationSettings& settings,
    ReprojectionResult& best)
{
    if (catalogStars.empty())
    {
        std::cerr << "ERROR: Catalog is empty.\n";
        return false;
    }

    CameraModel camera =
        BuildCameraModel(
            cameraSettings,
            identification);

    best =
        SelectBestCandidateByReprojection(
            identification.candidates,
            identification.cameraStars,
            identification.centroids,
            catalogStars,
            camera,
            settings.topK,
            settings.pixelThreshold,
            settings.reprojectionMagnitudeLimit);

    PrintReprojectionResult(best);

    if (!best.success)
    {
        std::cerr << "Reprojection verification failed.\n";
        return false;
    }

    return true;
}

// ============================================================
// Pipeline stage: attitude estimation output
// ============================================================

static void PrintFinalAttitudeResults(
    const ReprojectionResult& best,
    AttitudeResult& attitudeSVD,
    AttitudeResult& attitudeQUEST)
{
    attitudeSVD = best.attitude;

    attitudeQUEST =
        EstimateAttitudeQUEST(best.matchedStars);

    std::cout << "\n***********SVD Attitude Solver***********\n";
    PrintAttitudeResult(attitudeSVD);

    auto [rollSVD, pitchSVD, yawSVD] =
        QuaternionToEular(attitudeSVD.q);

    std::cout << "roll: " << rollSVD
        << " rad, pitch: " << pitchSVD
        << " rad, yaw: " << yawSVD
        << " rad\n";

    std::cout << "\n***********QUEST Attitude Solver***********\n";
    PrintAttitudeResult(attitudeQUEST);

    auto [rollQUEST, pitchQUEST, yawQUEST] =
        QuaternionToEular(attitudeQUEST.q);

    std::cout << "roll: " << rollQUEST
        << " rad, pitch: " << pitchQUEST
        << " rad, yaw: " << yawQUEST
        << " rad\n";
}

// ============================================================
// Pipeline stage: validation
// ============================================================

static void RunSyntheticValidation(
    const AppPaths& paths,
    const ReprojectionResult& best,
    const AttitudeResult& attitudeSVD,
    const AttitudeResult& attitudeQUEST,
    const Quaternion& qTruth)
{
    std::vector<TruthStar> truthStars =
        LoadTruthMetadataCSV(
            paths.TruthMetadataPath());

    bool hasTruthQuaternion = true;

    ValidationSummary validation =
        ValidateSyntheticResult(
            best,
            truthStars,
            attitudeSVD,
            attitudeQUEST,
            hasTruthQuaternion,
            qTruth);

    PrintValidationSummary(validation);
}

// ============================================================
// Pipeline stage: visualization
// ============================================================

static void SaveVisualizationImages(
    const AppPaths& paths,
    const StarIdentificationDebugData& identification,
    const ReprojectionResult* best)
{
    SaveDetectedStarLabelImage(
        paths.InputImagePath(),
        identification.centroids,
        paths.LabeledDetectedImagePath());

    if (best != nullptr && best->success)
    {
        SaveCandidateStarLabelImage(
            paths.InputImagePath(),
            identification.centroids,
            best->candidate,
            paths.LabeledCandidateImagePath());
    }
}

// ============================================================
// Full pipeline: option 1
// ============================================================

static int RunIdentificationAndAttitudePipeline(
    const AppPaths& paths,
    const CameraSettings& cameraSettings,
    const IdentificationSettings& identificationSettings,
    const Quaternion& qTruth,
    const std::vector<CatalogStar>& catalogStars)
{
    std::vector<StarPairAngle> pairDatabase =
        LoadPairDatabaseBinary(
            paths.pairDatabasePath);

    if (pairDatabase.empty())
    {
        std::cerr << "ERROR: Pair database loading failed.\n";
        return -1;
    }

    std::string binaryImagePath;
    std::string binaryImageName;

    if (!RunPreprocessing(
        paths,
        identificationSettings,
        binaryImagePath,
        binaryImageName))
    {
        return -1;
    }

    StarIdentificationDebugData identification;

    if (!RunStarIdentification(
        binaryImagePath,
        pairDatabase,
        cameraSettings,
        identificationSettings,
        identification))
    {
        return -1;
    }

    SaveVisualizationImages(
        paths,
        identification,
        nullptr);

    ReprojectionResult best;

    if (!RunReprojectionVerification(
        catalogStars,
        identification,
        cameraSettings,
        identificationSettings,
        best))
    {
        return -1;
    }

    SaveVisualizationImages(
        paths,
        identification,
        &best);

    AttitudeResult attitudeSVD;
    AttitudeResult attitudeQUEST;

    PrintFinalAttitudeResults(
        best,
        attitudeSVD,
        attitudeQUEST);

    RunSyntheticValidation(
        paths,
        best,
        attitudeSVD,
        attitudeQUEST,
        qTruth);

    return 0;
}

// ============================================================
// Full pipeline: option 2
// ============================================================

static int RunTestImageGenerationPipeline(
    const AppPaths& paths,
    const CameraSettings& cameraSettings,
    const SyntheticImageSettings& syntheticSettings,
    const Quaternion& qTruth,
    const std::vector<CatalogStar>& catalogStars)
{
    if (catalogStars.empty())
    {
        std::cerr << "ERROR: Catalog is empty.\n";
        return -1;
    }

    PrintTruthAttitude(
        syntheticSettings,
        qTruth);

    std::vector<CameraTruthStar> cameraFrameStars =
        ConvertCatalogIntoCameraFrame(
            catalogStars,
            qTruth,
            syntheticSettings.maxBrightness,
            cameraSettings.width,
            cameraSettings.height,
            cameraSettings.FocalLengthPixels(),
            syntheticSettings.generationMagnitudeLimit);

    std::vector<PixelCoordinate> pixelCoordinates =
        ConvertIntoPixelCoordinates(
            cameraFrameStars,
            cameraSettings.width,
            cameraSettings.height,
            cameraSettings.FocalLengthPixels());

    Generate3ColorsConstellation(
        pixelCoordinates,
        cameraSettings.width,
        cameraSettings.height,
        syntheticSettings.outputImageName,
        paths.imageFolder);

    WriteCoordinatesToFile(
        pixelCoordinates,
        syntheticSettings.outputMetadataName,
        paths.imageFolder);

    return 0;
}

int RunStarTrackerApplication()
{
    std::cout << "Catalog-Based Star Identification\n";

    AppPaths paths;
    CameraSettings cameraSettings;
    IdentificationSettings identificationSettings;
    SyntheticImageSettings syntheticSettings;

    Quaternion qTruth =
        BuildTruthQuaternion(
            syntheticSettings);

    if (!EnsureFolderExists(paths.outputFolder))
    {
        return -1;
    }

    std::vector<CatalogStar> catalogStars =
        LoadCatalogBinary(
            paths.catalogBinaryPath);

    if (catalogStars.empty())
    {
        std::cerr << "ERROR: Catalog loading failed.\n";
        return -1;
    }

    char userInput =
        AskUserOption();

    if (userInput == '1')
    {
        return RunIdentificationAndAttitudePipeline(
            paths,
            cameraSettings,
            identificationSettings,
            qTruth,
            catalogStars);
    }

    if (userInput == '2')
    {
        return RunTestImageGenerationPipeline(
            paths,
            cameraSettings,
            syntheticSettings,
            qTruth,
            catalogStars);
    }

    std::cerr << "Invalid option selected.\n";

    return -1;
}