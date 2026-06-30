#include <algorithm>
#include <iostream>
#include <cmath>
//Attitude Estimation algorithm 
//this session should be used after the star identification is done 
//and the matched stars are found

#include "AttitudeEstimation.h"

constexpr double ATT_PI = 3.14159265358979323846;

static cv::Vec3d NormalizeVector(const cv::Vec3d& v)
{
    double n = cv::norm(v);

    if (n == 0.0)
        return cv::Vec3d(0.0, 0.0, 0.0);

    return v / n;
}

Quaternion NormalizeQuaternion(const Quaternion& q)
{
    double norm = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);

    return {
        q.w / norm,
        q.x / norm,
        q.y / norm,
        q.z / norm
    };
}

//converter:: rotation matrix to quaternion
static Quaternion RotationMatrixToQuaternion(const cv::Matx33d& R)
{
    double trace =
        R(0, 0) + R(1, 1) + R(2, 2);

    double qw, qx, qy, qz;

    if (trace > 0.0)
    {
        double s = std::sqrt(trace + 1.0) * 2.0;

        qw = 0.25 * s;
        qx = (R(2, 1) - R(1, 2)) / s;
        qy = (R(0, 2) - R(2, 0)) / s;
        qz = (R(1, 0) - R(0, 1)) / s;
    }
    else if ((R(0, 0) > R(1, 1)) &&
        (R(0, 0) > R(2, 2)))
    {
        double s =
            std::sqrt(
                1.0 + R(0, 0) - R(1, 1) - R(2, 2)
            ) * 2.0;

        qw = (R(2, 1) - R(1, 2)) / s;
        qx = 0.25 * s;
        qy = (R(0, 1) + R(1, 0)) / s;
        qz = (R(0, 2) + R(2, 0)) / s;
    }
    else if (R(1, 1) > R(2, 2))
    {
        double s =
            std::sqrt(
                1.0 + R(1, 1) - R(0, 0) - R(2, 2)
            ) * 2.0;

        qw = (R(0, 2) - R(2, 0)) / s;
        qx = (R(0, 1) + R(1, 0)) / s;
        qy = 0.25 * s;
        qz = (R(1, 2) + R(2, 1)) / s;
    }
    else
    {
        double s =
            std::sqrt(
                1.0 + R(2, 2) - R(0, 0) - R(1, 1)
            ) * 2.0;

        qw = (R(1, 0) - R(0, 1)) / s;
        qx = (R(0, 2) + R(2, 0)) / s;
        qy = (R(1, 2) + R(2, 1)) / s;
        qz = 0.25 * s;
    }

    Quaternion q = Quaternion{ qw, qx, qy, qz };
    q = NormalizeQuaternion(q);

    // keep scalar part positive 
    if (q.w < 0.0)
    {
        q.w = -q.w;
        q.x = -q.x;
        q.y = -q.y;
        q.z = -q.z;
    }

    return q;
}

cv::Matx33d QuaternionToRotationMatrix(const Quaternion& q_in)
{
    Quaternion q = NormalizeQuaternion(q_in);

    double w = q.w;
    double x = q.x;
    double y = q.y;
    double z = q.z;

    cv::Matx33d R;

    R(0, 0) = 1.0 - 2.0 * (y * y + z * z);
    R(0, 1) = 2.0 * (x * y - z * w);
    R(0, 2) = 2.0 * (x * z + y * w);

    R(1, 0) = 2.0 * (x * y + z * w);
    R(1, 1) = 1.0 - 2.0 * (x * x + z * z);
    R(1, 2) = 2.0 * (y * z - x * w);

    R(2, 0) = 2.0 * (x * z - y * w);
    R(2, 1) = 2.0 * (y * z + x * w);
    R(2, 2) = 1.0 - 2.0 * (x * x + y * y);

    return R;
}

static double AngularErrorDeg(
    const cv::Vec3d& measured,
    const cv::Vec3d& predicted)
{
    double dot =
        measured[0] * predicted[0] +
        measured[1] * predicted[1] +
        measured[2] * predicted[2];

    dot = std::clamp(dot, -1.0, 1.0);

    double angleRad = std::acos(dot);

    return angleRad * 180.0 / ATT_PI;
}

cv::Matx33d CalculateAttitudeProfileMatrix(const std::vector<MatchedStar>& matchedStars)
{
    cv::Matx33d B = cv::Matx33d::zeros();

    for (const auto& star : matchedStars)
    {
        cv::Vec3d b =
            NormalizeVector(
                cv::Vec3d(star.bx, star.by, star.bz)
            );

        cv::Vec3d r =
            NormalizeVector(
                cv::Vec3d(star.rx, star.ry, star.rz)
            );

        // B Attitude profile matrix = sum( b_i * r_i^T ) 
        B += cv::Matx33d(
            b[0] * r[0], b[0] * r[1], b[0] * r[2],
            b[1] * r[0], b[1] * r[1], b[1] * r[2],
            b[2] * r[0], b[2] * r[1], b[2] * r[2]
        );
    }

    return B;
}

//svd based attitude solver
AttitudeResult EstimateAttitudeSVD(
    const std::vector<MatchedStar>& matchedStars)
{
    AttitudeResult result{};

    if (matchedStars.size() < 3)
    {
        std::cerr << "ERROR: At least 3 matched stars are required for attitude estimation.\n";
        return result;
    }

    cv::Matx33d B = CalculateAttitudeProfileMatrix(matchedStars);

    cv::Mat Bmat(3, 3, CV_64F);

    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            Bmat.at<double>(row, col) = B(row, col);
        }
    }

    cv::Mat W;
    cv::Mat U;
    cv::Mat Vt;

    cv::SVD::compute(Bmat, W, U, Vt);

    cv::Mat Rmat = U * Vt;

    // Determinant correction
    if (cv::determinant(Rmat) < 0.0)
    {
        U.at<double>(0, 2) *= -1.0;
        U.at<double>(1, 2) *= -1.0;
        U.at<double>(2, 2) *= -1.0;

        Rmat = U * Vt;
    }

    cv::Matx33d R;

    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            R(row, col) = Rmat.at<double>(row, col);
        }
    }

    result.R = R;
    result.q = RotationMatrixToQuaternion(R);
    result.usedStars = static_cast<int>(matchedStars.size());
    result.success = true;

    double sumSquaredError = 0.0;

    for (const auto& star : matchedStars)
    {
        cv::Vec3d bMeasured =
            NormalizeVector(
                cv::Vec3d(star.bx, star.by, star.bz)
            );

        cv::Vec3d rCatalog =
            NormalizeVector(
                cv::Vec3d(star.rx, star.ry, star.rz)
            );

        // predicted camera vector = R * catalog vector
        cv::Vec3d bPredicted = R * rCatalog;
        bPredicted = NormalizeVector(bPredicted);

        double errorDeg =
            AngularErrorDeg(bMeasured, bPredicted);

        sumSquaredError += errorDeg * errorDeg;
    }

    result.rmsResidualDeg =
        std::sqrt(sumSquaredError / matchedStars.size());

    return result;
}

//QUEST algorithm attitude solver
AttitudeResult EstimateAttitudeQUEST(
    const std::vector<MatchedStar>& matchedStars)
{
    AttitudeResult result{};

    if (matchedStars.size() < 3)
    {
        std::cerr << "ERROR: At least 3 matched stars are required.\n";
        return result;
    }

    cv::Matx33d B = CalculateAttitudeProfileMatrix(matchedStars);

    // Davenport K matrix
    cv::Matx33d S = B + B.t();

    double sigma =
        B(0, 0) + B(1, 1) + B(2, 2);

    cv::Vec3d Z(
        B(1, 2) - B(2, 1),
        B(2, 0) - B(0, 2),
        B(0, 1) - B(1, 0));

    cv::Mat K = cv::Mat::zeros(4, 4, CV_64F);

    // upper-left block
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            K.at<double>(i, j) =
                S(i, j) - (i == j ? sigma : 0.0);
        }
    }

    // upper-right
    K.at<double>(0, 3) = Z[0];
    K.at<double>(1, 3) = Z[1];
    K.at<double>(2, 3) = Z[2];

    // lower-left
    K.at<double>(3, 0) = Z[0];
    K.at<double>(3, 1) = Z[1];
    K.at<double>(3, 2) = Z[2];

    // lower-right
    K.at<double>(3, 3) = sigma;

    // biggest eigenvector = optimal quaternion
    cv::Mat eigenValues;
    cv::Mat eigenVectors;

    cv::eigen(K, eigenValues, eigenVectors);
    //std::cout << "Eigenvalues:\n" << eigenValues << std::endl;
    //std::cout << "Eigenvector row 0:\n";
    //for (int i = 0; i < 4; i++)
    //    std::cout << eigenVectors.at<double>(0, i) << " ";

    // eigenvectors as scalar last
    Quaternion q;

    q.x = eigenVectors.at<double>(0, 0);
    q.y = eigenVectors.at<double>(0, 1);
    q.z = eigenVectors.at<double>(0, 2);
    q.w = eigenVectors.at<double>(0, 3);

    q = NormalizeQuaternion(q);

    if (q.w < 0.0)
    {
        q.w = -q.w;
        q.x = -q.x;
        q.y = -q.y;
        q.z = -q.z;
    }

    //Quaternion to rotation matrix
    result.R = QuaternionToRotationMatrix(q).t();
    result.q = RotationMatrixToQuaternion(result.R);

    result.usedStars =
        static_cast<int>(matchedStars.size());

    result.success = true;

    // RMS residual
    double sumSquaredError = 0.0;

    for (const auto& star : matchedStars)
    {
        cv::Vec3d bMeasured =
            NormalizeVector(
                cv::Vec3d(star.bx, star.by, star.bz));

        cv::Vec3d rCatalog =
            NormalizeVector(
                cv::Vec3d(star.rx, star.ry, star.rz));

        //body to inertial frame rotation
        cv::Vec3d bPredicted =
            result.R * rCatalog;

        bPredicted =
            NormalizeVector(bPredicted);

        double errorDeg =
            AngularErrorDeg(
                bMeasured,
                bPredicted);

        sumSquaredError +=
            errorDeg * errorDeg;
    }

    result.rmsResidualDeg =
        std::sqrt(
            sumSquaredError /
            matchedStars.size());

    return result;
}

//visualization of the attitude estimation result
void PrintAttitudeResult(const AttitudeResult& result)
{
    if (!result.success)
    {
        std::cout << "\nAttitude estimation failed.\n";
        return;
    }

    std::cout << "Used matched stars = "
        << result.usedStars << "\n";

    std::cout << "\nRotation Matrix R "
        << "(catalog/inertial -> camera):\n";

    std::cout << result.R(0, 0) << "  "
        << result.R(0, 1) << "  "
        << result.R(0, 2) << "\n";

    std::cout << result.R(1, 0) << "  "
        << result.R(1, 1) << "  "
        << result.R(1, 2) << "\n";

    std::cout << result.R(2, 0) << "  "
        << result.R(2, 1) << "  "
        << result.R(2, 2) << "\n";

    std::cout << "\nQuaternion q_est [qw, qx, qy, qz]:\n";
    Quaternion q = result.q;
    std::cout << "["
        << q.w << ", "
        << q.x << ", "
        << q.y << ", "
        << q.z << "]\n";

    std::cout << "\nRMS star residual error = "
        << result.rmsResidualDeg
        << " deg\n";
}

double QuaternionAngleErrorDeg(
    const Quaternion& qEstimated,
    const Quaternion& qTruth)
{
    Quaternion qe = NormalizeQuaternion(qEstimated);
    Quaternion qt = NormalizeQuaternion(qTruth);

    double dot =
        qe.w * qt.w +
        qe.x * qt.x +
        qe.y * qt.y +
        qe.z * qt.z;

    //same attitude
    dot = std::abs(dot);

    dot = std::clamp(dot, -1.0, 1.0);

    double angleRad = 2.0 * std::acos(dot);

    return angleRad * 180.0 / ATT_PI;
}

std::tuple<double, double, double> QuaternionToEular(Quaternion& q_est)
{
    double roll, pitch, yaw;

    //roll (x-axis)
    double sinr_cosp = 2.0 * (q_est.w * q_est.x + q_est.y * q_est.z);
    double cosr_cosp = 1.0 - 2.0 * (q_est.x * q_est.x + q_est.y * q_est.y);
    roll = std::atan2(sinr_cosp, cosr_cosp);

    //pitch (y-axis)
    double sinp = 2.0 * (q_est.w * q_est.y - q_est.z * q_est.x);
    if (std::abs(sinp) >= 1)
        pitch = std::copysign(ATT_PI / 2, sinp); // 90 degrees clamp if out of range
    else
        pitch = std::asin(sinp);

    //yaw (z-axis)
    double siny_cosp = 2.0 * (q_est.w * q_est.z + q_est.x * q_est.y);
    double cosy_cosp = 1.0 - 2.0 * (q_est.y * q_est.y + q_est.z * q_est.z);
    yaw = std::atan2(siny_cosp, cosy_cosp);

    return std::make_tuple(roll, pitch, yaw);
}