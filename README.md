# Catalog-Based Star Tracker

## Overview
Short explanation of what this project does.

## Main Features
- Star image preprocessing
- Star centroid detection
- Camera unit vector conversion
- Hipparcos catalog matching
- Pyramid-based star identification
- Reprojection verification
- Attitude estimation using SVD / QUEST

## Project Pipeline
Input Image → Preprocessing → Star Detection → Pyramid Matching → Reprojection Verification → Attitude Estimation

## Requirements
- C++
- OpenCV
- Visual Studio / CMake
- Hipparcos binary catalog files

## How to Run
Step-by-step command or Visual Studio steps.

## Example Result
A C++/OpenCV catalog-based star tracker pipeline for synthetic star field generation, star detection, candidate labeling, catalog matching, and attitude estimation.

<p align="center">
  <img src="https://raw.githubusercontent.com/ShineEAunG/Catalog-Based-Star-Tracker/main/doc/image.png" alt="Terminal output" width="70">
</p>

<p align="center">
  <em>Example output in Terminal for catalog-based star tracker pipeline project.</em>
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/ShineEAunG/Catalog-Based-Star-Tracker/main/doc/Synthetic_Star_Field.png" alt="Generated synthetic star field" width="45%">
  <img src="https://raw.githubusercontent.com/ShineEAunG/Catalog-Based-Star-Tracker/main/doc/labeled_candidate_stars.png" alt="Labeled star identification result" width="45%">
</p>

## Future Work
- Raspberry Pi 5 implementation
- Real sky image testing
- Lens distortion correction
- Faster pair database search
