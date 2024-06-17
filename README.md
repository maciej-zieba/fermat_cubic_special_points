# Sextactic and Type-9 Points on the Fermat Cubic

**Authors:** Łukasz Merta and Maciej Ziȩba

**Date:** 14.06.2023

## Introduction

This repository contains code and results associated with the article titled "Sextactic and Type-9 Points on the Fermat Cubic and Associated Objects."

## Prerequisites

The code is written using Singular 4.1.1. To run the code, you need to install Singular and use code from current repository.

## Structure of the repository

- README.md - this file
- data/

  - 6_type_points.txt - coordinates of type-6 points
  - 9_type_points.txt - coordinates of type-9 points
  - 9_type_cycles.txt - coordinates of triplets of type-9 points created from tangents
  - 24_deg_curve.txt - equation of the 24-degree curve passing through type-9 points
  - 72_points_where_nine_conics_intersect.txt - coordinates of points where nine conics intersect

- src/

  - main.cpp - main file containing calls for procedures
  - procedures.cpp - procedures used in the main file
