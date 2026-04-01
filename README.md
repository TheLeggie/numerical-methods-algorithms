# Numerical Methods Algorithms

This repository contains MATLAB scripts implementing numerical methods. The algorithms were developed as part of Numerical Methods course at the Warsaw University of Technology.

## Numerical Integration

* **`triangular_quadrature.m`**: Calculates the approximate value of a double integral for a function over a triangle defined by three vertices.
It implements a composite 3rd-order quadrature (midpoint rule) by dividing the main triangle into congruent sub-triangles.

## Eigenvalues and QR Decomposition

* **`householder.m`**: Performs the QR decomposition of a square matrix using Householder reflections, returning an orthogonal matrix Q and an upper triangular matrix R.
* **`inv_power_method_with_qr.m`**: Finds the eigenvalue with the smallest absolute value of a square matrix, along with its corresponding eigenvector and relative error estimate.
It utilizes the inverse power method combined with the QR decomposition obtained via Householder reflections.
