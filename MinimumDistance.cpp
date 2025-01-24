/**
 * @file MinimumDistance.cpp
 * @author Euncheol Im
 * @contact euncheol.im@kist.re.kr
 * @brief Implementation of an algorithm to compute the minimum distance between two 3D line segments.
 * @version 0.1
 * @date 2025-01-24
 * 
 * @copyright Copyright (c) 2025
 * 
 * @details
 * This implementation is based on the algorithm described in the following paper:
 * 
 * [Paper Title] : ON FAST COMPUTATION OF DISTANCE BETWEEN LINE SEGMENTS
 * [Author(s)] : Vladimir J. LUMELSKY
 * [Publication Year] August 16th, 1985
 * [Journal or Conference Name] Information Processing Letters 21 (1985) 55-61, North-Holland
 */


#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <algorithm>

using Eigen::Vector3d;

// Custom clamp function
template <typename T>
T clamp(const T& value, const T& min, const T& max) {
    return std::max(min, std::min(max, value));
}
double computeMinDistance3D(
    const Vector3d& start1, const Vector3d& end1,
    const Vector3d& start2, const Vector3d& end2,
    double& t, double& u) {
    // Compute direction vectors and differences
    Vector3d d1 = end1 - start1; // Direction vector of segment 1 (AB)
    Vector3d d2 = end2 - start2; // Direction vector of segment 2 (CD)
    Vector3d r = start2 - start1;

    // Compute  D1, D2, R, S1, S2
    double D1 = d1.dot(d1);
    double D2 = d2.dot(d2);
    double R = d1.dot(d2);
    double S1 = d1.dot(r);
    double S2 = d2.dot(r);

    // denominator
    double denominator = D1 * D2 - R * R;

    // Step 1 IF parallel ?
    if (D1 == 0 && D2 > 0) { // AB degenerates into a point
        t = 0;
        u = clamp(-S2 / D2, 0.0, 1.0);
        return (start1 - (start2 + u * d2)).squaredNorm();
    }
    if (D2 == 0 && D1 > 0) { // CD degenerates into a point
        u = 0;
        t = clamp(S1 / D1, 0.0, 1.0);
        return (start2 - (start1 + t * d1)).squaredNorm();
    }
    if (D1 == 0 && D2 == 0) { // Both segments degenerate into points
        t = 0;
        u = 0;
        return r.squaredNorm();
    }
    if (denominator == 0) { // Segments are parallel
        t = 0;
        u = -S2 / D2;
        if (u < 0 || u > 1) {
            u = clamp(u, 0.0, 1.0); // Clamp u
            t = clamp((u * R + S1) / D1, 0.0, 1.0); // Recompute t (Step 4)
        }
        return (start1 + t * d1 - (start2 + u * d2)).squaredNorm();
    }

    // Step 2 eq. (11) --> t
    t = (S1 * D2 - S2 * R) / denominator;
    std::cout << "Step2 t: " << t<< ", u: "<<u << std::endl;
    t = clamp(t, 0.0, 1.0);
    std::cout << "Step2 t: " << t<< ", u: "<<u << std::endl;
    // Step 3 eq. (10) --> u
    // u = (-S2 * D1 + S1 * R) / denominator;
    u = (t*R - S2) / D2;
    std::cout << "Step3 t: " << t<< ", u: "<<u << std::endl;
    u = clamp(u, 0.0, 1.0);
    std::cout << "Step3 t: " << t<< ", u: "<<u << std::endl;
    // Step 4 eq. (10) --> t
    t = (u*R + S1)/D1;
    std::cout << "Step4 t: " << t<< ", u: "<<u << std::endl;
    t = clamp(t, 0.0, 1.0);
    std::cout << "Step4 t: " << t<< ", u: "<<u << std::endl;


    // Step 5: Compute the actual minimum distance
    Vector3d closest1 = start1 + t * d1;
    Vector3d closest2 = start2 + u * d2;
    return (closest1 - closest2).squaredNorm();
}

int main() {
    // Paper Examples
    Vector3d start1(0, 0, 0), end1(1, 2, 1); // Segment 1: A(0,0,0) to B(1,2,1)
    Vector3d start2(1, 0, 0), end2(2, 1, 0); // Segment 2: C(1,0,0) to D(2,1,0)

    double t, u;
    double minDistSquared = computeMinDistance3D(start1, end1, start2, end2, t, u);
    double minDist = std::sqrt(minDistSquared);

    std::cout << "Minimum distance squared: " << minDistSquared << std::endl;
    std::cout << "Minimum distance: " << minDist << std::endl;

    return 0;
}
