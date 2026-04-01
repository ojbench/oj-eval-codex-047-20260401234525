// Implementation for Controller::get_v_next()
// Notes:
// - We assume the judge provides declarations for Controller, Vec, and Monitor,
//   as described in the problem statement. Here we only define the member
//   function body, relying on those declarations.
// - Strategy: Move toward target within speed limit, with simple collision
//   avoidance by damping approach along line-of-centers if a predicted close
//   pass within the next time interval could occur (using others' current
//   velocities as estimates). This keeps the solution safe and robust.

#ifndef TIME_INTERVAL
#define TIME_INTERVAL 0.1
#endif

#include <cmath>
#include <algorithm>

// Forward declarations to satisfy compilers that require awareness of types.
// The real definitions are expected to come from the judge's framework.
struct Vec;
class Monitor;
class Controller;

// Helper inline utilities operating on Vec assuming x/y fields and basic
// arithmetic operators (+, -, *, /) are available from the framework.
inline double vec_dot(const Vec &a, const Vec &b) { return a.x * b.x + a.y * b.y; }
inline double vec_len2(const Vec &a) { return a.x * a.x + a.y * a.y; }
inline double vec_len(const Vec &a) { return std::sqrt(vec_len2(a)); }
inline Vec vec_clamp_len(const Vec &v, double maxlen) {
    double l2 = vec_len2(v);
    if (l2 <= 0.0) return v;
    double l = std::sqrt(l2);
    if (l <= maxlen) return v;
    return v * (maxlen / l);
}

// Define the method outside of class; the class must be declared by the framework.
// Signature assumed as: Vec Controller::get_v_next() const or non-const; we avoid const.
// We match common pattern: Vec Controller::get_v_next() {

Vec Controller::get_v_next() {
    // Read own state
    const Vec p = this->pos_cur;
    const Vec vcur = this->v_cur;
    const Vec ptar = this->pos_tar;
    const double r_self = this->r;
    const double vmax = this->v_max;

    const double eps = 1e-9;
    const double tau = TIME_INTERVAL;
    const double margin = 1e-3; // small safety padding

    // Desired velocity towards target
    Vec to_tar = ptar - p;
    double dist = vec_len(to_tar);
    Vec v_des;
    if (dist <= eps) {
        v_des = Vec{0.0, 0.0};
    } else {
        double desired_speed = std::min(vmax, dist / std::max(tau, eps));
        v_des = to_tar / dist * desired_speed;
    }

    // If recently speeding, be conservative (optional heuristic)
    if (monitor && monitor->get_speeding(this->id)) {
        v_des = v_des * 0.8;
    }

    // Collision avoidance: damp approach along line-of-centers for close neighbors
    Vec v = v_des;
    if (monitor) {
        int N = monitor->get_robot_number();
        // Iterate a couple passes to reduce predicted conflicts
        for (int pass = 0; pass < 2; ++pass) {
            for (int j = 0; j < N; ++j) {
                if (j == this->id) continue;
                Vec pj = monitor->get_pos_cur(j);
                Vec vj = monitor->get_v_cur(j);
                double rj = monitor->get_r(j);

                Vec rel_p = p - pj;
                double d0 = vec_len(rel_p);
                if (d0 <= eps) continue;

                // Skip far neighbors for efficiency
                double interact_range = r_self + rj + std::max(vmax, vec_len(vj)) * tau * 2.0 + 1e-6;
                if (d0 > interact_range) continue;

                Vec rel_v = v - vj;
                double A = vec_len2(rel_v);
                double tstar = 0.0;
                if (A > eps) {
                    double B = vec_dot(rel_p, rel_v);
                    tstar = -B / A;
                    if (tstar < 0.0) tstar = 0.0;
                    if (tstar > tau) tstar = tau;
                } else {
                    tstar = 0.0;
                }

                Vec at_min = rel_p + rel_v * tstar;
                double dmin = vec_len(at_min);
                double Rsum = r_self + rj + margin;
                if (dmin <= Rsum) {
                    // Compute unit along current separation
                    Vec u = rel_p / d0;
                    // Approach rate along u (negative means approaching)
                    double appr = vec_dot(rel_v, u);

                    // Target non-approach small positive to open distance
                    double target = std::max(0.0, 0.1 * vmax);
                    double delta = target - appr; // we add delta * u to our v
                    v = v + u * delta;

                    // Re-clamp to speed limit
                    v = vec_clamp_len(v, vmax);
                }

                // If starting too close, gently separate
                if (d0 < Rsum) {
                    Vec u = rel_p / d0;
                    v = v + u * std::min(0.2 * vmax, (Rsum - d0) / std::max(tau, eps));
                    v = vec_clamp_len(v, vmax);
                }
            }
        }
    }

    // Final clamp and return
    v = vec_clamp_len(v, vmax);
    return v;
}

