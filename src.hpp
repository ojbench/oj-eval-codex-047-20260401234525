// Implementation for Controller::get_v_next()
// Strategy: Move toward target within speed limit, with simple collision
// avoidance by damping approach along line-of-centers if a predicted close pass
// within the next time interval could occur. Uses only local helpers inside
// the function to avoid incomplete-type issues.

#ifndef TIME_INTERVAL
#define TIME_INTERVAL 0.1
#endif

#include <cmath>
#include <algorithm>

Vec Controller::get_v_next() {
    // Local helpers that assume Vec has public x, y and scalar ops.
    auto dot = [](const Vec &a, const Vec &b) { return a.x * b.x + a.y * b.y; };
    auto len2 = [&](const Vec &a) { return a.x * a.x + a.y * a.y; };
    auto len = [&](const Vec &a) { return std::sqrt(len2(a)); };
    auto clamp_len = [&](const Vec &v, double maxlen) {
        double l2 = len2(v);
        if (l2 <= 0.0) return v;
        double l = std::sqrt(l2);
        if (l <= maxlen) return v;
        return v * (maxlen / l);
    };

    const Vec p = this->pos_cur;
    const Vec vcur = this->v_cur;
    const Vec ptar = this->pos_tar;
    const double r_self = this->r;
    const double vmax = this->v_max;

    const double eps = 1e-9;
    const double tau = TIME_INTERVAL;
    const double margin = 1e-3;

    // Desired velocity towards target
    Vec to_tar = ptar - p;
    double dist = len(to_tar);
    Vec v_des = vcur * 0.0;
    if (dist > eps) {
        double desired_speed = std::min(vmax, dist / std::max(tau, eps));
        v_des = (to_tar / dist) * desired_speed;
    }

    // If recently speeding, be conservative
    if (monitor && monitor->get_speeding(this->id)) {
        v_des = v_des * 0.8;
    }

    // Collision avoidance: damp approach along line-of-centers for close neighbors
    Vec v = v_des;
    if (monitor) {
        int N = monitor->get_robot_number();
        for (int pass = 0; pass < 2; ++pass) {
            for (int j = 0; j < N; ++j) {
                if (j == this->id) continue;
                Vec pj = monitor->get_pos_cur(j);
                Vec vj = monitor->get_v_cur(j);
                double rj = monitor->get_r(j);

                Vec rel_p = p - pj;
                double d0 = len(rel_p);
                if (d0 <= eps) continue;

                double interact_range = r_self + rj + std::max(vmax, len(vj)) * tau * 2.0 + 1e-6;
                if (d0 > interact_range) continue;

                Vec rel_v = v - vj;
                double A = len2(rel_v);
                double tstar = 0.0;
                if (A > eps) {
                    double B = dot(rel_p, rel_v);
                    tstar = -B / A;
                    if (tstar < 0.0) tstar = 0.0;
                    if (tstar > tau) tstar = tau;
                }

                Vec at_min = rel_p + rel_v * tstar;
                double dmin = len(at_min);
                double Rsum = r_self + rj + margin;
                if (dmin <= Rsum) {
                    Vec u = rel_p / d0;
                    double appr = dot(rel_v, u);
                    double target = std::max(0.0, 0.1 * vmax);
                    double delta = target - appr;
                    v = v + u * delta;
                    v = clamp_len(v, vmax);
                }

                if (d0 < Rsum) {
                    Vec u = rel_p / d0;
                    v = v + u * std::min(0.2 * vmax, (Rsum - d0) / std::max(tau, eps));
                    v = clamp_len(v, vmax);
                }
            }
        }
    }

    v = clamp_len(v, vmax);
    return v;
}
