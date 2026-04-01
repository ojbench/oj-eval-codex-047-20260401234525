// Controller implementation for Problem 047
// Only defines class Controller. Uses Vec/Monitor from judge headers.

#ifndef SRC_HPP
#define SRC_HPP

#include <cmath>
#include <algorithm>

class Vec;      // provided by math.h (included before this header by judge)
class Monitor;  // provided by monitor.h (included before this header by judge)

#ifndef TIME_INTERVAL
#define TIME_INTERVAL 0.1
#endif

class Controller {
public:
    Controller(const Vec &pos_tar_, double v_max_, double r_, int id_, Monitor *monitor_)
        : pos_cur(), v_cur(), pos_tar(pos_tar_), r(r_), v_max(v_max_), id(id_), monitor(monitor_) {}

    void set_pos_cur(const Vec &p) { pos_cur = p; }
    void set_v_cur(const Vec &v) { v_cur = v; }

    Vec get_v_next();

private:
    // Helpers relying on Vec having public x/y
    static double dot(const Vec &a, const Vec &b) { return a.x * b.x + a.y * b.y; }
    static double len2(const Vec &a) { return a.x * a.x + a.y * a.y; }
    static double len(const Vec &a) { return std::sqrt(len2(a)); }
    static Vec clamp_len(const Vec &v, double maxlen);

private:
    Vec pos_cur;
    Vec v_cur;
    Vec pos_tar;
    double r = 0.0;
    double v_max = 0.0;
    int id = 0;
    Monitor *monitor = nullptr;
};

// Assume Vec supports +, -, *, / with scalars
inline Vec Controller::clamp_len(const Vec &v, double maxlen) {
    double l2 = len2(v);
    if (l2 <= 0.0) return v;
    double l = std::sqrt(l2);
    if (l <= maxlen) return v;
    return v * (maxlen / l);
}

inline Vec Controller::get_v_next() {
    const double eps = 1e-9;
    const double tau = TIME_INTERVAL;
    const double margin = 1e-3;
    const double arrive_thr = std::max(1e-3, 0.5 * std::max(0.0, r));

    Vec to_tar = pos_tar - pos_cur;
    double dist = len(to_tar);
    Vec v_des = v_cur * 0.0;
    if (dist <= arrive_thr) {
        v_des = v_des * 0.0; // stop near target
    } else if (dist > eps) {
        double desired_speed = std::min(v_max, dist / std::max(tau, eps));
        v_des = (to_tar / dist) * desired_speed;
    }

    if (monitor && monitor->get_speeding(id)) {
        v_des = v_des * 0.8;
    }

    Vec v = v_des;
    if (monitor) {
        int N = monitor->get_robot_number();
        for (int pass = 0; pass < 3; ++pass) {
            for (int j = 0; j < N; ++j) {
                if (j == id) continue;
                Vec pj = monitor->get_pos_cur(j);
                Vec vj = monitor->get_v_cur(j);
                double rj = monitor->get_r(j);

                Vec rel_p = pos_cur - pj;
                double d0 = len(rel_p);
                if (d0 <= eps) continue;

                double interact_range = r + rj + std::max(v_max, len(vj)) * tau * 2.0 + 1e-6;
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
                double Rsum = r + rj + margin;
                if (dmin <= Rsum) {
                    Vec u = rel_p / d0;
                    double appr = dot(rel_v, u);
                    double target = std::max(0.0, 0.1 * v_max);
                    double delta = target - appr;
                    v = v + u * delta;
                    // add tangential sidestep to break symmetry
                    Vec u_perp(-u.y, u.x);
                    double sidestep = std::min(0.3 * v_max, (Rsum - dmin) / std::max(tau, eps));
                    double sign = (id < j) ? 1.0 : -1.0;
                    v = v + u_perp * (sign * sidestep);
                    v = clamp_len(v, v_max);
                }

                if (d0 < Rsum) {
                    Vec u = rel_p / d0;
                    v = v + u * std::min(0.2 * v_max, (Rsum - d0) / std::max(tau, eps));
                    v = clamp_len(v, v_max);
                }
            }
        }
    }

    v = clamp_len(v, v_max);
    return v;
}

#endif // SRC_HPP
