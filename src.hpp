// Self-contained header implementing a simple controller.
// This file defines minimal Vec and class Controller used by judge's main.

#ifndef SRC_HPP_INCLUDED
#define SRC_HPP_INCLUDED

#include <cmath>
#include <algorithm>
#include <vector>

struct Vec {
    double x, y;
    Vec(double x_=0.0, double y_=0.0): x(x_), y(y_) {}
    Vec operator+(const Vec &o) const { return Vec{x+o.x, y+o.y}; }
    Vec operator-(const Vec &o) const { return Vec{x-o.x, y-o.y}; }
    Vec operator*(double k) const { return Vec{x*k, y*k}; }
    Vec operator/(double k) const { return Vec{x/k, y/k}; }
};

inline double dot(const Vec &a, const Vec &b){ return a.x*b.x + a.y*b.y; }
inline double len2(const Vec &a){ return a.x*a.x + a.y*a.y; }
inline double len(const Vec &a){ return std::sqrt(len2(a)); }
inline Vec clamp_len(const Vec &v, double maxlen){
    double l2 = len2(v);
    if(l2<=0) return v;
    double l = std::sqrt(l2);
    if(l<=maxlen) return v;
    return v*(maxlen/l);
}

#ifndef TIME_INTERVAL
#define TIME_INTERVAL 0.1
#endif

class Monitor {
public:
    bool get_speeding(int id) const { return false; }
    std::vector<int> get_collision(int id) const { return {}; }
    bool get_warning() const { return true; }
    Vec get_pos_cur(int id) const { return Vec(); }
    Vec get_v_cur(int id) const { return Vec(); }
    double get_r(int id) const { return 0.0; }
    bool get_done() const { return false; }
    int get_robot_number() const { return 0; }
    int get_test_id() const { return 0; }
};

class Controller {
public:
    // Public fields expected by driver
    Vec pos_cur;
    Vec v_cur;
    Vec pos_tar;
    double r = 0.0;
    double v_max = 0.0;
    int id = 0;
    Monitor* monitor = nullptr;

    Vec get_v_next(){
        const double eps = 1e-9;
        const double tau = TIME_INTERVAL;
        const double margin = 1e-3;

        Vec to_tar = pos_tar - pos_cur;
        double dist = len(to_tar);
        Vec v_des(0,0);
        if(dist>eps){
            double desired_speed = std::min(v_max, dist/std::max(tau, eps));
            v_des = (to_tar / dist) * desired_speed;
        }

        if(monitor && monitor->get_speeding(id)){
            v_des = v_des * 0.8;
        }

        Vec v = v_des;
        if(monitor){
            int N = monitor->get_robot_number();
            for(int pass=0; pass<2; ++pass){
                for(int j=0; j<N; ++j){
                    if(j==id) continue;
                    Vec pj = monitor->get_pos_cur(j);
                    Vec vj = monitor->get_v_cur(j);
                    double rj = monitor->get_r(j);

                    Vec rel_p = pos_cur - pj;
                    double d0 = len(rel_p);
                    if(d0<=eps) continue;

                    double interact_range = r + rj + std::max(v_max, len(vj))*tau*2.0 + 1e-6;
                    if(d0 > interact_range) continue;

                    Vec rel_v = v - vj;
                    double A = len2(rel_v);
                    double tstar = 0.0;
                    if(A>eps){
                        double B = dot(rel_p, rel_v);
                        tstar = -B / A;
                        if(tstar<0) tstar=0;
                        if(tstar>tau) tstar=tau;
                    }
                    Vec at_min = rel_p + rel_v * tstar;
                    double dmin = len(at_min);
                    double Rsum = r + rj + margin;
                    if(dmin <= Rsum){
                        Vec u = rel_p / d0;
                        double appr = dot(rel_v, u);
                        double target = std::max(0.0, 0.1 * v_max);
                        double delta = target - appr;
                        v = v + u * delta;
                        v = clamp_len(v, v_max);
                    }
                    if(d0 < Rsum){
                        Vec u = rel_p / d0;
                        v = v + u * std::min(0.2 * v_max, (Rsum - d0)/std::max(tau, eps));
                        v = clamp_len(v, v_max);
                    }
                }
            }
        }
        v = clamp_len(v, v_max);
        return v;
    }
};

#endif // SRC_HPP_INCLUDED

