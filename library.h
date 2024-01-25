#include <concepts>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#ifndef PLR_LIBRARY_H
#define PLR_LIBRARY_H

const double DELTA = 0.005;

template<typename T>
std::string to_string(T value) {
    const size_t value_size = sizeof(value) / sizeof(char);
    union {
        char buffer[value_size];
        T value;
    } obj;
    obj.value = value;
    return std::string {std::begin(obj.buffer), std::end(obj.buffer)};
}

template<typename T>
T to_type(std::string str) {
    const size_t value_size = sizeof(T) / sizeof(char);
    union {3
        char buffer[value_size];
        T value;
    } obj;
    size_t count = 0;
    for (auto i = std::begin(str);i != std::end(str);i++) {
        obj.buffer[count++] = *(i);
    }
    return obj.value;
}

// Represent a Point in a 2-D plane
template<typename T>
struct Point {
    static_assert(std::is_floating_point<T>(), "Only floating point is allowed to construct this struct.");
    T x;
    T y;

    Point() = default;

    Point(T x, T y) : x(x), y(y) {}

    Point<T> getUpperBound(T gamma) {
        return Point(this->x, this->y + gamma);
    }

    Point<T> getLowerBound(T gamma) {
        return Point(this->x, this->y - gamma);
    }
};

// Represent a line in 2D Plane
template<typename T>
struct Line {
    static_assert(std::is_floating_point<T>(), "Only floating point is allowed to construct this struct.");

    T a1; // The factor for x (or slope)
    T a2; // The factor for constant term

    Line() = default;

    Line(T a1, T a2) : a1(a1), a2(a2) {}

    Line(Point<T> a, Point<T> b) {
        this->a1 = (b.y - a.y) / (b.x - a.x);
        this->a2 = -this->a1 * b.x + b.y;
    }

    // Get the intersection pt within two lines
    // REQUIRED: Two lines should not have the same slope
    Point<T> getIntersection(Line<T> anotherLine) {
        assert(!(this->a1 == anotherLine.a1));
        return Point<T>((anotherLine.a2 - this->a2) / (this->a1 - anotherLine.a1),
                        (this->a1 * anotherLine.a2 - anotherLine.a1 * this->a2) / (this->a1 - anotherLine.a1));
    }

    // state whether the point is above the current line
    bool above(Point<T> p) {
        return p.y > this->a1 * p.x + this->a2;
    }

    // state whether the point is below the current line
    bool below(Point<T> p) {
        return p.y < this->a1 * p.x + this->a2;
    }
};

// Represent a segment in PLR model
template<typename N, typename D>
struct __attribute__((packed)) Segment {
    static Segment<N, D> NO_VALID_SEGMENT;
    static_assert(std::is_floating_point<D>(), "Floating point should be placed in second placement,");
    static_assert(std::is_integral<N>(), "Integer should be placed in first placement.");

    N x_start; // The intersection pt x (since we can translate this pt exactly)
    N x_end;
    D slope;
    D y; // The intersection pt y
    bool operator==(const Segment<N, D> &another) {
        return (this->x_start == another.x_start) && abs(this->slope - another.slope) < DELTA &&
               abs(this->y - another.y) < DELTA && (this->x_end == another.x_end);
    }

    bool operator!=(const Segment<N, D> &another) {
        return !this->operator==(another);
    }

    Segment() = default;

    Segment(const Segment<N, D> &) = default;

    Segment(N x_start, N x_end, D slope, D y) : x_start(x_start), x_end(x_end), slope(slope), y(y) {}
};

template<typename N, typename D>
Segment<N, D> Segment<N, D>::NO_VALID_SEGMENT = {0, 0, 0, 0};

enum GREEDY_PLR_STATE {
    NEED_2_PT = 0,
    NEED_1_PT,
    READY,
    FINISHED
};

// Greedy PLR Model
template<typename N, typename D>
class GreedyPLR {
    static_assert(std::is_floating_point<D>(), "Floating point should be placed in second placement,");
    static_assert(std::is_integral<N>(), "Integer should be placed in first placement.");
public:
    GreedyPLR(D gamma) : state(GREEDY_PLR_STATE::NEED_2_PT), gamma(gamma) {}

    // Process a point
    // REQUIRED: The PLR Model is not at the finishing state
    Segment<N, D> process(Point<double> pt) {
        assert(state != GREEDY_PLR_STATE::FINISHED);
        last_pt = pt;
        switch (state) {
            case GREEDY_PLR_STATE::NEED_2_PT:
                s0 = pt;
                state = GREEDY_PLR_STATE::NEED_1_PT;
                break;
            case GREEDY_PLR_STATE::NEED_1_PT:
                s1 = pt;
                setup_();
                state = GREEDY_PLR_STATE::READY;
                break;
            case GREEDY_PLR_STATE::READY:
                return process_(pt);
            default:
                assert(false); // non-reachable code, suppress warning
        }
        return Segment<N, D>::NO_VALID_SEGMENT;
    }

    // Finish the PLR Model
    // REQUIRED: Has not been called finish()
    Segment<N, D> finish() {
        assert(state != GREEDY_PLR_STATE::FINISHED);
        switch (state) {
            case GREEDY_PLR_STATE::NEED_2_PT:
                state = GREEDY_PLR_STATE::FINISHED;
                return Segment<N, D>::NO_VALID_SEGMENT;
            case GREEDY_PLR_STATE::NEED_1_PT:
                state = GREEDY_PLR_STATE::FINISHED;
                return Segment<N, D>{static_cast<N>(s0.x), static_cast<N>(s0.x) + 1, 0, s0.y};
            case GREEDY_PLR_STATE::READY:
                state = GREEDY_PLR_STATE::FINISHED;
                return current_segment();
            default:
                assert(false); // Unreachable code, suppress warning
        }
        // Unreachable state, no sure whether compiler require this
        return Segment<N, D>::NO_VALID_SEGMENT;
    }

private:
    GREEDY_PLR_STATE state;
    D gamma;
    Point<D> last_pt;
    Point<D> s0;
    Point<D> s1;
    Point<D> pt_intersection_;
    Line<D> rho_lower;
    Line<D> rho_upper;

    void setup_() {
        this->rho_lower = Line<D>(s0.getUpperBound(gamma), s1.getLowerBound(gamma));
        this->rho_upper = Line<D>(s0.getLowerBound(gamma), s1.getUpperBound(gamma));
        this->pt_intersection_ = this->rho_lower.getIntersection(rho_upper);
    }

    Segment<N, D> current_segment() {
        N segment_start = s0.x;
        N segment_stop = last_pt.x;
        D avg_slope = (rho_upper.a1 + rho_lower.a1) / 2;
        D intercept = -avg_slope * pt_intersection_.x + pt_intersection_.y;
        return Segment<N, D>{segment_start, segment_stop, avg_slope, intercept};
    }

    Segment<N, D> process_(Point<D> pt) {
        if (!(rho_lower.above(pt) && rho_upper.below(pt))) {
            auto prev_segment = current_segment();
            s0 = pt;
            state = GREEDY_PLR_STATE::NEED_1_PT;
            return prev_segment;
        }
        auto s_upper = pt.getUpperBound(gamma);
        auto s_lower = pt.getLowerBound(gamma);

        if (rho_upper.below(s_upper)) {
            rho_upper = Line<D>(pt_intersection_, s_upper);
        }
        if (rho_lower.above(s_lower)) {
            rho_lower = Line<D>(pt_intersection_, s_lower);
        }

        return Segment<N, D>::NO_VALID_SEGMENT;
    }
};


template<typename N, typename D>
class PLRDataRep {
public:
    void Decode(std::string encoded_str) {
        const size_t elementSize = sizeof(Segment<N, D>);
        size_t count = encoded_str.size() / elementSize;
        size_t ptr = 0;
        size_t sizeN = sizeof(N);
        size_t sizeD = sizeof(D);

        assert(encoded_str.size() % elementSize == 0);
        for (int i = 0; i < count; i++) {
            auto n1 = encoded_str.substr(ptr, sizeN);
            ptr += sizeN;
            auto n2 = encoded_str.substr(ptr, sizeN);
            ptr += sizeN;
            auto d1 = encoded_str.substr(ptr, sizeD);
            ptr += sizeD;
            auto d2 = encoded_str.substr(ptr, sizeD);
            ptr += sizeD;
            segments_.push_back(Segment<N, D>(to_type<N>(n1), to_type<N>(n2), to_type<D>(d1), to_type<D>(d2)));
        }
    }

    std::string Encode() {
        std::string res{};
        for (auto i: segments_) {
            N n1 = i.x_start;
            N n2 = i.x_end;
            D d1 = i.slope;
            D d2 = i.y;
            res += to_string(n1);
            res += to_string(n2);
            res += to_string(d1);
            res += to_string(d2);
        }
        segments_.clear();
        return res;
    }

    PLRDataRep() = default;

    void Add(Segment<N,D> seg) {
        segments_.push_back(seg);
    }


    PLRDataRep(std::string encoded_str) {
        Decode(encoded_str);
    }

    const std::vector<Segment<N,D>> getSegs() const {
        return segments_;
    }

private:
    std::vector<Segment<N, D>> segments_;

};


#endif //PLR_LIBRARY_H
