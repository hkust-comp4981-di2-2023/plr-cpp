#include <concepts>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>

#ifndef PLR_LIBRARY_H
#define PLR_LIBRARY_H

const double DELTA = 0.005;

// This function is similar to the to_type function
// except "a" will result in "a\0\0\0\0\0\0\0"
// It is a reversed version of to_type
template <typename N>
N stringToNumber(std::string str) {
    const size_t value_size = sizeof(N) / sizeof(char);
    union {
        char buffer[value_size];
        N value {};
    } obj;
    size_t count = value_size-1;
    for (auto i = std::begin(str); i != std::end(str); i++) {
        obj.buffer[count--] = *(i);
//        std::cout << count << std::endl;
        assert(count >= 0);
    }
    return obj.value;
}

// A utility function for encoding any value to string
// Encoding scheme will be memcpy()
template<typename T>
std::string to_string(T value) {
    const size_t value_size = sizeof(value) / sizeof(char);
    union {
        char buffer[value_size];
        T value;
    } obj;
    obj.value = value;
    return std::string{std::begin(obj.buffer), std::end(obj.buffer)};
}

// A utility function for decoding any value from string
// If we have a string "a", it will convert to 97 (reversed order)
// Encoding scheme will be memcpy()
template<typename T>
T to_type(std::string str) {
    const size_t value_size = sizeof(T) / sizeof(char);
    union {
        char buffer[value_size];
        T value {};
    } obj;
    size_t count = 0;
    for (auto i = std::begin(str); i != std::end(str); i++) {
        obj.buffer[count++] = *(i);
    }
//    std::cout << obj.buffer[7] << std::endl;
    return obj.value;
}

// Represent a Point in a 2-D plane
template<typename T>
struct Point {
    static_assert(std::is_floating_point<T>(), "Only floating point is allowed to construct this struct.");
    T x;
    T y;

    Point() = default;

    Point(T _x, T _y) : x(_x), y(_y) {}

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

    Line(T _a1, T _a2) : a1(_a1), a2(_a2) {}

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
    D slope;
    D y; // The intersection pt y
    bool operator==(const Segment<N, D> &another) {
        return (this->x_start == another.x_start) && abs(this->slope - another.slope) < DELTA &&
               abs(this->y - another.y) < DELTA;
    }

    bool operator!=(const Segment<N, D> &another) {
        return !this->operator==(another);
    }

    Segment() = default;

    Segment(const Segment<N, D> &) = default;

    Segment(N _x_start, D _slope, D _y) : x_start(_x_start), slope(_slope), y(_y) {}
};

template<typename N, typename D>
Segment<N, D> Segment<N, D>::NO_VALID_SEGMENT = {0, 0, 0};

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
    GreedyPLR(D _gamma) : state(GREEDY_PLR_STATE::NEED_2_PT), gamma(_gamma) {}

    // Process a point
    // REQUIRED: The PLR Model is not at the finishing state
    void process(Point<double> pt) {
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
                process_(pt);
                break;
            default:
                assert(false); // non-reachable code, suppress warning
        }
    }

    // Finish the PLR Model
    // REQUIRED: Has not been called finish()
    std::vector<Segment<N,D>> finish() {
        assert(state != GREEDY_PLR_STATE::FINISHED);
        switch (state) {
            case GREEDY_PLR_STATE::NEED_2_PT:
                state = GREEDY_PLR_STATE::FINISHED;
                break;
            case GREEDY_PLR_STATE::NEED_1_PT:
                state = GREEDY_PLR_STATE::FINISHED;
                processed_segments.push_back(Segment<N, D>{static_cast<N>(s0.x), 0, s0.y});
                break;
            case GREEDY_PLR_STATE::READY:
                state = GREEDY_PLR_STATE::FINISHED;
                processed_segments.push_back(current_segment());
                break;
            default:
                assert(false); // Unreachable code, suppress warning
        }
        // Unreachable state, no sure whether compiler require this
        return processed_segments;
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
    std::vector<Segment<N,D>> processed_segments;

    void setup_() {
        this->rho_lower = Line<D>(s0.getUpperBound(gamma), s1.getLowerBound(gamma));
        this->rho_upper = Line<D>(s0.getLowerBound(gamma), s1.getUpperBound(gamma));
        this->pt_intersection_ = this->rho_lower.getIntersection(rho_upper);
    }

    Segment<N, D> current_segment() {
        N segment_start = s0.x;
        D avg_slope = (rho_upper.a1 + rho_lower.a1) / 2;
        D intercept = -avg_slope * pt_intersection_.x + pt_intersection_.y;
        return Segment<N, D>{segment_start, avg_slope, intercept};
    }

    void process_(Point<D> pt) {
        if (!(rho_lower.above(pt) && rho_upper.below(pt))) {
            auto prev_segment = current_segment();
            s0 = pt;
            state = GREEDY_PLR_STATE::NEED_1_PT;
            processed_segments.push_back(prev_segment);
        }
        auto s_upper = pt.getUpperBound(gamma);
        auto s_lower = pt.getLowerBound(gamma);

        if (rho_upper.below(s_upper)) {
            rho_upper = Line<D>(pt_intersection_, s_upper);
        }
        if (rho_lower.above(s_lower)) {
            rho_lower = Line<D>(pt_intersection_, s_lower);
        }
    }
};

// A class which represents a trained PLR Model Data
// It can be constructed in two ways
// 1. By converting constructor from gamma (error bound)
// 2. By converting constructor from an encoded string
// REQUIRED: String must be encoded from Encode() function.
template<typename N, typename D>
class PLRDataRep {
public:
    void Decode(const std::string &encoded_str) {
        const size_t elementSize = sizeof(Segment<N, D>);
        size_t count = encoded_str.size() / elementSize;
        size_t ptr = 0;
        size_t sizeN = sizeof(N);
        size_t sizeD = sizeof(D);

        assert(encoded_str.size() % elementSize == sizeD);
        this->gamma_ = to_type<D>(encoded_str.substr(ptr, sizeD));
        ptr += sizeD;
        for (int i = 0; i < count; i++) {
            auto n1 = encoded_str.substr(ptr, sizeN);
            ptr += sizeN;
            auto d1 = encoded_str.substr(ptr, sizeD);
            ptr += sizeD;
            auto d2 = encoded_str.substr(ptr, sizeD);
            ptr += sizeD;
            segments_.push_back(Segment<N, D>(to_type<N>(n1), to_type<D>(d1), to_type<D>(d2)));
        }
    }

    std::string Encode() {
        std::stringstream ss;
        ss << to_string<D>(gamma_);
        for (auto i: segments_) {
            N n1 = i.x_start;
            D d1 = i.slope;
            D d2 = i.y;
            ss << to_string(n1);
            ss << to_string(d1);
            ss << to_string(d2);
        }
        segments_.clear();
        return std::move(ss.str());
    }

    PLRDataRep() = delete;

    PLRDataRep(D gamma) : gamma_(gamma), segments_() {}

    PLRDataRep(D gamma, const std::vector<Segment<N, D>> &another) : gamma_(gamma), segments_(another) {}

    void Add(Segment<N, D> seg) {
        segments_.push_back(seg);
    }


    PLRDataRep(std::string encoded_str) {
        Decode(encoded_str);
    }

    D GetGamma() const {
        return gamma_;
    }

    std::vector<Segment<N, D>> GetSegs() const {
        return segments_;
    }

// Return the range of the possible block
// [lower bound, upper bound] (error-bound included)
// with the key encoded as type N
// [2,1] pair indicates its error (or all [l,r] s.t. r < l) is error or invalid.
    std::pair<N, N> GetValue(N key) {
        if (segments_.empty()) {
            return std::pair<N, N>();
        }
        auto comparator = [](const Segment<N,D>& s1, const Segment<N,D>& s2) {
            return s1.x_start < s2.x_start;
        };
        auto it = std::upper_bound(segments_.begin(), segments_.end(), Segment<N,D>(key,0,0), comparator);
        if (it == segments_.begin()) {
            return std::pair<N,N>(2,1);
        }
        auto res = *(--it);
        auto tar = res.slope * key + res.y;
        return std::pair<N,N>(floor(tar-gamma_),ceil((tar+gamma_)));
    }

private:
    D gamma_;
    std::vector<Segment<N, D>> segments_;
};


#endif //PLR_LIBRARY_H
