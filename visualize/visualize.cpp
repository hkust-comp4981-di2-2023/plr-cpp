//
// Created by Oscar Tse on 23/1/2024.
//
#include "../library.h"
#include <vector>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>
using namespace std;


int main() {
    const int GENERATE_TIMES = 25;
    double block_num = 0;
    double key_num = 1;
    double key_jump;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    auto a = std::normal_distribution<double>(500, 100);
    vector<Point<double>> original {};
    srand(time(nullptr));
    // Generate data
    for (int i = 0; i < GENERATE_TIMES;i++) {
        original.push_back(Point<double>(key_num,block_num));
        block_num++;
        key_jump = a(generator);
        key_num += key_jump;
    }

    vector<Segment<uint64_t,double>> seg {};
    auto plr = GreedyPLR<uint64_t,double>(0.0005);
    for (int i = 0;i<original.size();i++) {
        auto s = plr.process(original[i]);
        if (s != Segment<uint64_t,double>::NO_VALID_SEGMENT) {
            seg.push_back(s);
        }
    }
    auto s = plr.finish();
    if (s != Segment<uint64_t,double>::NO_VALID_SEGMENT) {
        seg.push_back(s);
    }

    ofstream original_data;
    std::filesystem::path cwd = std::filesystem::current_path() / "original_data.csv";
    cout << filesystem::current_path() << endl;
    original_data.open(cwd.string());
    original_data << "Key Num" << "," << "Block Num" << endl;
//    cout << "Key Num" << ", " << "Block Num" << endl;
    for (int i = 0;i < original.size();i++) {
//        cout << original[i].x << ", " << original[i].y << endl;
        original_data << original[i].x << "," << original[i].y << endl;
    }

    original_data.close();

    ofstream plr_data;
    cwd = std::filesystem::current_path() / "plr_data.csv";
    plr_data.open(cwd.string());
    plr_data << "x_start,x_end,slope,y" << endl;
    for (int i = 0; i < seg.size();i++) {
        plr_data << seg[i].x_start << "," <<seg[i].x_end << "," << seg[i].slope << "," << seg[i].y << endl;
    }
    plr_data.close();
    return 0;
}