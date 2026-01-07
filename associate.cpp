#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

struct FileEntry {
    double timestamp;
    string filename;
};

// 读取文件函数
vector<FileEntry> readFile(const string& path) {
    vector<FileEntry> data;
    ifstream file(path);
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        string timeStr, name;
        ss >> timeStr >> name;
        data.push_back({stod(timeStr), name});
    }
    return data;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: ./associate <rgb_file> <depth_file>" << endl;
        return -1;
    }

    auto rgb = readFile(argv[1]);
    auto depth = readFile(argv[2]);

    // 简单的最近邻匹配逻辑 (类似 Python 脚本)
    double max_diff = 0.02;
    
    // 设置输出精度，保证时间戳不丢失精度
    cout << fixed << setprecision(6);

    for (const auto& r : rgb) {
        double best_diff = max_diff;
        string best_depth_file = "";
        double best_depth_time = 0.0;
        bool found = false;

        // 在深度图中寻找时间最接近的一帧
        // 注意：这里为了代码简单，用了暴力搜索。
        // 因为TUM数据集不大，速度依然很快(毫秒级)。
        for (const auto& d : depth) {
            double diff = abs(r.timestamp - d.timestamp);
            if (diff < best_diff) {
                best_diff = diff;
                best_depth_file = d.filename;
                best_depth_time = d.timestamp;
                found = true;
            }
        }

        if (found) {
            cout << r.timestamp << " " << r.filename << " " 
                 << best_depth_time << " " << best_depth_file << endl;
        }
    }

    return 0;
}
