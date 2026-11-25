#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <cstdint>
#include <cstring>

// Usage: program <input_file> [output_file]
// If output_file is omitted, the barcode is written to standard output.

inline void print_help_and_exit(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <input_file> [output_file]\n"
              << "\n"
              << "input_file:  either a binary ALPHAIN1 file (the .bin files in 'Point clouds/*/inputs'),\n"
              << "             or an ascii point cloud: the first line is the ambient dimension d;\n"
              << "             every following line is one point given by d coordinates,\n"
              << "             separated by spaces, tabs or commas.\n"
              << "             Empty lines and lines starting with '#' are ignored.\n"
              << "             The format is detected automatically from the file contents.\n"
              << "output_file: optional. Each line is one H1 bar 'birth^2 death^2'.\n"
              << "             If omitted, the bars are written to standard output.\n";
    std::exit(1);
}

inline void parse_command_line(int argc, char** argv, std::string& input_filename, std::string& output_filename) {
    if (argc < 2 || argc > 3) print_help_and_exit(argv[0]);
    const std::string first = argv[1];
    if (first == "-h" || first == "--help") print_help_and_exit(argv[0]);
    input_filename = first;
    output_filename = (argc == 3) ? argv[2] : "";
}

// Loads a point cloud from an ascii file (format described in print_help_and_exit).
// On success, points holds one vector of length ambient_dim per point.
// On failure, prints the reason to std::cerr and returns false.
inline bool load_ascii(const std::string& filename, std::vector<std::vector<double>>& points, int& ambient_dim) {
    std::ifstream input_stream(filename);
    if (!input_stream.is_open()) {
        std::cerr << "Error: could not open '" << filename << "'\n";
        return false;
    }

    points.clear();
    ambient_dim = 0;
    bool dim_read = false;
    std::string cur_line;
    long long line_no = 0;

    while (std::getline(input_stream, cur_line)) {
        line_no++;
        // Treat commas as whitespace, so csv files also work
        for (char& c : cur_line) if (c == ',') c = ' ';
        // Skip empty lines and comments
        const size_t first = cur_line.find_first_not_of(" \t\r\f\v");
        if (first == std::string::npos || cur_line[first] == '#') continue;

        std::stringstream ss(cur_line);

        // Read ambient dimension (first valid line)
        if (!dim_read) {
            std::string extra;
            if (!(ss >> ambient_dim) || ambient_dim <= 0 || (ss >> extra)) {
                std::cerr << "Error (line " << line_no << "): the first line must contain only the ambient dimension\n";
                return false;
            }
            dim_read = true;
            continue;
        }

        std::vector<double> point;
        point.reserve(ambient_dim);
        double value;
        while (ss >> value) point.push_back(value);
        if (!ss.eof()) {
            std::cerr << "Error (line " << line_no << "): could not read a number\n";
            return false;
        }
        if ((int)point.size() != ambient_dim) {
            std::cerr << "Error (line " << line_no << "): point has " << point.size()
                      << " coordinates, expected " << ambient_dim << "\n";
            return false;
        }
        points.push_back(std::move(point));
    }

    if (!dim_read) {
        std::cerr << "Error: '" << filename << "' contains no data\n";
        return false;
    }
    return true;
}

// Loads a binary ALPHAIN1 point cloud (the .bin files in 'Point clouds/*/inputs').
// Layout, little-endian:
//   bytes  0- 7  the characters "ALPHAIN1"
//   bytes  8-15  int64 number of points n
//   bytes 16-23  int64 ambient dimension d
//   bytes 24-31  int64 reserved (0 in all files so far; ignored)
//   bytes 32-    n*d doubles, point by point: x_1[0..d-1], x_2[0..d-1], ...
// On failure, prints the reason to std::cerr and returns false.
inline bool load_alphain_bin(const std::string& filename, std::vector<std::vector<double>>& points, int& ambient_dim) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "Error: could not open '" << filename << "'\n";
        return false;
    }
    char magic[8];
    int64_t header[3];
    if (!in.read(magic, 8) || std::memcmp(magic, "ALPHAIN1", 8) != 0
        || !in.read(reinterpret_cast<char*>(header), sizeof(header))) {
        std::cerr << "Error: '" << filename << "' is not an ALPHAIN1 file\n";
        return false;
    }
    const int64_t n = header[0], d = header[1];
    in.seekg(0, std::ios::end);
    const int64_t file_size = (int64_t)in.tellg();
    if (n <= 0 || d <= 0 || d > 1000000 || file_size != 32 + 8 * n * d) {
        std::cerr << "Error: '" << filename << "' has a bad header (n = " << n << ", d = " << d
                  << ", file size " << file_size << " bytes, expected " << 32 + 8 * n * d << ")\n";
        return false;
    }
    in.seekg(32);
    std::vector<double> data((size_t)(n * d));
    if (!in.read(reinterpret_cast<char*>(data.data()), (std::streamsize)(8 * n * d))) {
        std::cerr << "Error: could not read the coordinates from '" << filename << "'\n";
        return false;
    }
    ambient_dim = (int)d;
    points.assign((size_t)n, std::vector<double>((size_t)d));
    for (int64_t i = 0; i < n; i++)
        for (int64_t l = 0; l < d; l++)
            points[i][l] = data[(size_t)(i * d + l)];
    return true;
}

// Loads a point cloud in either format: ALPHAIN1 binary if the file starts with "ALPHAIN1", ascii otherwise.
inline bool load_point_cloud(const std::string& filename, std::vector<std::vector<double>>& points, int& ambient_dim) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "Error: could not open '" << filename << "'\n";
        return false;
    }
    char magic[8] = {};
    in.read(magic, 8);
    const bool is_bin = in.gcount() == 8 && std::memcmp(magic, "ALPHAIN1", 8) == 0;
    in.close();
    return is_bin ? load_alphain_bin(filename, points, ambient_dim)
                  : load_ascii(filename, points, ambient_dim);
}

// Writes the barcode, one bar 'birth death' per line, with full double precision.
// Bars of length zero are skipped, and a bar that never died is written with death 'inf'.
// If filename is empty, writes to standard output.
inline bool write_barcode(const std::string& filename, const std::vector<std::pair<double, double>>& barcode) {
    std::ofstream file;
    if (!filename.empty()) {
        file.open(filename);
        if (!file.is_open()) {
            std::cerr << "Error: could not open '" << filename << "' for writing\n";
            return false;
        }
    }
    std::ostream& out = filename.empty() ? std::cout : file;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto& [birth, death] : barcode) {
        if (death == birth) continue;
        if (death < birth) out << birth << " inf\n";   // deathTime is still -1: never killed
        else out << birth << " " << death << "\n";
    }
    return true;
}
