// ===========================================================================
// gudhi_alpha.cpp
//
// Alpha filtration of a point cloud computed with GUDHI (Gudhi::alpha_complex::
// Alpha_complex, built on CGAL's d-dimensional Delaunay triangulation, so ANY
// ambient dimension works), with the same interface, file formats and timing
// protocol as the program alphadual -- so that the two can be compared run by
// run and the results read with the same tools (alphadual_io.read_output).
//
// Usage:
//   gudhi_alpha INPUT OUTPUT --r R [--max-dim D] [--precision fast|safe|exact]
//               [--stats FILE.json] [--results-csv FILE.csv]
//   gudhi_alpha --version
//
//   R          radius r_max (squared radius r^2 = GUDHI's max_alpha_square);
//              "inf" for the complete filtration; "half-diameter" for
//              r = diameter/2 (the diameter is then computed inside the timed
//              part, exactly, by comparing all pairs of points).
//   D          largest simplex dimension written to OUTPUT (default 2).
//              GUDHI always builds the full Delaunay triangulation first, so
//              D does not reduce its work; it only selects the output.
//   precision  GUDHI's three modes (default safe, as in GUDHI):
//                fast  = CGAL::Epick_d, approximate filtration values
//                safe  = CGAL::Epeck_d, filtration values accurate to ~1e-10
//                exact = CGAL::Epeck_d, exactly rounded filtration values
//
// Input: .bin (alphadual input format "ALPHAIN1"), or text .csv/.txt/.xyz
// (one point per line, coordinates separated by commas or spaces; a first
// line holding only the number of points is skipped).
// Output: alphadual output format 2 ("ALPHOUT2": simplices up to dimension D
// in filtration order: by value, then dimension, then vertex ids) and a JSON
// statistics file (OUTPUT.json unless --stats).
//
// Timing (as in alphadual): one interval "alpha_time_s" from just before the
// input is opened to the finished ordered filtration (writing excluded); one
// extra clock reading after the input is read splits it into read_time_s and
// compute_time_s. With --r half-diameter, compute_time_s includes the
// diameter (also reported alone as diameter_time_s).
//
// No timer, thread or check runs during the computation. To stop runs that
// take too long, use a time limit from outside, e.g. "timeout 600
// ./gudhi_alpha ..." (run_gudhi_folder.sh --time-limit does that).
// ===========================================================================

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>
#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>
// gudhi/Version.h exists in GUDHI releases (it is generated when GUDHI itself
// is configured, so a raw source checkout may lack it).
#if defined(__has_include)
#if __has_include(<gudhi/Version.h>)
#include <gudhi/Version.h>
#endif
#endif
#ifndef GUDHI_VERSION
#define GUDHI_VERSION "unknown"
#endif

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/resource.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#endif

using std::vector;
using i64 = int64_t;
using i32 = int32_t;

[[noreturn]] static void die(const std::string& msg) {
    std::fprintf(stderr, "gudhi_alpha: %s\n", msg.c_str());
    std::exit(1);
}

static double now_s() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------
// Peak process memory (MB), read after the timed region
// ---------------------------------------------------------------------------
static double peak_memory_mb() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc)))
        return double(pmc.PeakPagefileUsage) / 1048576.0;
    return -1;
#else
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
    return double(ru.ru_maxrss) / 1048576.0;      // bytes on macOS
#else
    return double(ru.ru_maxrss) / 1024.0;         // kilobytes on Linux
#endif
#endif
}

// ---------------------------------------------------------------------------
// Run description (for the results CSV line)
// ---------------------------------------------------------------------------
struct RunInfo {
    std::string input, stats_path, results_csv, precision, r_text;
    i64 N = 0, m = 0;
    int max_dim = 2;
};
static RunInfo g_run;

static std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        o += c;
    }
    return o;
}

static std::string base_name(const std::string& p) {
    const size_t k = p.find_last_of("/\\");
    return k == std::string::npos ? p : p.substr(k + 1);
}

// One line of the results CSV (header written if the file is new/empty).
static void append_results_csv(const std::string& status, double read_s, double compute_s,
                               double total_s, double diam_s, const vector<i64>& counts) {
    if (g_run.results_csv.empty()) return;
    bool need_header = true;
    {
        std::ifstream t(g_run.results_csv);
        need_header = !t.good() || t.peek() == std::ifstream::traits_type::eof();
    }
    FILE* f = std::fopen(g_run.results_csv.c_str(), "a");
    if (!f) return;
    if (need_header) {
        std::fprintf(f, "file,status,n_points,ambient_dim,r,max_dim,precision,read_time_s,"
                        "compute_time_s,alpha_time_s,diameter_time_s,peak_memory_mb");
        for (int k = 0; k <= g_run.max_dim; ++k) std::fprintf(f, ",n_dim%d", k);
        std::fprintf(f, "\n");
    }
    auto num = [&](double v) {
        if (std::isnan(v)) std::fprintf(f, ",");
        else std::fprintf(f, ",%.9f", v);
    };
    std::fprintf(f, "%s,%s,%lld,%lld,%s,%d,%s", base_name(g_run.input).c_str(), status.c_str(),
                 (long long)g_run.N, (long long)g_run.m, g_run.r_text.c_str(), g_run.max_dim,
                 g_run.precision.c_str());
    num(read_s); num(compute_s); num(total_s); num(diam_s);
    std::fprintf(f, ",%.1f", peak_memory_mb());
    for (int k = 0; k <= g_run.max_dim; ++k) {
        if (k < int(counts.size())) std::fprintf(f, ",%lld", (long long)counts[k]);
        else std::fprintf(f, ",");
    }
    std::fprintf(f, "\n");
    std::fclose(f);
}

// ---------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------
static void read_input(const std::string& path, i64& N, i64& m, vector<double>& X) {
    std::string ext = path.size() >= 4 ? path.substr(path.size() - 4) : "";
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".bin") {
        FILE* f = std::fopen(path.c_str(), "rb");
        if (!f) die("cannot open " + path);
        char magic[8];
        i64 hdr[3];
        if (std::fread(magic, 1, 8, f) != 8 || std::memcmp(magic, "ALPHAIN1", 8) != 0)
            die(path + " is not an alphadual input file");
        if (std::fread(hdr, sizeof(i64), 3, f) != 3) die("truncated header in " + path);
        N = hdr[0];
        m = hdr[1];
        X.resize(size_t(N * m));
        if (std::fread(X.data(), sizeof(double), X.size(), f) != X.size()) die("truncated " + path);
        std::fclose(f);
        return;
    }
    std::ifstream in(path);
    if (!in) die("cannot open " + path);
    std::string line;
    vector<vector<double>> rows;
    while (std::getline(in, line)) {
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t' || c == '\r') c = ' ';
        std::istringstream ss(line);
        vector<double> row;
        double v;
        while (ss >> v) row.push_back(v);
        if (!row.empty()) rows.push_back(std::move(row));
    }
    if (rows.size() > 1 && rows[0].size() == 1 && rows[1].size() > 1)
        rows.erase(rows.begin());                    // header line: number of points
    if (rows.size() < 2) die("fewer than 2 points in " + path);
    N = i64(rows.size());
    m = i64(rows[0].size());
    X.resize(size_t(N * m));
    for (i64 i = 0; i < N; ++i) {
        if (i64(rows[i].size()) != m) die("line " + std::to_string(i + 1) + " of " + path +
                                          " has a different number of coordinates");
        std::copy(rows[i].begin(), rows[i].end(), X.begin() + i * m);
    }
}

// Exact diameter: largest distance over all pairs (single-threaded).
static double diameter(const vector<double>& X, i64 N, i64 m) {
    double best = 0.0;
    for (i64 i = 0; i < N; ++i)
        for (i64 j = i + 1; j < N; ++j) {
            double s = 0.0;
            for (i64 d = 0; d < m; ++d) {
                const double t = X[i * m + d] - X[j * m + d];
                s += t * t;
            }
            best = std::max(best, s);
        }
    return std::sqrt(best);
}

// ---------------------------------------------------------------------------
// GUDHI computation -> ordered filtration (value, dim, vertices)
// ---------------------------------------------------------------------------
struct Rec {
    double v;
    int d;
    vector<i32> s;
};

template <class Kernel>
static void gudhi_filtration(const vector<double>& X, i64 N, i64 m, double a1, int max_dim,
                             bool exact, vector<Rec>& out) {
    using Point = typename Kernel::Point_d;
    vector<Point> pts;
    pts.reserve(size_t(N));
    for (i64 i = 0; i < N; ++i) pts.emplace_back(X.begin() + i * m, X.begin() + (i + 1) * m);
    Gudhi::alpha_complex::Alpha_complex<Kernel> ac(pts);   // Delaunay triangulation
    Gudhi::Simplex_tree<> st;
    if (!ac.create_complex(st, a1, exact)) die("GUDHI could not create the alpha complex");
    for (auto sh : st.complex_simplex_range()) {
        const int d = st.dimension(sh);
        if (d > max_dim) continue;
        const double v = st.filtration(sh);
        if (v > a1) continue;
        Rec r;
        r.v = v;
        r.d = d;
        for (auto vh : st.simplex_vertex_range(sh)) r.s.push_back(i32(vh));  // = point index
        std::sort(r.s.begin(), r.s.end());
        out.push_back(std::move(r));
    }
    // Same order as alphadual: value, then dimension, then vertex ids.
    std::sort(out.begin(), out.end(), [](const Rec& p, const Rec& q) {
        if (p.v != q.v) return p.v < q.v;
        if (p.d != q.d) return p.d < q.d;
        return p.s < q.s;
    });
}

static void write_output(const std::string& path, const vector<Rec>& recs, int max_dim) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) die("cannot write " + path);
    const i64 width = max_dim + 1, M = i64(recs.size());
    vector<i64> counts(size_t(width), 0);
    for (const Rec& r : recs) ++counts[r.d];
    std::fwrite("ALPHOUT2", 1, 8, f);
    std::fwrite(&M, sizeof(i64), 1, f);
    std::fwrite(&width, sizeof(i64), 1, f);
    std::fwrite(counts.data(), sizeof(i64), counts.size(), f);
    vector<i32> simp(size_t(M * width), -1);
    vector<int8_t> dims(static_cast<size_t>(M));
    vector<double> vals(static_cast<size_t>(M));
    for (i64 i = 0; i < M; ++i) {
        std::copy(recs[i].s.begin(), recs[i].s.end(), simp.begin() + i * width);
        dims[i] = int8_t(recs[i].d);
        vals[i] = recs[i].v;
    }
    std::fwrite(simp.data(), sizeof(i32), simp.size(), f);
    std::fwrite(dims.data(), 1, dims.size(), f);
    std::fwrite(vals.data(), sizeof(double), vals.size(), f);
    if (std::fclose(f) != 0) die("error writing " + path);
}

#define GA_STR2(x) #x
#define GA_STR(x) GA_STR2(x)
#if defined(__VERSION__)
#define GA_COMPILER __VERSION__
#elif defined(_MSC_FULL_VER)
#define GA_COMPILER "MSVC " GA_STR(_MSC_FULL_VER)
#else
#define GA_COMPILER "unknown"
#endif
#ifndef GA_BUILD_FLAGS
#define GA_BUILD_FLAGS "unknown"
#endif

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::printf("gudhi_alpha 1.0 | GUDHI %s | CGAL %s | compiler %s | flags %s\n",
                    GUDHI_VERSION, GA_STR(CGAL_VERSION), GA_COMPILER, GA_BUILD_FLAGS);
        return 0;
    }
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: gudhi_alpha INPUT OUTPUT --r R|inf|half-diameter [--max-dim D] "
                     "[--precision fast|safe|exact] "
                     "[--stats FILE] [--results-csv FILE]\n");
        return 2;
    }
    g_run.input = argv[1];
    const std::string out_path = argv[2];
    g_run.stats_path = out_path + ".json";
    g_run.precision = "safe";
    double r = NAN;
    bool half_diameter = false;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) die("missing value after " + a);
            return argv[++i];
        };
        if (a == "--r") {
            g_run.r_text = next();
            if (g_run.r_text == "half-diameter") half_diameter = true;
            else r = std::strtod(g_run.r_text.c_str(), nullptr);
        } else if (a == "--max-dim") g_run.max_dim = std::atoi(next().c_str());
        else if (a == "--precision") g_run.precision = next();
        else if (a == "--stats") g_run.stats_path = next();
        else if (a == "--results-csv") g_run.results_csv = next();
        else die("unknown option " + a);
    }
    if (std::isnan(r) && !half_diameter) die("give --r R (a number, inf or half-diameter)");
    if (g_run.max_dim < 0) die("--max-dim must be >= 0");
    if (g_run.precision != "fast" && g_run.precision != "safe" && g_run.precision != "exact")
        die("--precision must be fast, safe or exact");

    // ===================== timed region: read + compute =====================
    const double t0 = now_s();
    vector<double> X;
    read_input(g_run.input, g_run.N, g_run.m, X);
    const double t_read = now_s();
    double diam_s = NAN, rr = r;
    if (half_diameter) {
        rr = diameter(X, g_run.N, g_run.m) / 2.0;
        diam_s = now_s() - t_read;
    }
    const double a1 = std::isinf(rr) ? std::numeric_limits<double>::infinity() : rr * rr;
    vector<Rec> recs;
    if (g_run.precision == "fast")
        gudhi_filtration<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>(X, g_run.N, g_run.m, a1,
                                                                      g_run.max_dim, false, recs);
    else
        gudhi_filtration<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>(
            X, g_run.N, g_run.m, a1, g_run.max_dim, g_run.precision == "exact", recs);
    const double t1 = now_s();
    // ====================== end of timed region =============================

    write_output(out_path, recs, g_run.max_dim);
    const double t2 = now_s();

    vector<i64> counts(size_t(g_run.max_dim + 1), 0);
    for (const Rec& rc : recs) ++counts[rc.d];
    FILE* f = std::fopen(g_run.stats_path.c_str(), "w");
    if (!f) die("cannot write " + g_run.stats_path);
    std::fprintf(f, "{\n  \"implementation\": \"gudhi\",\n  \"status\": \"ok\",\n");
    std::fprintf(f, "  \"precision\": \"%s\",\n  \"input\": \"%s\",\n", g_run.precision.c_str(),
                 json_escape(g_run.input).c_str());
    std::fprintf(f, "  \"N\": %lld,\n  \"m\": %lld,\n  \"max_dim\": %d,\n", (long long)g_run.N,
                 (long long)g_run.m, g_run.max_dim);
    if (std::isinf(rr)) std::fprintf(f, "  \"r\": \"inf\",\n  \"a1\": \"inf\",\n");
    else std::fprintf(f, "  \"r\": %.17g,\n  \"a1\": %.17g,\n", rr, a1);
    std::fprintf(f, "  \"alpha_time_s\": %.9f,\n  \"read_time_s\": %.9f,\n", t1 - t0, t_read - t0);
    std::fprintf(f, "  \"compute_time_s\": %.9f,\n  \"write_time_s\": %.9f,\n", t1 - t_read, t2 - t1);
    if (half_diameter) std::fprintf(f, "  \"diameter_time_s\": %.9f,\n", diam_s);
    // peak memory: read once, after the timed region
    std::fprintf(f, "  \"peak_memory_mb\": %.1f,\n  \"n_simplices\": [", peak_memory_mb());
    for (int k = 0; k <= g_run.max_dim; ++k)
        std::fprintf(f, "%s%lld", k ? ", " : "", (long long)counts[k]);
    std::fprintf(f, "],\n");
    auto cnt = [&](int k) { return (long long)(k <= g_run.max_dim ? counts[k] : 0); };
    std::fprintf(f, "  \"n_vertices\": %lld,\n  \"n_edges\": %lld,\n  \"n_triangles\": %lld,\n",
                 cnt(0), cnt(1), cnt(2));
    if (recs.empty()) std::fprintf(f, "  \"max_value\": null\n}\n");
    else std::fprintf(f, "  \"max_value\": %.17g\n}\n", recs.back().v);
    std::fclose(f);
    append_results_csv("ok", t_read - t0, t1 - t_read, t1 - t0, diam_s, counts);
    return 0;
}
