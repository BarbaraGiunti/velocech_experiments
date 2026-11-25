#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epeck_d.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
// Types
// -----------------------------------------------------------------------------

using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
using Point = Kernel::Point_d;
using Point_vector = std::vector<Point>;

using Alpha_complex = Gudhi::alpha_complex::Alpha_complex<Kernel>;
using Simplex_tree = Gudhi::Simplex_tree<>;

using Clock = std::chrono::steady_clock;


// -----------------------------------------------------------------------------
// Read point cloud
//
// Input format:
//
//     N D
//     x11 x12 ... x1D
//     x21 x22 ... x2D
//     ...
//     xN1 xN2 ... xND
//
// Example:
//
//     7 2
//     1.0 1.0
//     7.0 0.0
//     4.0 6.0
//     ...
// -----------------------------------------------------------------------------

Point_vector read_point_cloud(
    const std::string& filename,
    std::size_t& num_points,
    std::size_t& dimension)
{
    std::ifstream input(filename);

    if (!input) {
        throw std::runtime_error(
            "Cannot open input file: " + filename
        );
    }

    input >> num_points >> dimension;

    if (!input) {
        throw std::runtime_error(
            "Invalid header. Expected: N dimension"
        );
    }

    if (num_points == 0) {
        throw std::runtime_error(
            "The point cloud contains zero points."
        );
    }

    if (dimension == 0) {
        throw std::runtime_error(
            "The ambient dimension must be greater than zero."
        );
    }

    Point_vector points;
    points.reserve(num_points);

    for (std::size_t i = 0; i < num_points; ++i) {

        std::vector<double> coordinates(dimension);

        for (std::size_t j = 0; j < dimension; ++j) {

            if (!(input >> coordinates[j])) {

                std::ostringstream message;

                message
                    << "Unexpected end of input while reading "
                    << "point " << i
                    << ", coordinate " << j;

                throw std::runtime_error(message.str());
            }
        }

        points.emplace_back(
            coordinates.begin(),
            coordinates.end()
        );
    }

    return points;
}


// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main(int argc, char** argv)
{
    try {

        if (argc != 2) {

            std::cerr
                << "Usage:\n"
                << "  " << argv[0] << " INPUT_FILE\n";

            return EXIT_FAILURE;
        }

        const std::string input_file = argv[1];

        // ---------------------------------------------------------------------
        // Read input.
        //
        // This happens OUTSIDE the benchmark timer.
        // ---------------------------------------------------------------------

        std::size_t num_points;
        std::size_t dimension;

        Point_vector points =
            read_point_cloud(
                input_file,
                num_points,
                dimension
            );


        // ---------------------------------------------------------------------
        // Full alpha-complex construction.
        //
        // The timer includes:
        //
        //   1. Alpha_complex construction
        //   2. Delaunay triangulation construction performed by GUDHI
        //   3. create_complex()
        //   4. alpha filtration construction
        //
        // Input-file parsing is NOT included.
        // ---------------------------------------------------------------------

        auto start = Clock::now();

        Alpha_complex alpha_complex(points);

        Simplex_tree simplex_tree;

        const bool success =
            alpha_complex.create_complex(simplex_tree);

        auto stop = Clock::now();


        // ---------------------------------------------------------------------
        // Check result.
        // ---------------------------------------------------------------------

        if (!success) {

            std::cerr
                << "ERROR: GUDHI failed to construct "
                << "the alpha complex.\n";

            return EXIT_FAILURE;
        }


        // ---------------------------------------------------------------------
        // Timing.
        // ---------------------------------------------------------------------

        const double elapsed =
            std::chrono::duration<double>(
                stop - start
            ).count();


        // ---------------------------------------------------------------------
        // Machine-readable output.
        //
        // Format:
        //
        //   mode,num_points,dimension,time_seconds,num_simplices
        // ---------------------------------------------------------------------

        std::cout
            << std::setprecision(17)
            << "full,"
            << num_points
            << ","
            << dimension
            << ","
            << elapsed
            << ","
            << simplex_tree.num_simplices()
            << '\n';


        return EXIT_SUCCESS;
    }
    catch (const std::exception& e) {

        std::cerr
            << "ERROR: "
            << e.what()
            << '\n';

        return EXIT_FAILURE;
    }
}