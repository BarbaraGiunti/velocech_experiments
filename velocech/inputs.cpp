void parse_command_line( int argc, char** argv, std::string& input_filename, std::string& output_filename) {

    if( argc < 2 ) print_help_and_exit();

    input_filename = argv[ argc - 2 ];
    output_filename = argv[ argc - 1 ];

}


// Loads the point cloud from given file in ascii format 
// Format: first line is number of points, and each line after that represents a point as a list of coordinates
// Ignores empty lines and lines starting with a '#'.
bool load_ascii(const std::string& filename, std::vector<std::vector<double>>& points){

    std::ifstream input_stream(filename);
    if (!input_stream.is_open()){ return false;}

    std::string cur_line;
    bool dim_read = false;

    while (std::getline(input_stream, cur_line)) {
        // Trim trailing whitespace
        cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);

        // Skip empty lines and comments
        if (cur_line.empty() || cur_line[0] == '#')
            continue;

        std::stringstream ss(cur_line);

        // Read ambient dimension (first valid line)
        if (!dim_read) {
            ss >> ambient_dim;
            if (ambient_dim <= 0) {return false;}
            dim_read = true;
        }
        else {
            std::vector<double> point;
            double value;

            // Read coordinates
            while (ss >> value) {
                point.push_back(value);
            }

            // Validate dimension
            if ((int)point.size() != ambient_dim) {
                std::cerr << "Error: point dimension mismatch\n";
                return false;
            }

        points.push_back(point);
    }

    input_stream.close();
    return dim_read;
}


\\read_successful = load_ascii( input_filename, point_cloud );

\\Usage
std::vector<std::vector<double>> points;

if (load_ascii("data.txt", points)) {
    std::cout << "Loaded " << points.size() << " points of dimension " << dim << "\n";
}
