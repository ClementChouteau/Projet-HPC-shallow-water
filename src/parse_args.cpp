#include <boost/program_options.hpp>
#include <iostream>
#include <shalw.h>

void parse_args(int argc, char **argv) {
  namespace po = boost::program_options;

  po::options_description desc("Options ");

  desc.add_options()
    ("help,h", "Display this message")
    ("size_x,x", po::value<int>()->default_value(256), "Size of the first dimension of the domain")
    ("size_y,y", po::value<int>()->default_value(256), "Size of the second dimension of the domain")
    ("nb_steps,t", po::value<int>()->default_value(1000), "Number of time steps")
    ("dt", po::value<double>()->default_value(300), "dt")
    ("dx", po::value<double>()->default_value(1000), "dx")
    ("dy", po::value<double>()->default_value(1000), "dy")
    ("pcor", po::value<double>()->default_value(1.e-5), "pcor")
    ("grav", po::value<double>()->default_value(0.01), "grav")
    ("dissip", po::value<double>()->default_value(0.00001), "dissip")
    ("alpha", po::value<double>()->default_value(0.15), "alpha")
    ("hmoy", po::value<double>()->default_value(100), "hmoy")
    ("hinit", po::value<double>()->default_value(15), "Height of the initial state")
    ("export", "Export state of hFil")
    ("export-path", po::value<std::string>()->default_value("."), "Path for the export");

  po::variables_map vars;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vars);
  po::notify(vars);

  if (vars.count("help")) {
    std::cout << desc;
    exit(EXIT_SUCCESS);
  }

  size_x = vars["size_x"].as<int>();
  size_y = vars["size_y"].as<int>();
  nb_steps = vars["nb_steps"].as<int>();
  dt = vars["dt"].as<double>();
  dx = vars["dx"].as<double>();
  dy = vars["dy"].as<double>();
  pcor = vars["pcor"].as<double>();
  grav = vars["grav"].as<double>();
  dissip = vars["dissip"].as<double>();
  alpha = vars["alpha"].as<double>();
  hmoy = vars["hmoy"].as<double>();
  height = vars["hinit"].as<double>();
  file_export = false;
  if (vars.count("export"))
    file_export = true;
  export_path = vars["export-path"].as<std::string>();
}
