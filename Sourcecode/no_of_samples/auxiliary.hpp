//
//  auxiliary.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef auxiliary_hpp
#define auxiliary_hpp

#include "data_and_parameters.hpp"

#include <string>

// *** TIME CONVERSION METHODS ***
inline int day_of_hour (int t) { return t / no_of_hours_per_day; }
inline int first_hour_of_day (int d) { return d * no_of_hours_per_day; }
inline int last_hour_of_day (int d) { return (d + 1) * no_of_hours_per_day - 1; }

// error messaging
void error (const char *method, const char *message);

// time measurement
// both functions can be called in a nested manner: tic(); tic(); toc(); toc();
void tic();
double toc();

// methods for reading in parameters from the HPC
bool extract_param_int (const char *param_env[], const char *param_name, int &param_value);
bool extract_param_double (const char *param_env[], const char *param_name, double &param_value);
bool extract_param_string (const char *param_env[], const char *param_name, std::string &param_value);

#endif /* auxiliary_hpp */
