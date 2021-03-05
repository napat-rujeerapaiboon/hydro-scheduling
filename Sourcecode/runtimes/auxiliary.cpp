//
//  auxiliary.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "auxiliary.hpp"

#include <stack>
#include <chrono>
#include <iostream>

using namespace std;
using namespace chrono;

// error messaging
void error (const char *method, const char *message) {
    cerr << endl;
    cerr << "Error [" << method << "]: " << message << endl;
    cerr << "Program aborted." << endl;
    exit (1);
}

// time measurement
stack<time_point<steady_clock> > recorded_times;

void tic() {
    time_point<steady_clock> current_time = steady_clock::now();
    recorded_times.push (current_time);
}

double toc() {
    time_point<steady_clock> last_time = recorded_times.top();
    recorded_times.pop();
    return ((duration<double>)(steady_clock::now() - last_time)).count();
}

// methods for reading in parameters from the HPC
bool extract_param_int (const char *param_env[], const char *param_name, int &param_value) {
    for (int i = 0; param_env[i] != nullptr; ++i) {
        string str = param_env[i];

        if (str.find (param_name) != string::npos) {
            param_value = atoi (str.substr (str.find_last_of ("=") + 1).c_str());
            cout << "PARAMETER SETTING: " << param_name << " = " << param_value << ";" << endl;
            return true;
        }
    }
    
    return false;
}

bool extract_param_double (const char *param_env[], const char *param_name, double &param_value) {
    for (int i = 0; param_env[i] != nullptr; ++i) {
        string str = param_env[i];
        
        if (str.find (param_name) != string::npos) {
            param_value = atof (str.substr (str.find_last_of ("=") + 1).c_str());
            cout << "PARAMETER SETTING: " << param_name << " = " << param_value << ";" << endl;
            return true;
        }
    }
    
    return false;
}

bool extract_param_string (const char *param_env[], const char *param_name, string &param_value) {
    for (int i = 0; param_env[i] != nullptr; ++i) {
        string str = param_env[i];

        if (str.find (param_name) != string::npos) {
            param_value = str.substr (str.find_last_of ("=") + 1).c_str();
            cout << "PARAMETER SETTING: " << param_name << " = " << param_value << ";" << endl;
            return true;
        }
    }
    
    return false;
}
