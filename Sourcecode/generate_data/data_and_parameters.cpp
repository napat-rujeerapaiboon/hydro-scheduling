//
//  data_and_parameters.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "data_and_parameters.hpp"

using namespace std;

// initial_reservoir_level is the water level at the beginning of the current day
double          initial_reservoir_level[no_of_reservoirs];

vector<vector<double> > spot_price__in_sample;          // [no_of_insample_scenarios][no_of_hours_per_year]
vector<vector<double> > reserve_up_price__in_sample;    // [no_of_insample_scenarios][no_of_hours_per_year]
vector<vector<double> > reserve_down_price__in_sample;  // [no_of_insample_scenarios][no_of_hours_per_year]

vector<double>  spot_price__out_of_sample;              // [no_of_hours_per_year]
vector<double>  reserve_up_price__out_of_sample;        // [no_of_hours_per_year]
vector<double>  reserve_down_price__out_of_sample;      // [no_of_hours_per_year]

// out-of-sample inflows start at the current day of the year, looking one day ahead from the current day
// in-sample inflows start at the current day of the year, looking one year ahead from the current day
// the first 24 hours of in-sample inflows [0]...[no_of_hours_per_day - 1] coincide with the out-of-sample inflows
vector<vector<vector<double> > > inflow__in_sample;     // [no_of_insample_scenarios]
                                                        // [no_of_reservoirs][no_of_hours_per_year]
vector<vector<double> >          inflow__out_of_sample; // [no_of_reservoirs][no_of_hours_per_day]


vector<bool>    reserve_up_call__out_of_sample;         // [no_of_hours_per_year]
vector<bool>    reserve_down_call__out_of_sample;       // [no_of_hours_per_year]

int             no_of_insample_scenarios__model = 5;    // how many of these in-sample scenarios do we use in the model?
