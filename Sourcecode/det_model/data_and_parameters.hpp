//
//  data_and_parameters.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef data_and_parameters_hpp
#define data_and_parameters_hpp

#include <string>
#include <vector>

// *** DIRECTORY MANAGEMENT ***
const bool          working_on_hpc = true;

const std::string   directory_hpc  = "/rds/general/user/wwiesema/ephemeral/energy";
const std::string   directory_home = "/Users/wolframwiesemann/Documents/Work/Publications/Current Papers/Energy Kilian/Sourcecodes/Data";

const std::string   directory = (working_on_hpc) ? directory_hpc : directory_home;

// *** PROBLEM DATA ***
const int           no_of_reservoirs = 3;
const int           no_of_arcs       = 3;
const int           topology_matrix[no_of_reservoirs][no_of_arcs] = { {-1,  0,  0},
                                                                      {+1, -1,  0},
                                                                      { 0, +1, -1} };


const double        ub_generation[no_of_arcs] = {40.68, 41.40, 50.40};              // in 1000s of m^3
const double        ub_pumping[no_of_arcs]    = {28.541, 0.0, 0.0};                 // in 1000s of m^3

const double        generator_efficiency[no_of_arcs]    = {0.00066812, 0.001032401, 0.000540668};       // in MWh/m^3
const double        inverse_pump_efficiency[no_of_arcs] = {0.000935485, 0.0, 0.0};                      // in MWh/m^3


const double        ub_reservoir_level[no_of_reservoirs] = {18500.0, 230.0, 4.0};                       // in 1000s of m^3
const double        lb_reservoir_level[no_of_reservoirs] = {0.0, 0.0, 0.0};                             // in 1000s of m^3

// initial_reservoir_level is the water level at the beginning of the current day
// reservoir_level__beginning_of_year is the water level at the beginning of the year -- used for overall initialisation
//                  as well as for the yearly cycle constraint
extern double       initial_reservoir_level[no_of_reservoirs];
const double        reservoir_level__beginning_of_year[no_of_reservoirs] = {0.75 * 18500.0, 0.75 * 230.0, 0.75 * 4.0};

const double        prob_of_reserve_up_call   = 0.01;
const double        prob_of_reserve_down_call = 0.01;


extern std::vector<std::vector<double> >    spot_price__in_sample;              // [no_of_insample_scenarios][no_of_hours_per_year]
extern std::vector<std::vector<double> >    reserve_up_price__in_sample;        // [no_of_insample_scenarios][no_of_hours_per_year]
extern std::vector<std::vector<double> >    reserve_down_price__in_sample;      // [no_of_insample_scenarios][no_of_hours_per_year]

extern std::vector<double>                  spot_price__out_of_sample;          // [no_of_hours_per_day]
extern std::vector<double>                  reserve_up_price__out_of_sample;    // [no_of_hours_per_day]
extern std::vector<double>                  reserve_down_price__out_of_sample;  // [no_of_hours_per_day]


// out-of-sample inflows start at the current day of the year, looking one day ahead from the current day
// in-sample inflows start at the current day of the year, looking one year ahead from the current day
// the first 24 hours of in-sample inflows [0]...[no_of_hours_per_day - 1] coincide with the out-of-sample inflows
extern std::vector<std::vector<std::vector<double> > >  inflow__in_sample;      // [no_of_insample_scenarios]
                                                                                // [no_of_reservoirs][no_of_hours_per_year]
extern std::vector<std::vector<double> >                inflow__out_of_sample;  // [no_of_reservoirs][no_of_hours_per_day]


extern std::vector<bool> reserve_up_call__out_of_sample;                        // [no_of_hours_per_day]
extern std::vector<bool> reserve_down_call__out_of_sample;                      // [no_of_hours_per_day]


const int       no_of_weeks_per_year = 52;
const int       no_of_days_per_week = 7;
const int       no_of_days_per_year  = no_of_weeks_per_year * no_of_days_per_week;
const int       no_of_hours_per_day  = 24;
const int       no_of_hours_per_year = no_of_days_per_year * no_of_hours_per_day;


const int       no_of_insample_scenarios__data = 50;        // how many in-sample scenarios do we have as data?
extern int      no_of_insample_scenarios__model;            // how many of these in-sample scenarios do we use in the model?

extern double   deterministic__percentage_of_spot__factor;

#endif /* data_and_parameters_hpp */
