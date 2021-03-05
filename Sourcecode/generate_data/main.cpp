//
//  main.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 07/07/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "gurobi_c++.h"

#include "results.hpp"
#include "auxiliary.hpp"
#include "input_data_io.hpp"
#include "master_primal.hpp"
#include "master_dual.hpp"
#include "daily_problem.hpp"
#include "data_and_parameters.hpp"

using namespace std;

void optimized_water_values__1_5_10_samples (int curr_day, bool include_reserve, tResultsFile &results) {
    // main loop: solve the problem for {1 sample}, {5 samples}, {10 samples}
    for (int method = 0; method < 3; ++method) {
        // set up initial_reservoir_levels
        if (curr_day == 0) {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
        } else {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
        }

        // calculate water values
        switch (method) {
            case 0:     no_of_insample_scenarios__model = 1;
                        break;
            case 1:     no_of_insample_scenarios__model = 5;
                        break;
            case 2:     no_of_insample_scenarios__model = 10;
                        break;
            default:    error ("main", "undefined case");
                        break;
        }
        
        double objval = solve_primal_master_problem__stochastic (include_reserve, &results.daily_records[curr_day][method].water_values);
        results.daily_records[curr_day][method].master_objval = objval;
        
        // solve problem for the first hour of the day
        vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
        solve_first_hour_of_daily_problem__stochastic (results.daily_records[curr_day][method].water_values, include_reserve, &bid__spot, &bid__reserve_up, &bid__reserve_down);
        
        // solve problem for each hour of the day
        vector<double> generate, pump, spill;
        generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        for (int theta = 0; theta < no_of_hours_per_day; ++theta)
            solve_each_hour_of_daily_problem__stochastic (theta, results.daily_records[curr_day][method].water_values,
                                            bid__spot, bid__reserve_up, bid__reserve_down,
                                            generate, pump, spill);

        // store all results: initial_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r)
            results.daily_records[curr_day][method].initial_reservoir_levels[r] = initial_reservoir_level[r];
        
        // store all results: water_values --> already filled in above

        // store all results: hourly_records
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            results.daily_records[curr_day][method].hourly_records[t].spot_bid = bid__spot[t];
            results.daily_records[curr_day][method].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

            results.daily_records[curr_day][method].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
            if (reserve_up_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = 0.0;

            results.daily_records[curr_day][method].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
            if (reserve_down_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = 0.0;
        }
        
        // store all results: terminal_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r) {
            double final_reservoir_level = initial_reservoir_level[r];
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                final_reservoir_level += inflow__out_of_sample[r][t];
                for (int a = 0; a < no_of_arcs; ++a)
                    final_reservoir_level += topology_matrix[r][a] * (generate[a * no_of_hours_per_day + t] -
                                                                      pump[a * no_of_hours_per_day + t] +
                                                                      spill[a * no_of_hours_per_day + t]);
            }
            results.daily_records[curr_day][method].terminal_reservoir_levels[r] = final_reservoir_level;
        }
    }
}

void optimized_water_values__primal_dual_mixed (int curr_day, bool include_reserve, tResultsFile &results) {
    // main loop: solve the problem for {primal}, {dual}, {mixed}
    for (int method = 0; method < 3; ++method) {
        // set up initial_reservoir_levels
        if (curr_day == 0) {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
        } else {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
        }

        // calculate water values
        switch (method) {
            case 0:     solve_primal_master_problem__stochastic (include_reserve, &results.daily_records[curr_day][0].water_values);
                        break;
            case 1:     solve_dual_master_problem__stochastic (include_reserve, &results.daily_records[curr_day][1].water_values);
                        break;
            case 2:     for (int r = 0; r < no_of_reservoirs; ++r)
                            results.daily_records[curr_day][2].water_values[r] = (results.daily_records[curr_day][0].water_values[r] +
                                                                                  results.daily_records[curr_day][1].water_values[r]) / 2.0;
                        break;
            default:    error ("main", "undefined case");
                        break;
        }
        
        // solve problem for the first hour of the day
        vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
        solve_first_hour_of_daily_problem__stochastic (results.daily_records[curr_day][method].water_values, include_reserve, &bid__spot, &bid__reserve_up, &bid__reserve_down);
        
        // solve problem for each hour of the day
        vector<double> generate, pump, spill;
        generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        for (int theta = 0; theta < no_of_hours_per_day; ++theta)
            solve_each_hour_of_daily_problem__stochastic (theta, results.daily_records[curr_day][method].water_values,
                                            bid__spot, bid__reserve_up, bid__reserve_down,
                                            generate, pump, spill);

        // store all results: initial_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r)
            results.daily_records[curr_day][method].initial_reservoir_levels[r] = initial_reservoir_level[r];
        
        // store all results: water_values --> already filled in above

        // store all results: hourly_records
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            results.daily_records[curr_day][method].hourly_records[t].spot_bid = bid__spot[t];
            results.daily_records[curr_day][method].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

            results.daily_records[curr_day][method].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
            if (reserve_up_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = 0.0;

            results.daily_records[curr_day][method].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
            if (reserve_down_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = 0.0;
        }
        
        // store all results: terminal_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r) {
            double final_reservoir_level = initial_reservoir_level[r];
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                final_reservoir_level += inflow__out_of_sample[r][t];
                for (int a = 0; a < no_of_arcs; ++a)
                    final_reservoir_level += topology_matrix[r][a] * (generate[a * no_of_hours_per_day + t] -
                                                                      pump[a * no_of_hours_per_day + t] +
                                                                      spill[a * no_of_hours_per_day + t]);
            }
            results.daily_records[curr_day][method].terminal_reservoir_levels[r] = final_reservoir_level;
        }
    }
}

// water_values[method][reservoir]
void constant_water_values__primal_dual_mixed (int curr_day, bool include_reserve, const vector<vector<double> > &water_values, tResultsFile &results) {
    // main loop: solve the problem for {primal}, {dual}, {mixed}
    for (int method = 0; method < 3; ++method) {
        // set up initial_reservoir_levels
        if (curr_day == 0) {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
        } else {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
        }

        // copy water values
        results.daily_records[curr_day][method].water_values = water_values[method];
        
        // solve problem for the first hour of the day
        vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
        solve_first_hour_of_daily_problem__stochastic (results.daily_records[curr_day][method].water_values, include_reserve, &bid__spot, &bid__reserve_up, &bid__reserve_down);
        
        // solve problem for each hour of the day
        vector<double> generate, pump, spill;
        generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        for (int theta = 0; theta < no_of_hours_per_day; ++theta)
            solve_each_hour_of_daily_problem__stochastic (theta, results.daily_records[curr_day][method].water_values,
                                            bid__spot, bid__reserve_up, bid__reserve_down,
                                            generate, pump, spill);

        // store all results: initial_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r)
            results.daily_records[curr_day][method].initial_reservoir_levels[r] = initial_reservoir_level[r];
        
        // store all results: water_values --> already filled in above

        // store all results: hourly_records
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            results.daily_records[curr_day][method].hourly_records[t].spot_bid = bid__spot[t];
            results.daily_records[curr_day][method].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

            results.daily_records[curr_day][method].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
            if (reserve_up_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = 0.0;

            results.daily_records[curr_day][method].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
            if (reserve_down_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = 0.0;
        }
        
        // store all results: terminal_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r) {
            double final_reservoir_level = initial_reservoir_level[r];
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                final_reservoir_level += inflow__out_of_sample[r][t];
                for (int a = 0; a < no_of_arcs; ++a)
                    final_reservoir_level += topology_matrix[r][a] * (generate[a * no_of_hours_per_day + t] -
                                                                      pump[a * no_of_hours_per_day + t] +
                                                                      spill[a * no_of_hours_per_day + t]);
            }
            results.daily_records[curr_day][method].terminal_reservoir_levels[r] = final_reservoir_level;
        }
    }
}

void optimized_water_values__deterministic__robust_bounds (int curr_day, bool include_reserve, tResultsFile &results) {
    // main loop: solve the problem for {deterministic master, deterministic sub}, {deterministic master, stochastic sub}, {stochastic master, deterministic sub}
    for (int method = 0; method < 3; ++method) {
        // set up initial_reservoir_levels
        if (curr_day == 0) {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
        } else {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
        }

        // calculate water values
        if ((method == 0) || (method == 1))
            solve_primal_master_problem__deterministic__robust_bounds (include_reserve, &results.daily_records[curr_day][method].water_values);
        else
            solve_primal_master_problem__stochastic (include_reserve, &results.daily_records[curr_day][method].water_values);
        
        // solve problem for the first hour of the day
        double current_water_levels = 0.0, max_water_levels = 0.0;
        for (int r = 0; r < no_of_reservoirs; ++r) {
            current_water_levels += initial_reservoir_level[r];
            max_water_levels += ub_reservoir_level[r];
        }
        bool deterministic_sufficient_water_available = (current_water_levels >= 0.05 * max_water_levels);
        
        vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
        if ((method == 0) || (method == 2))
            solve_first_hour_of_daily_problem__deterministic__robust_bounds (results.daily_records[curr_day][method].water_values,
                                                              (include_reserve && deterministic_sufficient_water_available),
                                                              &bid__spot, &bid__reserve_up, &bid__reserve_down);
        else
            solve_first_hour_of_daily_problem__stochastic (results.daily_records[curr_day][method].water_values,
                                               include_reserve, &bid__spot, &bid__reserve_up, &bid__reserve_down);

        // solve problem for each hour of the day
        vector<double> generate, pump, spill;
        generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
        for (int theta = 0; theta < no_of_hours_per_day; ++theta) {
            if ((method == 0) || (method == 2))
                solve_each_hour_of_daily_problem__deterministic__robust_bounds (theta, results.daily_records[curr_day][method].water_values,
                    bid__spot, bid__reserve_up, bid__reserve_down,
                    generate, pump, spill,
                    &results.daily_records[curr_day][method].hourly_records[theta].missed_reverse_up_commitment,
                    &results.daily_records[curr_day][method].hourly_records[theta].missed_reverse_down_commitment);
            else
                solve_each_hour_of_daily_problem__stochastic (theta, results.daily_records[curr_day][method].water_values,
                    bid__spot, bid__reserve_up, bid__reserve_down,
                    generate, pump, spill);
        }

        // store all results: initial_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r)
            results.daily_records[curr_day][method].initial_reservoir_levels[r] = initial_reservoir_level[r];
        
        // store all results: water_values --> already filled in above

        // store all results: hourly_records
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            results.daily_records[curr_day][method].hourly_records[t].spot_bid = bid__spot[t];
            results.daily_records[curr_day][method].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

            results.daily_records[curr_day][method].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
            if (reserve_up_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = 0.0;

            results.daily_records[curr_day][method].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
            if (reserve_down_call__out_of_sample[t])
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
            else
                results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = 0.0;
        }
        
        // store all results: terminal_reservoir_levels
        for (int r = 0; r < no_of_reservoirs; ++r) {
            double final_reservoir_level = initial_reservoir_level[r];
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                final_reservoir_level += inflow__out_of_sample[r][t];
                for (int a = 0; a < no_of_arcs; ++a)
                    final_reservoir_level += topology_matrix[r][a] * (generate[a * no_of_hours_per_day + t] -
                                                                      pump[a * no_of_hours_per_day + t] +
                                                                      spill[a * no_of_hours_per_day + t]);
            }
            results.daily_records[curr_day][method].terminal_reservoir_levels[r] = final_reservoir_level;
        }
    }
}

bool optimized_water_values__deterministic__percentage_of_spot (int curr_day, bool include_reserve, tResultsFile &results) {
    const int method = 0;
    
    // set up initial_reservoir_levels
    if (curr_day == 0) {
        for (int r = 0; r < no_of_reservoirs; ++r)
            initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
    } else {
        for (int r = 0; r < no_of_reservoirs; ++r)
            initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
    }

    // calculate water values
    solve_primal_master_problem__deterministic__percentage_of_spot (include_reserve, &results.daily_records[curr_day][method].water_values);

    // solve problem for the first hour of the day
    vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
    solve_first_hour_of_daily_problem__deterministic__percentage_of_spot (results.daily_records[curr_day][method].water_values,
                                                      include_reserve,
                                                      &bid__spot, &bid__reserve_up, &bid__reserve_down);

    // solve problem for each hour of the day
    vector<double> generate, pump, spill;
    generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    for (int theta = 0; theta < no_of_hours_per_day; ++theta) {
        double v = solve_each_hour_of_daily_problem__deterministic__percentage_of_spot (theta, results.daily_records[curr_day][method].water_values,
            bid__spot, bid__reserve_up, bid__reserve_down,
            generate, pump, spill);
        if (v < -5000.0) return false;
    }

    // store all results: initial_reservoir_levels
    for (int r = 0; r < no_of_reservoirs; ++r)
        results.daily_records[curr_day][method].initial_reservoir_levels[r] = initial_reservoir_level[r];
    
    // store all results: water_values --> already filled in above

    // store all results: hourly_records
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        results.daily_records[curr_day][method].hourly_records[t].spot_bid = bid__spot[t];
        results.daily_records[curr_day][method].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

        results.daily_records[curr_day][method].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
        if (reserve_up_call__out_of_sample[t])
            results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
        else
            results.daily_records[curr_day][method].hourly_records[t].reserve_up_revenues = 0.0;

        results.daily_records[curr_day][method].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
        if (reserve_down_call__out_of_sample[t])
            results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
        else
            results.daily_records[curr_day][method].hourly_records[t].reserve_down_revenues = 0.0;
    }
    
    // store all results: terminal_reservoir_levels
    for (int r = 0; r < no_of_reservoirs; ++r) {
        double final_reservoir_level = initial_reservoir_level[r];
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            final_reservoir_level += inflow__out_of_sample[r][t];
            for (int a = 0; a < no_of_arcs; ++a)
                final_reservoir_level += topology_matrix[r][a] * (generate[a * no_of_hours_per_day + t] -
                                                                  pump[a * no_of_hours_per_day + t] +
                                                                  spill[a * no_of_hours_per_day + t]);
        }
        results.daily_records[curr_day][method].terminal_reservoir_levels[r] = final_reservoir_level;
    }
    
    return true;
}

int main (int argc, const char *argv[], const char *param_env[]) {
    tic();

    // find out-of-sample number to work on
    int curr_out_of_sample_no = -1;
    
    if (working_on_hpc) {
        if (!extract_param_int (param_env, "PBS_ARRAY_INDEX", curr_out_of_sample_no))
            error ("main", "Could not find PBS_ARRAY_INDEX.");
    } else curr_out_of_sample_no = 0;
    
    generate_data (curr_out_of_sample_no);

    return 0;
    
    // read out results file; terminate if we are done
    tResultsFile results;
    int curr_day = -1;

    stringstream filename;
    filename << directory << "/Results/outer_sample." << curr_out_of_sample_no << ".txt";
    curr_day = results.read (filename.str().c_str());
    
    if (curr_day == no_of_days_per_year) {
        cout << "Entire year calculated -- nothing left to do!" << endl;
        return 0;
    }
    
    // read & verify data
    filename.str (string());
    filename << directory << "/Outer Sample " << curr_out_of_sample_no << "/Day_" << curr_day << ".txt";
    read_data (curr_out_of_sample_no, curr_day);
    
    verify_data();

    // run the model(s)
    //optimized_water_values__1_5_10_samples (curr_day, true, results);
    //optimized_water_values__deterministic__robust_bounds (curr_day, true, results);
    if (!optimized_water_values__deterministic__percentage_of_spot (curr_day, true, results))
        if (!optimized_water_values__deterministic__percentage_of_spot (curr_day, false, results))
            error ("main", "Daily subproblems infeasible even without reserve market involvement");
    
    /*vector<double> water_values__primal;
    water_values__primal.push_back (0.13164186043956);
    water_values__primal.push_back (0.0934426087912088);
    water_values__primal.push_back (0.0380003);

    vector<double> water_values__dual;
    water_values__dual.push_back (0.1027690);
    water_values__dual.push_back (0.0692433);
    water_values__dual.push_back (0.0285435935467033);

    vector<double> water_values__mixed;
    water_values__mixed.push_back (0.117205444230769);
    water_values__mixed.push_back (0.0813430);
    water_values__mixed.push_back (0.0332719620879121);

    vector<vector<double> > water_values;
    water_values.push_back (water_values__primal);
    water_values.push_back (water_values__dual);
    water_values.push_back (water_values__mixed);

    constant_water_values (curr_day, water_values, results);*/
    
    // write results files
    filename.str (string());
    filename << directory << "/Results/outer_sample." << curr_out_of_sample_no << ".txt";
    cout << "writing to \"" << filename.str().c_str() << "\"." << endl;
    results.write (filename.str().c_str(), curr_day);
    
    results.write_tables (curr_day);
    
    // that's it!
    cout << "*** OVERALL TIME: " << toc() << " SECONDS. ***" << endl;
    return 0;
}

