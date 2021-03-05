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

void optimized_water_values (int curr_day, bool include_reserve, tResultsFile &results) {
    // set up initial_reservoir_levels
    if (curr_day == 0) {
        for (int r = 0; r < no_of_reservoirs; ++r)
            initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
    } else {
        for (int r = 0; r < no_of_reservoirs; ++r)
            initial_reservoir_level[r] = results.daily_records[curr_day - 1].terminal_reservoir_levels[r];
    }

    // calculate water values
    double objval = solve_primal_master_problem__stochastic (include_reserve, &results.daily_records[curr_day].water_values);
    results.daily_records[curr_day].master_objval = objval;
    
    // solve problem for the first hour of the day
    vector<double> bid__spot, bid__reserve_up, bid__reserve_down;
    solve_first_hour_of_daily_problem__stochastic (results.daily_records[curr_day].water_values, include_reserve, &bid__spot, &bid__reserve_up, &bid__reserve_down);
    
    // solve problem for each hour of the day
    vector<double> generate, pump, spill;
    generate.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    pump.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    spill.resize (no_of_arcs * no_of_hours_per_day, 0.0);
    for (int theta = 0; theta < no_of_hours_per_day; ++theta)
        solve_each_hour_of_daily_problem__stochastic (theta, results.daily_records[curr_day].water_values,
                                        bid__spot, bid__reserve_up, bid__reserve_down,
                                        generate, pump, spill);

    // store all results: initial_reservoir_levels
    for (int r = 0; r < no_of_reservoirs; ++r)
        results.daily_records[curr_day].initial_reservoir_levels[r] = initial_reservoir_level[r];
    
    // store all results: water_values --> already filled in above

    // store all results: hourly_records
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        results.daily_records[curr_day].hourly_records[t].spot_bid = bid__spot[t];
        results.daily_records[curr_day].hourly_records[t].spot_revenues = spot_price__out_of_sample[t] * bid__spot[t];

        results.daily_records[curr_day].hourly_records[t].reserve_up_bid = bid__reserve_up[t];
        if (reserve_up_call__out_of_sample[t])
            results.daily_records[curr_day].hourly_records[t].reserve_up_revenues = reserve_up_price__out_of_sample[t] * bid__reserve_up[t];
        else
            results.daily_records[curr_day].hourly_records[t].reserve_up_revenues = 0.0;

        results.daily_records[curr_day].hourly_records[t].reserve_down_bid = bid__reserve_down[t];
        if (reserve_down_call__out_of_sample[t])
            results.daily_records[curr_day].hourly_records[t].reserve_down_revenues = reserve_down_price__out_of_sample[t] * bid__reserve_down[t];
        else
            results.daily_records[curr_day].hourly_records[t].reserve_down_revenues = 0.0;
        
        results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
        results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
        results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
        for (int a = 0; a < no_of_arcs; ++a) {
            results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
            results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
            results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
        }
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
        results.daily_records[curr_day].terminal_reservoir_levels[r] = final_reservoir_level;
    }
}

#if 0

void optimized_water_values__N_samples (int curr_day, bool include_reserve, tResultsFile &results) {
    // main loop: solve the problem for {1 sample}, {5 samples}, {10 samples}
    int method = 0;
    {
        // set up initial_reservoir_levels
        if (curr_day == 0) {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = reservoir_level__beginning_of_year[r];
        } else {
            for (int r = 0; r < no_of_reservoirs; ++r)
                initial_reservoir_level[r] = results.daily_records[curr_day - 1][method].terminal_reservoir_levels[r];
        }

        // calculate water values
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

            results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
            for (int a = 0; a < no_of_arcs; ++a) {
                results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
            }
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

            results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
            for (int a = 0; a < no_of_arcs; ++a) {
                results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
            }
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

            results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
            for (int a = 0; a < no_of_arcs; ++a) {
                results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
            }
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

            results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
            results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
            for (int a = 0; a < no_of_arcs; ++a) {
                results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
                results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
            }
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

        results.daily_records[curr_day].hourly_records[t].generate.resize (no_of_arcs, 0.0);
        results.daily_records[curr_day].hourly_records[t].pump.resize (no_of_arcs, 0.0);
        results.daily_records[curr_day].hourly_records[t].spill.resize (no_of_arcs, 0.0);
        for (int a = 0; a < no_of_arcs; ++a) {
            results.daily_records[curr_day].hourly_records[t].generate[a] = generate[a * no_of_hours_per_day + t];
            results.daily_records[curr_day].hourly_records[t].pump[a] = pump[a * no_of_hours_per_day + t];
            results.daily_records[curr_day].hourly_records[t].spill[a] = spill[a * no_of_hours_per_day + t];
        }
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

#endif

int main (int argc, const char *argv[], const char *param_env[]) {
    tic();

    // find out-of-sample number to work on
    int curr_out_of_sample_no = -1;
    
    if (working_on_hpc) {
        if (!extract_param_int (param_env, "PBS_ARRAY_INDEX", curr_out_of_sample_no))
            error ("main", "Could not find PBS_ARRAY_INDEX.");
    } else curr_out_of_sample_no = 0;
    
    // read #-of-samples to work on
    if (argc != 2)
        error ("main", "Expected single command line argument with number of samples to operate on.");
    
    no_of_insample_scenarios__model = atoi (argv[1]);
    if ((no_of_insample_scenarios__model < 1) || (no_of_insample_scenarios__model > no_of_insample_scenarios__data))
        error ("main", "Incorrect number of samples to operate on.");

    // read out results file; terminate if we are done
    tResultsFile results;
    int curr_day = -1;

    stringstream filename;
    filename << directory << "/Results/noreserve_" << no_of_insample_scenarios__model << "_" << curr_out_of_sample_no << ".txt";
    curr_day = results.read (filename.str().c_str());
    
    if (curr_day == no_of_days_per_year) {
        cout << "Entire year calculated -- nothing left to do!" << endl;
        return 0;
    }
    
    // read & verify data
    string tmp_dir;
    if (!extract_param_string (param_env, "TMPDIR", tmp_dir))
        error ("main", "Could not find PBS_ARRAY_INDEX.");
    
    stringstream command;
    command << "tar -xf \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day_" << curr_day << ".tar.gz\" -C $TMPDIR";
    cout << "*** " << command.str().c_str() << endl;
    system (command.str().c_str());
    
    stringstream curr_dir;
    curr_dir << tmp_dir.c_str() << "/rds/general/user/wwiesema/ephemeral/energy/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day;
    read_data (curr_dir.str().c_str());
    verify_data();

    command.str (string());
    command << "rm -rf \"$TMPDIR/rds\"";
    system (command.str().c_str());

    // run the model
    optimized_water_values (curr_day, false, results);
    
    // write results files
    filename.str (string());
    filename << directory << "/Results/noreserve_" << no_of_insample_scenarios__model << "_" << curr_out_of_sample_no << ".txt";
    cout << "writing to \"" << filename.str().c_str() << "\"." << endl;
    results.write (filename.str().c_str(), curr_day);
    
    //results.write_tables (curr_day);
    
    // that's it!
    cout << "*** OVERALL TIME: " << toc() << " SECONDS. ***" << endl;
    return 0;
}

