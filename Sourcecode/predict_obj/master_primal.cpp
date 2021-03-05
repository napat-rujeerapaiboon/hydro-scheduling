//
//  master_primal.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "master_primal.hpp"

#include <fstream>

#include "gurobi_c++.h"

#include "auxiliary.hpp"
#include "data_and_parameters.hpp"

using namespace std;

double solve_primal_master_problem__stochastic (bool include_reserve, vector<double> *water_values, const char *filename) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    tic();
    cout << "Setting up decision variables..." << flush;
    
    GRBVar *bid__spot         = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);
    GRBVar *bid__reserve_up   = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);
    GRBVar *bid__reserve_down = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);

    GRBVar *generate = model.addVars (no_of_insample_scenarios__model * no_of_arcs * no_of_hours_per_year);
    GRBVar *pump     = model.addVars (no_of_insample_scenarios__model * no_of_arcs * no_of_hours_per_year);
    GRBVar *spill    = model.addVars (no_of_insample_scenarios__model * no_of_arcs * no_of_hours_per_year);

    GRBVar *initial_water_level = model.addVars (no_of_reservoirs);             // to obtain shadow prices
    GRBVar *water_level = model.addVars (no_of_insample_scenarios__model * no_of_reservoirs * no_of_days_per_year);
    GRBVar *water_level__LDR__intercept = model.addVars (no_of_reservoirs * no_of_days_per_year);
    GRBVar *water_level__LDR__slope_cum_inflow_curr_day = model.addVars (no_of_reservoirs * no_of_days_per_year * no_of_reservoirs);
    GRBVar *water_level__LDR__slope_cum_inflow_past = model.addVars (no_of_reservoirs * no_of_days_per_year * no_of_reservoirs);
    GRBVar *water_level__LDR__slope_avg_spot_curr_day = model.addVars (no_of_reservoirs * no_of_days_per_year);
    
    cout << "done: " << toc() << " seconds." << endl;
    
    GRBLinExpr expr;
    
    // set up objective function
    tic();
    cout << "Setting up objective function..." << flush;
    
    expr = 0;
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            expr += (spot_price__in_sample[s][t] * bid__spot[s * no_of_hours_per_year + t] +
                     prob_of_reserve_up_call * reserve_up_price__in_sample[s][t] * bid__reserve_up[s * no_of_hours_per_year + t] +
                     prob_of_reserve_down_call * reserve_down_price__in_sample[s][t] * bid__reserve_down[s * no_of_hours_per_year + t])
                    / (double)(no_of_insample_scenarios__model);
    model.setObjective (expr, GRB_MAXIMIZE);
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up variable bounds
    tic();
    cout << "Setting up variable bounds..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            bid__spot[s * no_of_hours_per_year + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int s = 0; s < no_of_insample_scenarios__model; ++s)
            for (int t = 0; t < no_of_hours_per_year; ++t) {
                bid__reserve_up[s * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, 0.0);
                bid__reserve_down[s * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, 0.0);
            }
    }

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int a = 0; a < no_of_arcs; ++a) {
                generate[(s * no_of_arcs + a) * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
                pump[(s * no_of_arcs + a) * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
            }
    
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int d = 0; d < no_of_days_per_year; ++d) {
            for (int s = 0; s < no_of_insample_scenarios__model; ++s) {
                water_level[(s * no_of_reservoirs + r) * no_of_days_per_year + d].set (GRB_DoubleAttr_LB, lb_reservoir_level[r]);
                water_level[(s * no_of_reservoirs + r) * no_of_days_per_year + d].set (GRB_DoubleAttr_UB, ub_reservoir_level[r]);
            }

            water_level__LDR__intercept[r * no_of_days_per_year + d].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            for (int rr = 0; rr < no_of_reservoirs; ++rr) {
                water_level__LDR__slope_cum_inflow_curr_day[(r * no_of_days_per_year + d) * no_of_reservoirs + rr].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                water_level__LDR__slope_cum_inflow_past[(r * no_of_days_per_year + d) * no_of_reservoirs + rr].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            }
            water_level__LDR__slope_avg_spot_curr_day[r * no_of_days_per_year + d].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
        }
        
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: initial water level (for shadow prices)
    tic();
    cout << "Setting up initial water level constraints..." << flush;

    vector<GRBConstr> initial_water_level_constraints;
    for (int r = 0; r < no_of_reservoirs; ++r)
        initial_water_level_constraints.push_back (model.addConstr (initial_water_level[r] == initial_reservoir_level[r]));
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: robust market fulfillment
    tic();
    cout << "Setting up robust market fulfillment constraints..." << flush;

    double maximum_pumping_quantity = 0.0;
    for (int a = 0; a < no_of_arcs; ++a)
        maximum_pumping_quantity += inverse_pump_efficiency[a] * ub_pumping[a];
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            expr = 0;
            for (int a = 0; a < no_of_arcs; ++a)
                expr += generator_efficiency[a] * generate[(s * no_of_arcs + a) * no_of_hours_per_year + t] -
                        inverse_pump_efficiency[a] * pump[(s * no_of_arcs + a) * no_of_hours_per_year + t];
            model.addConstr (bid__spot[s * no_of_hours_per_year + t] +
                             bid__reserve_up[s * no_of_hours_per_year + t] == expr);
            
            model.addConstr (bid__spot[s * no_of_hours_per_year + t] -
                             bid__reserve_down[s * no_of_hours_per_year + t] +
                             maximum_pumping_quantity >= 0.0);
        }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: hourly water level dynamics
    tic();
    cout << "Setting up hourly water level dynamics constraints..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int r = 0; r < no_of_reservoirs; ++r)
            for (int t = 0; t < no_of_hours_per_year; ++t) {
                if (day_of_hour (t) == 0)
                    expr = initial_water_level[r];
                else
                    expr = water_level[(s * no_of_reservoirs + r) * no_of_days_per_year + day_of_hour (t) - 1];
                                                
                for (int tt = first_hour_of_day (day_of_hour (t)); tt <= t; ++tt) {
                    expr += inflow__in_sample[s][r][tt];
                    for (int a = 0; a < no_of_arcs; ++a)
                        expr += topology_matrix[r][a] *
                                (generate[(s * no_of_arcs + a) * no_of_hours_per_year + tt] -
                                 pump[(s * no_of_arcs + a) * no_of_hours_per_year + tt] +
                                 spill[(s * no_of_arcs + a) * no_of_hours_per_year + tt]);
                }
                
                model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
                
                if (t == last_hour_of_day (day_of_hour (t)))
                    model.addConstr (expr >= water_level[(s * no_of_reservoirs + r) * no_of_days_per_year + day_of_hour (t)]);
            }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: yearly cycle
    /*tic();
    cout << "Setting up yearly cycle constraints..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int r = 0; r < no_of_reservoirs; ++r)
            model.addConstr (water_level[(s * no_of_reservoirs + r) * no_of_days_per_year +
                                         (no_of_days_per_year - 1)] >= reservoir_level__beginning_of_year[r]);
    
    cout << "done: " << toc() << " seconds." << endl;*/

    // set up constraints: linear decision rules for water level targets
    tic();
    cout << "Setting up linear decision rules for water level targets..." << flush;
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int r = 0; r < no_of_reservoirs; ++r)
            for (int d = 0; d < no_of_days_per_year; ++d) {
                double coeff__avg_spot_curr_day = 0.0;
                for (int t = first_hour_of_day (d); t < first_hour_of_day (d + 1); ++t)
                    coeff__avg_spot_curr_day += spot_price__in_sample[s][t];
                
                expr = water_level__LDR__intercept[r * no_of_days_per_year + d] +
                       coeff__avg_spot_curr_day * water_level__LDR__slope_avg_spot_curr_day[r * no_of_days_per_year + d];
                
                for (int rr = 0; rr < no_of_reservoirs; ++rr) {
                    double coeff__cum_inflow_curr_day = 0.0;
                    for (int t = first_hour_of_day (d); t < first_hour_of_day (d + 1); ++t)
                        coeff__cum_inflow_curr_day += inflow__in_sample[s][rr][t];
                                            
                    double coeff__cum_inflow_past = 0.0;
                    for (int t = 0; t < first_hour_of_day (d); ++t)
                        coeff__cum_inflow_past += inflow__in_sample[s][rr][t];
                    
                    expr += coeff__cum_inflow_curr_day *
                            water_level__LDR__slope_cum_inflow_curr_day[(r * no_of_days_per_year + d) * no_of_reservoirs + rr] +
                            coeff__cum_inflow_past *
                            water_level__LDR__slope_cum_inflow_past[(r * no_of_days_per_year + d) * no_of_reservoirs + rr];
                }
                model.addConstr (water_level[(s * no_of_reservoirs + r) * no_of_days_per_year + d] == expr);
            }
    
    cout << "done: " << toc() << " seconds." << endl;
    
    // solve the problem
    model.optimize();

    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_primal_master_problem", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    for (int r = 0; r < no_of_reservoirs; ++r)
        cout << "WATER VALUE " << r << ": " << initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi) << endl;
    
    if (water_values != nullptr)
        for (int r = 0; r < no_of_reservoirs; ++r)
            (*water_values)[r] = initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi);

    ofstream fout (filename);
    for (int s = 0; s < no_of_insample_scenarios__model; ++s) {
      for (int t = 0; t < no_of_hours_per_year; ++t)
	fout << spot_price__in_sample[s][t] * bid__spot[s * no_of_hours_per_year + t].get (GRB_DoubleAttr_X) << " "
	     << prob_of_reserve_up_call * reserve_up_price__in_sample[s][t] * bid__reserve_up[s * no_of_hours_per_year + t].get (GRB_DoubleAttr_X) << " "
	     << prob_of_reserve_down_call * reserve_down_price__in_sample[s][t] * bid__reserve_down[s * no_of_hours_per_year + t].get (GRB_DoubleAttr_X) << endl;
      fout << endl;
    }

    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;

    delete []initial_water_level;
    delete []water_level;
    delete []water_level__LDR__intercept;
    delete []water_level__LDR__slope_cum_inflow_curr_day;
    delete []water_level__LDR__slope_cum_inflow_past;
    delete []water_level__LDR__slope_avg_spot_curr_day;
    
    // that's it!
    return obj_value;
}

double solve_primal_master_problem__deterministic__robust_bounds (bool include_reserve, vector<double> *water_values) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    tic();
    cout << "Setting up decision variables..." << flush;
    
    GRBVar *bid__spot         = model.addVars (no_of_hours_per_year);
    GRBVar *bid__reserve_up   = model.addVars (no_of_hours_per_year);
    GRBVar *bid__reserve_down = model.addVars (no_of_hours_per_year);

    GRBVar *generate = model.addVars (no_of_arcs * no_of_hours_per_year);
    GRBVar *pump     = model.addVars (no_of_arcs * no_of_hours_per_year);
    GRBVar *spill    = model.addVars (no_of_arcs * no_of_hours_per_year);

    GRBVar *initial_water_level = model.addVars (no_of_reservoirs);             // to obtain shadow prices
    GRBVar *water_level = model.addVars (no_of_reservoirs * no_of_days_per_year);
    
    cout << "done: " << toc() << " seconds." << endl;
    
    GRBLinExpr expr;
    
    // set up objective function
    tic();
    cout << "Setting up objective function..." << flush;
    
    expr = 0;
    for (int s = 0; s < no_of_insample_scenarios__data; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            expr += (spot_price__in_sample[s][t] * bid__spot[t] +
                     prob_of_reserve_up_call * reserve_up_price__in_sample[s][t] * bid__reserve_up[t] +
                     prob_of_reserve_down_call * reserve_down_price__in_sample[s][t] * bid__reserve_down[t])
                    / (double)(no_of_insample_scenarios__data);
    model.setObjective (expr, GRB_MAXIMIZE);
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up variable bounds
    tic();
    cout << "Setting up variable bounds..." << flush;

    for (int t = 0; t < no_of_hours_per_year; ++t)
        bid__spot[t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            bid__reserve_up[t].set (GRB_DoubleAttr_UB, 0.0);
            bid__reserve_down[t].set (GRB_DoubleAttr_UB, 0.0);
        }
    }

    for (int t = 0; t < no_of_hours_per_year; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            generate[a * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            pump[a * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int d = 0; d < no_of_days_per_year; ++d) {
            water_level[r * no_of_days_per_year + d].set (GRB_DoubleAttr_LB, lb_reservoir_level[r]);
            water_level[r * no_of_days_per_year + d].set (GRB_DoubleAttr_UB, ub_reservoir_level[r]);
        }
        
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: initial water level (for shadow prices)
    tic();
    cout << "Setting up initial water level constraints..." << flush;

    vector<GRBConstr> initial_water_level_constraints;
    for (int r = 0; r < no_of_reservoirs; ++r)
        initial_water_level_constraints.push_back (model.addConstr (initial_water_level[r] == initial_reservoir_level[r]));
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: robust market fulfillment
    tic();
    cout << "Setting up robust market fulfillment constraints..." << flush;

    double maximum_pumping_quantity = 0.0, maximum_generation_quantity = 0.0;
    for (int a = 0; a < no_of_arcs; ++a) {
        maximum_pumping_quantity += inverse_pump_efficiency[a] * ub_pumping[a];
        maximum_generation_quantity += generator_efficiency[a] * ub_generation[a];
    }
    
    for (int t = 0; t < no_of_hours_per_year; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * generate[a * no_of_hours_per_year + t] -
                    inverse_pump_efficiency[a] * pump[a * no_of_hours_per_year + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t]
                                      - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
        
        model.addConstr (bid__spot[t] - bid__reserve_down[t] + maximum_pumping_quantity >= 0.0);
        model.addConstr (bid__spot[t] + bid__reserve_up[t] <= maximum_generation_quantity);

        //model.addConstr (bid__reserve_up[t] + bid__reserve_down[t] <= maximum_reserve_investment * bid__spot[t]);
    }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: hourly water level dynamics
    tic();
    cout << "Setting up hourly water level dynamics constraints..." << flush;

    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            if (day_of_hour (t) == 0)
                expr = initial_water_level[r];
            else
                expr = water_level[r * no_of_days_per_year + day_of_hour (t) - 1];
                                            
            for (int tt = first_hour_of_day (day_of_hour (t)); tt <= t; ++tt) {
                for (int s = 0; s < no_of_insample_scenarios__data; ++s)
                    expr += inflow__in_sample[s][r][tt] / (double)(no_of_insample_scenarios__data);
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] *
                            (generate[a * no_of_hours_per_year + tt] -
                             pump[a * no_of_hours_per_year + tt] +
                             spill[a * no_of_hours_per_year + tt]);
            }
            
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
            
            if (t == last_hour_of_day (day_of_hour (t)))
                model.addConstr (expr >= water_level[r * no_of_days_per_year + day_of_hour (t)]);
        }
    
    cout << "done: " << toc() << " seconds." << endl;
    
    // solve the problem
    model.optimize();

    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_primal_master_problem__deterministic", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    for (int r = 0; r < no_of_reservoirs; ++r)
        cout << "WATER VALUE " << r << ": " << initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi) << endl;
    
    if (water_values != nullptr)
        for (int r = 0; r < no_of_reservoirs; ++r)
            (*water_values)[r] = initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi);

    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;

    delete []initial_water_level;
    delete []water_level;
    
    // that's it!
    return obj_value;
}

double solve_primal_master_problem__deterministic__percentage_of_spot (bool include_reserve, vector<double> *water_values) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    tic();
    cout << "Setting up decision variables..." << flush;
    
    GRBVar *bid__spot         = model.addVars (no_of_hours_per_year);
    GRBVar *bid__reserve_up   = model.addVars (no_of_hours_per_year);
    GRBVar *bid__reserve_down = model.addVars (no_of_hours_per_year);

    GRBVar *generate = model.addVars (no_of_arcs * no_of_hours_per_year);
    GRBVar *pump     = model.addVars (no_of_arcs * no_of_hours_per_year);
    GRBVar *spill    = model.addVars (no_of_arcs * no_of_hours_per_year);

    GRBVar *initial_water_level = model.addVars (no_of_reservoirs);             // to obtain shadow prices
    GRBVar *water_level = model.addVars (no_of_reservoirs * no_of_days_per_year);
    
    cout << "done: " << toc() << " seconds." << endl;
    
    GRBLinExpr expr;
    
    // set up objective function
    tic();
    cout << "Setting up objective function..." << flush;
    
    expr = 0;
    for (int s = 0; s < no_of_insample_scenarios__data; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            expr += (spot_price__in_sample[s][t] * bid__spot[t] +
                     prob_of_reserve_up_call * reserve_up_price__in_sample[s][t] * bid__reserve_up[t] +
                     prob_of_reserve_down_call * reserve_down_price__in_sample[s][t] * bid__reserve_down[t])
                    / (double)(no_of_insample_scenarios__data);
    model.setObjective (expr, GRB_MAXIMIZE);
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up variable bounds
    tic();
    cout << "Setting up variable bounds..." << flush;

    for (int t = 0; t < no_of_hours_per_year; ++t)
        bid__spot[t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            bid__reserve_up[t].set (GRB_DoubleAttr_UB, 0.0);
            bid__reserve_down[t].set (GRB_DoubleAttr_UB, 0.0);
        }
    }

    for (int t = 0; t < no_of_hours_per_year; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            generate[a * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            pump[a * no_of_hours_per_year + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int d = 0; d < no_of_days_per_year; ++d) {
            water_level[r * no_of_days_per_year + d].set (GRB_DoubleAttr_LB, lb_reservoir_level[r]);
            water_level[r * no_of_days_per_year + d].set (GRB_DoubleAttr_UB, ub_reservoir_level[r]);
        }
        
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: initial water level (for shadow prices)
    tic();
    cout << "Setting up initial water level constraints..." << flush;

    vector<GRBConstr> initial_water_level_constraints;
    for (int r = 0; r < no_of_reservoirs; ++r)
        initial_water_level_constraints.push_back (model.addConstr (initial_water_level[r] == initial_reservoir_level[r]));
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: robust market fulfillment
    tic();
    cout << "Setting up robust market fulfillment constraints..." << flush;

    for (int t = 0; t < no_of_hours_per_year; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * generate[a * no_of_hours_per_year + t] -
                    inverse_pump_efficiency[a] * pump[a * no_of_hours_per_year + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t]
                                      - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
        
        model.addConstr (bid__reserve_up[t] + bid__reserve_down[t] <= deterministic__percentage_of_spot__factor * bid__spot[t]);
    }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: hourly water level dynamics
    tic();
    cout << "Setting up hourly water level dynamics constraints..." << flush;

    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            if (day_of_hour (t) == 0)
                expr = initial_water_level[r];
            else
                expr = water_level[r * no_of_days_per_year + day_of_hour (t) - 1];
                                            
            for (int tt = first_hour_of_day (day_of_hour (t)); tt <= t; ++tt) {
                for (int s = 0; s < no_of_insample_scenarios__data; ++s)
                    expr += inflow__in_sample[s][r][tt] / (double)(no_of_insample_scenarios__data);
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] *
                            (generate[a * no_of_hours_per_year + tt] -
                             pump[a * no_of_hours_per_year + tt] +
                             spill[a * no_of_hours_per_year + tt]);
            }
            
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
            
            if (t == last_hour_of_day (day_of_hour (t)))
                model.addConstr (expr >= water_level[r * no_of_days_per_year + day_of_hour (t)]);
        }
    
    cout << "done: " << toc() << " seconds." << endl;
    
    // solve the problem
    model.optimize();

    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_primal_master_problem__deterministic", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    for (int r = 0; r < no_of_reservoirs; ++r)
        cout << "WATER VALUE " << r << ": " << initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi) << endl;
    
    if (water_values != nullptr)
        for (int r = 0; r < no_of_reservoirs; ++r)
            (*water_values)[r] = initial_water_level_constraints[r].get (GRB_DoubleAttr_Pi);

    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;

    delete []initial_water_level;
    delete []water_level;
    
    // that's it!
    return obj_value;
}
