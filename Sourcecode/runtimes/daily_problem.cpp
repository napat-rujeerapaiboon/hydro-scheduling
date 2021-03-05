//
//  daily_problem.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "daily_problem.hpp"

#include <sstream>

#include "gurobi_c++.h"

#include "auxiliary.hpp"
#include "data_and_parameters.hpp"

using namespace std;

double solve_first_hour_of_daily_problem__stochastic (const vector<double> &water_values, bool include_reserve, vector<double> *result_s, vector<double> *result_u, vector<double> *result_v) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *bid__spot         = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_up   = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_down = model.addVars (no_of_hours_per_day);

    GRBVar *generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *spill       = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *water_flow  = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBVar *exp_generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *exp_pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *exp_spill       = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *exp_water_flow  = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBLinExpr expr, expr2;
    
    // set up objective function
    expr = 0;
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr += spot_price__in_sample[0][t] * bid__spot[t] +
                prob_of_reserve_up_call * reserve_up_price__in_sample[0][t] * bid__reserve_up[t] +
                prob_of_reserve_down_call * reserve_down_price__in_sample[0][t] * bid__reserve_down[t];
    
        for (int r = 0; r < no_of_reservoirs; ++r) {
            expr2 = inflow__in_sample[0][r][t];
            for (int a = 0; a < no_of_arcs; ++a)
                expr2 += topology_matrix[r][a] * exp_water_flow[a * no_of_hours_per_day + t];
            expr += water_values[r] * expr2;
        }
    }

    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        bid__spot[t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            bid__reserve_up[t].set (GRB_DoubleAttr_UB, 0.0);
            bid__reserve_down[t].set (GRB_DoubleAttr_UB, 0.0);
        }
    }
    
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);

            exp_generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            exp_pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            water_flow[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            exp_water_flow[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
        }

    // set up constraints: f-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            model.addConstr (water_flow[a * no_of_hours_per_day + t] == generate[a * no_of_hours_per_day + t] -
                                                                        pump[a * no_of_hours_per_day + t] +
                                                                        spill[a * no_of_hours_per_day + t]);
            model.addConstr (exp_water_flow[a * no_of_hours_per_day + t] == exp_generate[a * no_of_hours_per_day + t] -
                                                                            exp_pump[a * no_of_hours_per_day + t] +
                                                                            exp_spill[a * no_of_hours_per_day + t]);
        }
    
    // set up constraints: suv-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + bid__reserve_up[t] == expr);
        
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += inverse_pump_efficiency[a] * ub_pumping[a];
        model.addConstr (bid__spot[t] - bid__reserve_down[t] >= -expr);
        
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * exp_generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * exp_pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t] - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
    }
    
    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] * water_flow[a * no_of_hours_per_day + tt];
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);

            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] * exp_water_flow[a * no_of_hours_per_day + tt];
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();
    
    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_first_hour_of_daily_problem__stochastic", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    if (result_s != nullptr) {
        (*result_s).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_s).push_back (bid__spot[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_u != nullptr) {
        (*result_u).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_u).push_back (bid__reserve_up[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_v != nullptr) {
        (*result_v).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_v).push_back (bid__reserve_down[t].get (GRB_DoubleAttr_X));
    }

    runtime_slave = model.get (GRB_DoubleAttr_Runtime);
    
    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;
    delete []water_flow;

    delete []exp_generate;
    delete []exp_pump;
    delete []exp_spill;
    delete []exp_water_flow;

    // that's it!
    return obj_value;
}

double solve_first_hour_of_daily_problem__deterministic__robust_bounds (const vector<double> &water_values, bool include_reserve, vector<double> *result_s, vector<double> *result_u, vector<double> *result_v) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *bid__spot         = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_up   = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_down = model.addVars (no_of_hours_per_day);

    GRBVar *generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *spill       = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *water_flow  = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBLinExpr expr, expr2;
    
    // set up objective function
    expr = 0;
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr += spot_price__in_sample[0][t] * bid__spot[t] +
                prob_of_reserve_up_call * reserve_up_price__in_sample[0][t] * bid__reserve_up[t] +
                prob_of_reserve_down_call * reserve_down_price__in_sample[0][t] * bid__reserve_down[t];
    
        for (int r = 0; r < no_of_reservoirs; ++r) {
            expr2 = inflow__in_sample[0][r][t];
            for (int a = 0; a < no_of_arcs; ++a)
                expr2 += topology_matrix[r][a] * water_flow[a * no_of_hours_per_day + t];
            expr += water_values[r] * expr2;
        }
    }

    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        bid__spot[t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            bid__reserve_up[t].set (GRB_DoubleAttr_UB, 0.0);
            bid__reserve_down[t].set (GRB_DoubleAttr_UB, 0.0);
        }
    }
    
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a)
            water_flow[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    // set up constraints: f-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a)
            model.addConstr (water_flow[a * no_of_hours_per_day + t] == generate[a * no_of_hours_per_day + t] -
                                                                        pump[a * no_of_hours_per_day + t] +
                                                                        spill[a * no_of_hours_per_day + t]);
    
    // set up constraints: suv-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t] - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
        
        double maximum_pumping_quantity = 0.0, maximum_generation_quantity = 0.0;
        for (int a = 0; a < no_of_arcs; ++a) {
            maximum_pumping_quantity += inverse_pump_efficiency[a] * ub_pumping[a];
            maximum_generation_quantity += generator_efficiency[a] * ub_generation[a];
        }

        model.addConstr (bid__spot[t] - bid__reserve_down[t] + maximum_pumping_quantity >= 0.0);
        model.addConstr (bid__spot[t] + bid__reserve_up[t] <= maximum_generation_quantity);
        
        //model.addConstr (bid__reserve_up[t] + bid__reserve_down[t] <= maximum_reserve_investment * bid__spot[t]);
    }
    
    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] * water_flow[a * no_of_hours_per_day + tt];
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();
    
    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_first_hour_of_daily_problem__deterministic", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    if (result_s != nullptr) {
        (*result_s).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_s).push_back (bid__spot[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_u != nullptr) {
        (*result_u).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_u).push_back (bid__reserve_up[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_v != nullptr) {
        (*result_v).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_v).push_back (bid__reserve_down[t].get (GRB_DoubleAttr_X));
    }
    
    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;
    delete []water_flow;

    // that's it!
    return obj_value;
}

double solve_first_hour_of_daily_problem__deterministic__percentage_of_spot (const vector<double> &water_values, bool include_reserve, vector<double> *result_s, vector<double> *result_u, vector<double> *result_v) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *bid__spot         = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_up   = model.addVars (no_of_hours_per_day);
    GRBVar *bid__reserve_down = model.addVars (no_of_hours_per_day);

    GRBVar *generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *spill       = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *water_flow  = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBLinExpr expr, expr2;
    
    // set up objective function
    expr = 0;
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr += spot_price__in_sample[0][t] * bid__spot[t] +
                prob_of_reserve_up_call * reserve_up_price__in_sample[0][t] * bid__reserve_up[t] +
                prob_of_reserve_down_call * reserve_down_price__in_sample[0][t] * bid__reserve_down[t];
    
        for (int r = 0; r < no_of_reservoirs; ++r) {
            expr2 = inflow__in_sample[0][r][t];
            for (int a = 0; a < no_of_arcs; ++a)
                expr2 += topology_matrix[r][a] * water_flow[a * no_of_hours_per_day + t];
            expr += water_values[r] * expr2;
        }
    }

    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        bid__spot[t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    if (!include_reserve) {
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            bid__reserve_up[t].set (GRB_DoubleAttr_UB, 0.0);
            bid__reserve_down[t].set (GRB_DoubleAttr_UB, 0.0);
        }
    }
    
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a)
            water_flow[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);

    // set up constraints: f-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a)
            model.addConstr (water_flow[a * no_of_hours_per_day + t] == generate[a * no_of_hours_per_day + t] -
                                                                        pump[a * no_of_hours_per_day + t] +
                                                                        spill[a * no_of_hours_per_day + t]);
    
    // set up constraints: suv-constraints
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t] - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
        
        double maximum_pumping_quantity = 0.0, maximum_generation_quantity = 0.0;
        for (int a = 0; a < no_of_arcs; ++a) {
            maximum_pumping_quantity += inverse_pump_efficiency[a] * ub_pumping[a];
            maximum_generation_quantity += generator_efficiency[a] * ub_generation[a];
        }

        model.addConstr (bid__reserve_up[t] + bid__reserve_down[t] <= deterministic__percentage_of_spot__factor * bid__spot[t]);
    }
    
    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = 0; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a)
                    expr += topology_matrix[r][a] * water_flow[a * no_of_hours_per_day + tt];
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();
    
    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_first_hour_of_daily_problem__deterministic__percentage_of_spot", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    if (result_s != nullptr) {
        (*result_s).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_s).push_back (bid__spot[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_u != nullptr) {
        (*result_u).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_u).push_back (bid__reserve_up[t].get (GRB_DoubleAttr_X));
    }
    
    if (result_v != nullptr) {
        (*result_v).clear();
        for (int t = 0; t < no_of_hours_per_day; ++t)
            (*result_v).push_back (bid__reserve_down[t].get (GRB_DoubleAttr_X));
    }
    
    // de-allocate variables
    delete []bid__spot;
    delete []bid__reserve_up;
    delete []bid__reserve_down;

    delete []generate;
    delete []pump;
    delete []spill;
    delete []water_flow;

    // that's it!
    return obj_value;
}

double solve_each_hour_of_daily_problem__stochastic (int theta, const vector<double> &water_values, const vector<double> &bid__spot, const vector<double> &bid__reserve_up, const vector<double> &bid__reserve_down, vector<double> &res_generate, vector<double> &res_pump, vector<double> &res_spill) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *var_generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_spill       = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBVar *var_exp_generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_exp_pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_exp_spill       = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBLinExpr expr, expr2;

    // set up objective function
    expr = 0;
    
    for (int r = 0; r < no_of_reservoirs; ++r) {
        for (int a = 0; a < no_of_arcs; ++a) {
            expr += water_values[r] * topology_matrix[r][a] * (var_generate[a * no_of_hours_per_day + theta] -
                                                               var_pump[a * no_of_hours_per_day + theta] +
                                                               var_spill[a * no_of_hours_per_day + theta]);
            
            for (int t = theta + 1; t < no_of_hours_per_day; ++t)
                expr += water_values[r] * topology_matrix[r][a] * (var_exp_generate[a * no_of_hours_per_day + theta] -
                                                                   var_exp_pump[a * no_of_hours_per_day + theta] +
                                                                   var_exp_spill[a * no_of_hours_per_day + theta]);
        }
    }

    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            var_generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            var_pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);

            var_exp_generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            var_exp_pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }

    // set up constraints: suv-constraint for theta
    expr = bid__spot[theta];
    if (reserve_up_call__out_of_sample[theta]) expr += bid__reserve_up[theta];
    if (reserve_down_call__out_of_sample[theta]) expr -= bid__reserve_down[theta];
    
    expr2 = 0;
    for (int a = 0; a < no_of_arcs; ++a)
        expr2 += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + theta] -
                 inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + theta];
    model.addConstr (expr == expr2);
    
    // set up constraints: suv-constraint for t > theta
    for (int t = theta + 1; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + bid__reserve_up[t] == expr);
        
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * var_exp_generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * var_exp_pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t]
                                      - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
    }

    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = theta; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a) {
                    if (tt >= theta)
                        expr += topology_matrix[r][a] * (var_generate[a * no_of_hours_per_day + tt] -
                                                         var_pump[a * no_of_hours_per_day + tt] +
                                                         var_spill[a * no_of_hours_per_day + tt]);
                    else
                        expr += topology_matrix[r][a] * (res_generate[a * no_of_hours_per_day + tt] -
                                                         res_pump[a * no_of_hours_per_day + tt] +
                                                         res_spill[a * no_of_hours_per_day + tt]);
                }
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);

            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a) {
                    if (tt >= theta)
                        expr += topology_matrix[r][a] * (var_exp_generate[a * no_of_hours_per_day + tt] -
                                                         var_exp_pump[a * no_of_hours_per_day + tt] +
                                                         var_exp_spill[a * no_of_hours_per_day + tt]);
                    else
                        expr += topology_matrix[r][a] * (res_generate[a * no_of_hours_per_day + tt] -
                                                         res_pump[a * no_of_hours_per_day + tt] +
                                                         res_spill[a * no_of_hours_per_day + tt]);
                }
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();
    
    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_each_hour_of_daily_problem__stochastic", "Optimization problem not solved to optimality.");
    }

    double obj_value = model.get(GRB_DoubleAttr_ObjVal);

    runtime_slave += model.get (GRB_DoubleAttr_Runtime);
    
    // record the generation, pump and spill decisions
    for (int a = 0; a < no_of_arcs; ++a) {
        res_generate[a * no_of_hours_per_day + theta] = var_generate[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
        res_pump[a * no_of_hours_per_day + theta] = var_pump[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
        res_spill[a * no_of_hours_per_day + theta] = var_spill[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
    }
    
    // de-allocate variables
    delete []var_generate;
    delete []var_pump;
    delete []var_spill;

    delete []var_exp_generate;
    delete []var_exp_pump;
    delete []var_exp_spill;

    // that's it!
    return obj_value;
}

// problem = 1: original obj fct; attempt to honour all reserve market commitments; problem may be infeasible
// problem = 2: keep all reserve market commitments flexible, maximize here-and-now commitment; should always be feasible
// problem = 3: original obj fct; fix here-and-now commitment to (val_u, val_v), keep future commitments flexible; should always be feasible
pair<bool, double> solve_each_hour_of_daily_problem__deterministic__robust_bounds__subproblem (int theta, const vector<double> &water_values,
                            const vector<double> &bid__spot, const vector<double> &bid__reserve_up, const vector<double> &bid__reserve_down,
                            vector<double> &res_generate, vector<double> &res_pump, vector<double> &res_spill,
                            int problem, double *val_u, double *val_v) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *var_generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_spill       = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBVar *var_u           = model.addVars (no_of_hours_per_day);
    GRBVar *var_v           = model.addVars (no_of_hours_per_day);

    GRBLinExpr expr, expr2;

    // set up objective function
    expr = 0;
    
    if (problem == 2) {
        expr = var_u[theta] + var_v[theta];
    } else {
        for (int r = 0; r < no_of_reservoirs; ++r)
            for (int a = 0; a < no_of_arcs; ++a)
                for (int t = theta; t < no_of_hours_per_day; ++t)
                    expr += water_values[r] * topology_matrix[r][a] *
                            (var_generate[a * no_of_hours_per_day + theta] -
                             var_pump[a * no_of_hours_per_day + theta] +
                             var_spill[a * no_of_hours_per_day + theta]);
    }

    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            var_generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            var_pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }
    
    // problem = 1: original obj fct; attempt to honour all reserve market commitments; problem may be infeasible
    // problem = 2: keep all reserve market commitments flexible, maximize here-and-now commitment; should always be feasible
    // problem = 3: original obj fct; fix here-and-now commitment to (val_u, val_v), keep future commitments flexible; should always be feasible
    if (problem == 1) {
        for (int t = theta; t < no_of_hours_per_day; ++t) {
            var_u[t].set (GRB_DoubleAttr_LB, bid__reserve_up[t]);
            var_v[t].set (GRB_DoubleAttr_LB, bid__reserve_down[t]);

            var_u[t].set (GRB_DoubleAttr_UB, bid__reserve_up[t]);
            var_v[t].set (GRB_DoubleAttr_UB, bid__reserve_down[t]);
        }
    } else if (problem == 2) {
        for (int t = theta; t < no_of_hours_per_day; ++t) {
            var_u[t].set (GRB_DoubleAttr_UB, bid__reserve_up[t]);
            var_v[t].set (GRB_DoubleAttr_UB, bid__reserve_down[t]);
        }
    } else {
        var_u[theta].set (GRB_DoubleAttr_LB, *val_u);
        var_v[theta].set (GRB_DoubleAttr_LB, *val_v);

        for (int t = theta; t < no_of_hours_per_day; ++t) {
            var_u[t].set (GRB_DoubleAttr_UB, bid__reserve_up[t]);
            var_v[t].set (GRB_DoubleAttr_UB, bid__reserve_down[t]);
        }
    }
    
    // set up constraints: suv-constraint for theta
    expr = bid__spot[theta];
    if (reserve_up_call__out_of_sample[theta]) expr += var_u[theta];
    if (reserve_down_call__out_of_sample[theta]) expr -= var_v[theta];
    
    expr2 = 0;
    for (int a = 0; a < no_of_arcs; ++a)
        expr2 += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + theta] -
                 inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + theta];
    model.addConstr (expr == expr2);
    
    // set up constraints: suv-constraint for t > theta
    for (int t = theta + 1; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * var_u[t]
                                      - prob_of_reserve_down_call * var_v[t] == expr);
    }

    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = theta; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a) {
                    if (tt >= theta)
                        expr += topology_matrix[r][a] * (var_generate[a * no_of_hours_per_day + tt] -
                                                         var_pump[a * no_of_hours_per_day + tt] +
                                                         var_spill[a * no_of_hours_per_day + tt]);
                    else
                        expr += topology_matrix[r][a] * (res_generate[a * no_of_hours_per_day + tt] -
                                                         res_pump[a * no_of_hours_per_day + tt] +
                                                         res_spill[a * no_of_hours_per_day + tt]);
                }
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();
    stringstream fname;
    fname << "problem_" << problem << ".lp";
    model.write (fname.str().c_str());

    if ((model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) && (model.get (GRB_IntAttr_Status) != GRB_INFEASIBLE)) {
        model.write ("error.lp");
        error ("solve_each_hour_of_daily_problem__deterministic__subproblem", "Optimization problem neither proven infeasible nor solved to optimality.");
    }

    pair<bool, double> result;
    result.first = (model.get (GRB_IntAttr_Status) == GRB_OPTIMAL);
    
    if (result.first) {
        result.second = model.get(GRB_DoubleAttr_ObjVal);
    
        // record the generation, pump and spill decisions
        for (int a = 0; a < no_of_arcs; ++a) {
            res_generate[a * no_of_hours_per_day + theta] = var_generate[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
            res_pump[a * no_of_hours_per_day + theta] = var_pump[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
            res_spill[a * no_of_hours_per_day + theta] = var_spill[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
        }
    } else result.second = -10000.0;
    
    if (problem == 2) {
        *val_u = var_u[theta].get (GRB_DoubleAttr_X);
        *val_v = var_v[theta].get (GRB_DoubleAttr_X);
    }
    
    // de-allocate variables
    delete []var_generate;
    delete []var_pump;
    delete []var_spill;

    // that's it!
    return result;
}

double solve_each_hour_of_daily_problem__deterministic__robust_bounds (int theta, const vector<double> &water_values, const vector<double> &bid__spot, const vector<double> &bid__reserve_up, const vector<double> &bid__reserve_down, vector<double> &res_generate, vector<double> &res_pump, vector<double> &res_spill,
    double *missed_commitment_u, double *missed_commitment_v) {
    // problem = 1: original obj fct; attempt to honour all reserve market commitments; problem may be infeasible
    cout << "[[ SOLVING PROBLEM 1 FOR HOUR " << theta << " ]]" << endl;
    pair<bool, double> result;
    double val_u = -1.0, val_v = -1.0;
    result = solve_each_hour_of_daily_problem__deterministic__robust_bounds__subproblem
                (theta, water_values, bid__spot, bid__reserve_up, bid__reserve_down,
                 res_generate, res_pump, res_spill,
                 1, &val_u, &val_v);
    if (result.first) {
        *missed_commitment_u = 0.0;
        *missed_commitment_v = 0.0;
        return result.second;
    }
    
    // problem = 2: keep all reserve market commitments flexible, maximize here-and-now commitment; should always be feasible
    cout << "[[ SOLVING PROBLEM 2 FOR HOUR " << theta << " ]]" << endl;
    result = solve_each_hour_of_daily_problem__deterministic__robust_bounds__subproblem
                (theta, water_values, bid__spot, bid__reserve_up, bid__reserve_down,
                 res_generate, res_pump, res_spill,
                 2, &val_u, &val_v);
    if (!result.first) error ("solve_each_hour_of_daily_problem__deterministic", "Expected subproblem 2 to be solvable.");
    
    *missed_commitment_u = bid__reserve_up[theta] - val_u;
    *missed_commitment_v = bid__reserve_down[theta] - val_v;

    cout << "*** " << *missed_commitment_u << " " << *missed_commitment_v << endl;
    val_u -= 0.00001; if (val_u < 0.0) val_u = 0.0;
    val_v -= 0.00001; if (val_v < 0.0) val_v = 0.0;
    cout << "*** " << val_u << " " << val_v << endl;

    // Problem 3: optimize original objective while honouring the reserve market commitments from Problem 2
    cout << "[[ SOLVING PROBLEM 3 FOR HOUR " << theta << " ]]" << endl;
    result = solve_each_hour_of_daily_problem__deterministic__robust_bounds__subproblem
                (theta, water_values, bid__spot, bid__reserve_up, bid__reserve_down,
                 res_generate, res_pump, res_spill,
                 3, &val_u, &val_v);
    if (!result.first) error ("solve_each_hour_of_daily_problem__deterministic", "Expected subproblem 3 to be solvable.");
    
    return result.second;
}

double solve_each_hour_of_daily_problem__deterministic__percentage_of_spot (int theta, const vector<double> &water_values, const vector<double> &bid__spot, const vector<double> &bid__reserve_up, const vector<double> &bid__reserve_down, vector<double> &res_generate, vector<double> &res_pump, vector<double> &res_spill) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    //model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    //model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    GRBVar *var_generate    = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_pump        = model.addVars (no_of_arcs * no_of_hours_per_day);
    GRBVar *var_spill       = model.addVars (no_of_arcs * no_of_hours_per_day);

    GRBLinExpr expr, expr2;

    // set up objective function
    expr = 0;
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int a = 0; a < no_of_arcs; ++a)
            for (int t = theta; t < no_of_hours_per_day; ++t)
                expr += water_values[r] * topology_matrix[r][a] *
                        (var_generate[a * no_of_hours_per_day + theta] -
                         var_pump[a * no_of_hours_per_day + theta] +
                         var_spill[a * no_of_hours_per_day + theta]);
    model.setObjective (expr, GRB_MAXIMIZE);

    // set up variable bounds
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int a = 0; a < no_of_arcs; ++a) {
            var_generate[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_generation[a]);
            var_pump[a * no_of_hours_per_day + t].set (GRB_DoubleAttr_UB, ub_pumping[a]);
        }
    
    // set up constraints: suv-constraint for theta
    expr = bid__spot[theta];
    if (reserve_up_call__out_of_sample[theta]) expr += bid__reserve_up[theta];
    if (reserve_down_call__out_of_sample[theta]) expr -= bid__reserve_down[theta];
    
    expr2 = 0;
    for (int a = 0; a < no_of_arcs; ++a)
        expr2 += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + theta] -
                 inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + theta];
    model.addConstr (expr == expr2);
    
    // set up constraints: suv-constraint for t > theta
    for (int t = theta + 1; t < no_of_hours_per_day; ++t) {
        expr = 0;
        for (int a = 0; a < no_of_arcs; ++a)
            expr += generator_efficiency[a] * var_generate[a * no_of_hours_per_day + t] -
                    inverse_pump_efficiency[a] * var_pump[a * no_of_hours_per_day + t];
        model.addConstr (bid__spot[t] + prob_of_reserve_up_call * bid__reserve_up[t]
                                      - prob_of_reserve_down_call * bid__reserve_down[t] == expr);
    }

    // set up constraints: hourly water level dynamics
    for (int r = 0; r < no_of_reservoirs; ++r)
        for (int t = theta; t < no_of_hours_per_day; ++t) {
            expr = initial_reservoir_level[r];
            for (int tt = 0; tt <= t; ++tt) {
                expr += inflow__in_sample[0][r][tt];
                for (int a = 0; a < no_of_arcs; ++a) {
                    if (tt >= theta)
                        expr += topology_matrix[r][a] * (var_generate[a * no_of_hours_per_day + tt] -
                                                         var_pump[a * no_of_hours_per_day + tt] +
                                                         var_spill[a * no_of_hours_per_day + tt]);
                    else
                        expr += topology_matrix[r][a] * (res_generate[a * no_of_hours_per_day + tt] -
                                                         res_pump[a * no_of_hours_per_day + tt] +
                                                         res_spill[a * no_of_hours_per_day + tt]);
                }
            }
            model.addRange (expr, lb_reservoir_level[r], ub_reservoir_level[r]);
        }

    // solve the problem
    model.optimize();

    if ((model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) && (model.get (GRB_IntAttr_Status) != GRB_INFEASIBLE)) {
        model.write ("error.lp");
        error ("solve_each_hour_of_daily_problem__deterministic__robust_bounds", "Optimization problem neither proven infeasible nor solved to optimality.");
    }

    if (model.get (GRB_IntAttr_Status) == GRB_OPTIMAL) {
        // record the generation, pump and spill decisions
        for (int a = 0; a < no_of_arcs; ++a) {
            res_generate[a * no_of_hours_per_day + theta] = var_generate[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
            res_pump[a * no_of_hours_per_day + theta] = var_pump[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
            res_spill[a * no_of_hours_per_day + theta] = var_spill[a * no_of_hours_per_day + theta].get (GRB_DoubleAttr_X);
        }

        // de-allocate variables
        delete []var_generate;
        delete []var_pump;
        delete []var_spill;

        return model.get(GRB_DoubleAttr_ObjVal);
    } else return -10000.0;
}
