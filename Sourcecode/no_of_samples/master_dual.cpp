//
//  master_dual.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "master_dual.hpp"

#include "gurobi_c++.h"

#include "auxiliary.hpp"
#include "data_and_parameters.hpp"

using namespace std;

double solve_dual_master_problem__stochastic (bool include_reserve, vector<double> *water_values) {
    // set up model
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel (*env);

    model.set (GRB_IntParam_Method, 2);
    model.set (GRB_IntParam_Threads, 4);
    model.set (GRB_IntParam_Crossover, 0);

    // set up decision variables
    tic();
    cout << "Setting up decision variables..." << flush;
    
    GRBVar *alpha     = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);
    GRBVar *alpha_hat = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);
    
    GRBVar *beta     = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);
    GRBVar *beta_hat = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);
    
    GRBVar *gamma     = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);
    GRBVar *gamma_hat = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_arcs);

    GRBVar *kappa       = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);
    GRBVar *kappa_hat   = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);
    GRBVar *kappa_prime = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year);

    GRBVar *nu_plus      = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_reservoirs);
    GRBVar *nu_hat_plus  = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_reservoirs);
    GRBVar *nu_minus     = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_reservoirs);
    GRBVar *nu_hat_minus = model.addVars (no_of_insample_scenarios__model * no_of_hours_per_year * no_of_reservoirs);

    GRBVar *lambda_0 = model.addVars (no_of_reservoirs);
    GRBVar *lambda_D = model.addVars (no_of_reservoirs);
    GRBVar *lambda   = model.addVars (no_of_insample_scenarios__model * no_of_days_per_year * no_of_reservoirs);
    
    GRBVar *Phi     = model.addVars (no_of_days_per_year * no_of_reservoirs * no_of_reservoirs);
    GRBVar *Psi     = model.addVars (no_of_days_per_year * no_of_reservoirs * no_of_reservoirs);
    GRBVar *epsilon = model.addVars (no_of_days_per_year * no_of_reservoirs);
    GRBVar *theta   = model.addVars (no_of_days_per_year * no_of_reservoirs);
    
    GRBVar *mu_plus  = model.addVars (no_of_insample_scenarios__model * no_of_days_per_year * no_of_reservoirs);
    GRBVar *mu_minus = model.addVars (no_of_insample_scenarios__model * no_of_days_per_year * no_of_reservoirs);

    cout << "done: " << toc() << " seconds." << endl;

    GRBLinExpr expr;

    // set up objective function
    tic();
    cout << "Setting up objective function..." << flush;
    
    expr = 0;
    for (int s = 0; s < no_of_insample_scenarios__model; ++s) {
        for (int d = 0; d < no_of_days_per_year; ++d)
            for (int r = 0; r < no_of_reservoirs; ++r)
                expr += (d != 0 ? lb_reservoir_level[r] : 0.0) * mu_minus[(s * no_of_days_per_year + d) * no_of_reservoirs + r] +
                        (d != 0 ? ub_reservoir_level[r] : 0.0) * mu_plus[(s * no_of_days_per_year + d) * no_of_reservoirs + r];
        
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            for (int r = 0; r < no_of_reservoirs; ++r) {
                if (day_of_hour (t) < no_of_days_per_year - 1)
                    expr += inflow__in_sample[s][r][t] * lambda[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r];
                else
                    expr += inflow__in_sample[s][r][t] * lambda_D[r];
            }
            
            for (int a = 0; a < no_of_arcs; ++a)
                expr += ub_generation[a] * (alpha[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                            alpha_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a]) +
                        ub_pumping[a] * (beta[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                         beta_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                         inverse_pump_efficiency[a] * kappa_prime[s * no_of_hours_per_year + t]);
            
            for (int r = 0; r < no_of_reservoirs; ++r) {
                double cum_inflows = 0.0;
                for (int tt = first_hour_of_day (day_of_hour (t)); tt <= t; ++tt)
                    cum_inflows += inflow__in_sample[s][r][tt];
                
                expr += (nu_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r] +
                         nu_hat_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r]) *
                        (lb_reservoir_level[r] - cum_inflows) +
                        (nu_plus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r] +
                         nu_hat_plus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r]) *
                        (ub_reservoir_level[r] - cum_inflows);
            }
        }
    }
    expr /= (double)(no_of_insample_scenarios__model);
    model.setObjective (expr, GRB_MINIMIZE);
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up variable bounds
    tic();
    cout << "Setting up variable bounds..." << flush;
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int a = 0; a < no_of_arcs; ++a) {
                gamma[(s * no_of_hours_per_year + t) * no_of_arcs + a].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                gamma[(s * no_of_hours_per_year + t) * no_of_arcs + a].set (GRB_DoubleAttr_UB, 0.0);

                gamma_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                gamma_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a].set (GRB_DoubleAttr_UB, 0.0);
            }

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            kappa[s * no_of_hours_per_year + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            kappa_hat[s * no_of_hours_per_year + t].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
        }

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int r = 0; r < no_of_reservoirs; ++r) {
                nu_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                nu_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r].set (GRB_DoubleAttr_UB, 0.0);

                nu_hat_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                nu_hat_minus[(s * no_of_hours_per_year + t) * no_of_reservoirs + r].set (GRB_DoubleAttr_UB, 0.0);
            }
    
    for (int r = 0; r < no_of_reservoirs; ++r) {
        lambda_0[r].set (GRB_DoubleAttr_UB, 0.0);
        lambda_D[r].set (GRB_DoubleAttr_UB, 0.0);
    }
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int d = 0; d < no_of_days_per_year; ++d)
            for (int r = 0; r < no_of_reservoirs; ++r)
                lambda[(s * no_of_days_per_year + d) * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
    
    for (int d = 0; d < no_of_days_per_year; ++d)
        for (int r = 0; r < no_of_reservoirs; ++r) {
            for (int rr = 0; rr < no_of_reservoirs; ++rr) {
                Phi[(d * no_of_reservoirs + r) * no_of_reservoirs + rr].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                Psi[(d * no_of_reservoirs + r) * no_of_reservoirs + rr].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            }
            
            epsilon[d * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
            theta[d * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
        }
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int d = 0; d < no_of_days_per_year; ++d)
            for (int r = 0; r < no_of_reservoirs; ++r) {
                mu_minus[(s * no_of_days_per_year + d) * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, -GRB_INFINITY);
                mu_minus[(s * no_of_days_per_year + d) * no_of_reservoirs + r].set (GRB_DoubleAttr_UB, 0.0);
            }
 
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int r = 0; r < no_of_reservoirs; ++r) {
            mu_plus[(s * no_of_days_per_year + 0) * no_of_reservoirs + r].set (GRB_DoubleAttr_UB, 0.0);
            mu_minus[(s * no_of_days_per_year + 0) * no_of_reservoirs + r].set (GRB_DoubleAttr_LB, 0.0);
        }
    
    cout << "done: " << toc() << " seconds." << endl;
    
    // set up constraints: suv-constraints
    tic();
    cout << "Setting up suv constraints..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t) {
            model.addConstr (kappa[s * no_of_hours_per_year + t] +
                             kappa_hat[s * no_of_hours_per_year + t] -
                             kappa_prime[s * no_of_hours_per_year + t] == spot_price__in_sample[s][t]);
            if (include_reserve) {
                model.addConstr (kappa[s * no_of_hours_per_year + t] +
                                 prob_of_reserve_up_call * kappa_hat[s * no_of_hours_per_year + t] >=
                                 prob_of_reserve_up_call * reserve_up_price__in_sample[s][t]);
                model.addConstr (kappa_prime[s * no_of_hours_per_year + t] -
                                 prob_of_reserve_down_call * kappa_hat[s * no_of_hours_per_year + t] >=
                                 prob_of_reserve_down_call * reserve_down_price__in_sample[s][t]);
            }
        }

    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: gp-constraints
    tic();
    cout << "Setting up gp constraints..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int a = 0; a < no_of_arcs; ++a) {
                model.addConstr (alpha[(s * no_of_hours_per_year + t) * no_of_arcs + a] -
                                 gamma[(s * no_of_hours_per_year + t) * no_of_arcs + a] -
                                 generator_efficiency[a] * kappa[s * no_of_hours_per_year + t] >= 0.0);
                model.addConstr (alpha_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a] -
                                 gamma_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a] -
                                 generator_efficiency[a] * kappa_hat[s * no_of_hours_per_year + t] >= 0.0);

                model.addConstr (beta[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                 gamma[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                 inverse_pump_efficiency[a] * kappa[s * no_of_hours_per_year + t] >= 0.0);
                model.addConstr (beta_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                 gamma_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a] +
                                 inverse_pump_efficiency[a] * kappa_hat[s * no_of_hours_per_year + t] >= 0.0);
            }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: f-constraints
    tic();
    cout << "Setting up f constraints..." << flush;
    
    GRBLinExpr expr2;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int a = 0; a < no_of_arcs; ++a) {
                expr = gamma[(s * no_of_hours_per_year + t) * no_of_arcs + a];
                for (int r = 0; r < no_of_reservoirs; ++r)
                    for (int tt = t; tt <= last_hour_of_day (day_of_hour (t)); ++tt)
                        expr += topology_matrix[r][a] *
                                (nu_minus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r] +
                                 nu_plus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r]);
                model.addConstr (expr == 0.0);

                expr2 = 0;
                
                expr = gamma_hat[(s * no_of_hours_per_year + t) * no_of_arcs + a];
                for (int r = 0; r < no_of_reservoirs; ++r) {
                    for (int tt = t; tt <= last_hour_of_day (day_of_hour (t)); ++tt)
                        expr += topology_matrix[r][a] *
                                (nu_hat_minus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r] +
                                 nu_hat_plus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r]);
                    
                    if (day_of_hour (t) < no_of_days_per_year - 1)
                        expr2 += topology_matrix[r][a] * lambda[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r];
                    else
                        expr2 += topology_matrix[r][a] * lambda_D[r];
                }
                
                model.addConstr (expr == expr2);
            }

    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: w-constraints
    tic();
    cout << "Setting up w constraints..." << flush;

    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int r = 0; r < no_of_reservoirs; ++r) {
                expr = mu_plus[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r] +
                       mu_minus[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r];
                
                for (int tt = first_hour_of_day (day_of_hour (t)); tt <= last_hour_of_day (day_of_hour (t)); ++tt)
                    expr += nu_minus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r] +
                            nu_plus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r] +
                            nu_hat_minus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r] +
                            nu_hat_plus[(s * no_of_hours_per_year + tt) * no_of_reservoirs + r];
                
                if (day_of_hour (t) == 0)
                    model.addConstr (expr == lambda[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r] -
                                             lambda_0[r]);
                else if (day_of_hour (t) == no_of_days_per_year - 1)
                    model.addConstr (expr == lambda_D[r] -
                                             lambda[(s * no_of_days_per_year + day_of_hour (t) - 1) * no_of_reservoirs + r]);
                else
                    model.addConstr (expr == lambda[(s * no_of_days_per_year + day_of_hour (t)) * no_of_reservoirs + r] -
                                             lambda[(s * no_of_days_per_year + day_of_hour (t) - 1) * no_of_reservoirs + r]);
            }
    
    cout << "done: " << toc() << " seconds." << endl;

    // set up constraints: LDR constraints
    tic();
    cout << "Setting up LDR constraints..." << flush;
    
    for (int s = 0; s < no_of_insample_scenarios__model; ++s)
        for (int d = 0; d < no_of_days_per_year; ++d)
            for (int r = 0; r < no_of_reservoirs; ++r) {
                expr = epsilon[d * no_of_reservoirs + r];
                
                for (int rr = 0; rr < no_of_reservoirs; ++rr) {
                    double coeff__cum_inflow_curr_day = 0.0;
                    for (int t = first_hour_of_day (d); t <= last_hour_of_day (d); ++t)
                        coeff__cum_inflow_curr_day += inflow__in_sample[s][rr][t];
                    expr += Phi[(d * no_of_reservoirs + r) * no_of_reservoirs + rr] * coeff__cum_inflow_curr_day;
                                        
                    double coeff__cum_inflow_past = 0.0;
                    if (d > 0)
                        for (int t = 0; t <= last_hour_of_day (d - 1); ++t)
                            coeff__cum_inflow_past += inflow__in_sample[s][rr][t];
                    expr += Psi[(d * no_of_reservoirs + r) * no_of_reservoirs + rr] * coeff__cum_inflow_curr_day;
                }

                double coeff__avg_spot_curr_day = 0.0;
                for (int t = first_hour_of_day (d); t <= last_hour_of_day (d); ++t)
                    coeff__avg_spot_curr_day += spot_price__in_sample[s][t] / no_of_hours_per_day;
                expr += theta[d * no_of_reservoirs + r] * coeff__avg_spot_curr_day;
                
                model.addConstr (lambda[(s * no_of_days_per_year + d) * no_of_reservoirs + r] == expr);
            }

    cout << "done: " << toc() << " seconds." << endl;

    // set up final water value constraints if desired (dual pendant to the primal cycle constraints)
    /*if (final_water_values != nullptr)
        for (int r = 0; r < no_of_reservoirs; ++r)
            model.addConstr (lambda_D[r] == (*final_water_values)[r]);*/
    
    // solve the problem
    model.optimize();
    
    if (model.get (GRB_IntAttr_Status) != GRB_OPTIMAL) {
        model.write ("error.lp");
        error ("solve_dual_master_problem", "Optimization problem not solved to optimality.");
    }
    
    double obj_value = model.get(GRB_DoubleAttr_ObjVal);
    
    if (water_values != nullptr)
        for (int r = 0; r < no_of_reservoirs; ++r)
            (*water_values)[r] = lambda[(0 * no_of_days_per_year + 0) * no_of_reservoirs + r].get (GRB_DoubleAttr_X);

    // de-allocate variables
    delete []alpha;
    delete []alpha_hat;
    
    delete []beta;
    delete []beta_hat;
    
    delete []gamma;
    delete []gamma_hat;

    delete []kappa;
    delete []kappa_hat;
    delete []kappa_prime;

    delete []nu_plus;
    delete []nu_hat_plus;
    delete []nu_minus;
    delete []nu_hat_minus;

    delete []lambda_0;
    delete []lambda_D;
    delete []lambda;
    
    delete []Phi;
    delete []Psi;
    delete []epsilon;
    delete []theta;
    
    delete []mu_plus;
    delete []mu_minus;

    // that's it!
    return obj_value;
}
