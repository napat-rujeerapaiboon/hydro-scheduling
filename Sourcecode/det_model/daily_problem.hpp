//
//  daily_problem.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef daily_problem_hpp
#define daily_problem_hpp

#include <vector>

double solve_first_hour_of_daily_problem__stochastic (const std::vector<double> &water_values, bool include_reserve, std::vector<double> *result_s = nullptr, std::vector<double> *result_u = nullptr, std::vector<double> *result_v = nullptr);

double solve_first_hour_of_daily_problem__deterministic__robust_bounds (const std::vector<double> &water_values, bool include_reserve, std::vector<double> *result_s = nullptr, std::vector<double> *result_u = nullptr, std::vector<double> *result_v = nullptr);

double solve_first_hour_of_daily_problem__deterministic__percentage_of_spot (const std::vector<double> &water_values, bool include_reserve, std::vector<double> *result_s = nullptr, std::vector<double> *result_u = nullptr, std::vector<double> *result_v = nullptr);

double solve_each_hour_of_daily_problem__stochastic (int theta, const std::vector<double> &water_values, const std::vector<double> &bid__spot, const std::vector<double> &bid__reserve_up, const std::vector<double> &bid__reserve_down, std::vector<double> &res_generate, std::vector<double> &res_pump, std::vector<double> &res_spill);

double solve_each_hour_of_daily_problem__deterministic__robust_bounds (int theta, const std::vector<double> &water_values, const std::vector<double> &bid__spot, const std::vector<double> &bid__reserve_up, const std::vector<double> &bid__reserve_down, std::vector<double> &res_generate, std::vector<double> &res_pump, std::vector<double> &res_spill, double *missed_commitment_u, double *missed_commitment_v);

double solve_each_hour_of_daily_problem__deterministic__percentage_of_spot (int theta, const std::vector<double> &water_values, const std::vector<double> &bid__spot, const std::vector<double> &bid__reserve_up, const std::vector<double> &bid__reserve_down, std::vector<double> &res_generate, std::vector<double> &res_pump, std::vector<double> &res_spill);

#endif /* daily_problem_hpp */
