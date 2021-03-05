//
//  master_primal.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef master_primal_hpp
#define master_primal_hpp

#include <vector>

double solve_primal_master_problem__stochastic (bool include_reserve, std::vector<double> *water_values = nullptr);
double solve_primal_master_problem__deterministic__robust_bounds (bool include_reserve, std::vector<double> *water_values = nullptr);
double solve_primal_master_problem__deterministic__percentage_of_spot (bool include_reserve, std::vector<double> *water_values = nullptr);

#endif /* master_primal_hpp */
