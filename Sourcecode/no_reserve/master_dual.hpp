//
//  master_dual.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef master_dual_hpp
#define master_dual_hpp

#include <vector>

double solve_dual_master_problem__stochastic (bool include_reserve, std::vector<double> *water_values = nullptr);

#endif /* master_dual_hpp */
