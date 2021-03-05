//
//  input_data_io.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef input_data_io_hpp
#define input_data_io_hpp

// generate the input data -- only to be run once for a given problem instance
void generate_data (int curr_out_of_sample_no);

// read the input data
void read_data (const char *curr_dir);

// conduct some basic sensibility checks on the data
void verify_data();

#endif /* input_data_io_hpp */
