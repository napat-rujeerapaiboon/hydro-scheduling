//
//  results.hpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#ifndef results_hpp
#define results_hpp

#include <vector>
#include <fstream>

/*
   *******************************
   *** FILE FORMAT FOR RESULTS ***
   *******************************

   for each day:

   1) day # (first day starts with 0)

   2) for each of {primal water values}, {dual water values}, {mixed water values}:

       2)a) initial reservoir levels
 
       2)b) master problem optimal objective value
 
       2)c) water values for the current day

       2)d) for each hour of the current day:

           2)d)1) spot bid for that hour, spot revenues for that hour
           2)d)2) reserve-up bid for that hour, reserve-up revenues for that hour
           2)d)3) reserve-down bid for that hour, reserve-down revenues for that hour
           2)d)4) missed reserve-up commitment, missed reserve-down commitment for that hour
           2)d)5) generate (for all arcs), pump (for all arcs), spill (for all arcs)

       2)e) terminal reservoir levels
*/

// data structure for a single hour of a day; 2)d) above
struct tResults_HourlyRecord {
    double spot_bid,                        spot_revenues;
    double reserve_up_bid,                  reserve_up_revenues;
    double reserve_down_bid,                reserve_down_revenues;
    double missed_reverse_up_commitment,    missed_reverse_down_commitment;
    std::vector<double> generate, pump, spill;
    
    tResults_HourlyRecord() : spot_bid (-1.0), spot_revenues (-1.0),
                              reserve_up_bid (-1.0), reserve_up_revenues (-1.0),
                              reserve_down_bid (-1.0), reserve_down_revenues (-1.0),
                              missed_reverse_up_commitment (-1.0), missed_reverse_down_commitment (-1.0) {}
    
    void read (std::ifstream &fin);
    void write (std::ofstream &fout);
};

// data structure for a single day and method (i.e., primal/dual/mixed); 2) above
struct tResults_DailyRecord {
    std::vector<double> initial_reservoir_levels;                           // [no_of_reservoirs]
    std::vector<double> water_values;                                       // [no_of_reservoirs]
    std::vector<tResults_HourlyRecord> hourly_records;                      // [no_of_hours_per_day]
    std::vector<double> terminal_reservoir_levels;                           // [no_of_reservoirs]
    double master_objval;
    
    tResults_DailyRecord();

    void read (std::ifstream &fin);
    void write (std::ofstream &fout);
};

// data structure for all days in the result
struct tResultsFile {
    std::vector<std::vector<tResults_DailyRecord> > daily_records;          // [no_of_days_per_year][primal/dual/mixed]
    
    tResultsFile();

    int read (const char *filename);                                        // returns day that should be solved next
    void write (const char *filename, int curr_day);

    // {primal/dual/mixed}_dailyreservoirlevels.csv     -- rows = days, columns = reservoirs
    // {primal/dual/mixed}_dailywatervalues.csv         -- rows = days, columns = reservoirs
    // {primal/dual/mixed}_hourlybids.csv               -- rows = hours, columns = cumulative bids per market
    // {primal/dual/mixed}_hourlyrevenues.csv           -- rows = hours, columns = cumulative revenues per market (spot/up/down),
    //                                                  -- joint cumulative revenues
    void write_tables (int curr_day);
};


#endif /* results_hpp */
