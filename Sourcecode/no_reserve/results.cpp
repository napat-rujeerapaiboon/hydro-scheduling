//
//  results.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "results.hpp"

#include <sstream>

#include "auxiliary.hpp"
#include "data_and_parameters.hpp"

using namespace std;

void tResults_HourlyRecord::read (ifstream &fin) {
    fin >> spot_bid;                        fin >> spot_revenues;
    fin >> reserve_up_bid;                  fin >> reserve_up_revenues;
    fin >> reserve_down_bid;                fin >> reserve_down_revenues;
    fin >> missed_reverse_up_commitment;    fin >> missed_reverse_down_commitment;
    
    generate.resize (no_of_arcs);
    for (int a = 0; a < no_of_arcs; ++a)
        fin >> generate[a];

    pump.resize (no_of_arcs);
    for (int a = 0; a < no_of_arcs; ++a)
        fin >> pump[a];

    spill.resize (no_of_arcs);
    for (int a = 0; a < no_of_arcs; ++a)
        fin >> spill[a];
}

void tResults_HourlyRecord::write (ofstream &fout) {
    fout << "\t\t";
    fout << spot_bid << " ";                        fout << spot_revenues << "\t";
    fout << reserve_up_bid << " ";                  fout << reserve_up_revenues << "\t";
    fout << reserve_down_bid << " ";                fout << reserve_down_revenues << "\t";
    fout << missed_reverse_up_commitment << " ";    fout << missed_reverse_down_commitment << endl;
    
    fout << "\t\t";
    for (int a = 0; a < no_of_arcs; ++a) {
        fout << generate[a];
        if (a < no_of_arcs - 1) fout << " "; else fout << "\t";
    }

    for (int a = 0; a < no_of_arcs; ++a) {
        fout << pump[a];
        if (a < no_of_arcs - 1) fout << " "; else fout << "\t";
    }

    for (int a = 0; a < no_of_arcs; ++a) {
        fout << spill[a];
        if (a < no_of_arcs - 1) fout << " "; else fout << endl;
    }
}

tResults_DailyRecord::tResults_DailyRecord() {
    initial_reservoir_levels.resize (no_of_reservoirs, -1.0);
    master_objval = -1.0;
    water_values.resize (no_of_reservoirs, -1.0);
    hourly_records.resize (no_of_hours_per_day);
    terminal_reservoir_levels.resize (no_of_reservoirs, -1.0);
}

void tResults_DailyRecord::read (ifstream &fin) {
    for (int r = 0; r < no_of_reservoirs; ++r)
        fin >> initial_reservoir_levels[r];

    fin >> master_objval;
    
    for (int r = 0; r < no_of_reservoirs; ++r)
        fin >> water_values[r];

    for (int t = 0; t < no_of_hours_per_day; ++t)
        hourly_records[t].read (fin);

    for (int r = 0; r < no_of_reservoirs; ++r)
        fin >> terminal_reservoir_levels[r];
}

void tResults_DailyRecord::write (ofstream &fout) {
    fout << "\t";
    for (int r = 0; r < no_of_reservoirs; ++r)
        fout << initial_reservoir_levels[r] << " ";
    fout << endl;
    
    fout << "\t" << master_objval << endl;

    fout << "\t";
    for (int r = 0; r < no_of_reservoirs; ++r)
        fout << water_values[r] << " ";
    fout << endl;

    for (int t = 0; t < no_of_hours_per_day; ++t)
        hourly_records[t].write (fout);

    fout << "\t";
    for (int r = 0; r < no_of_reservoirs; ++r)
        fout << terminal_reservoir_levels[r] << " ";
    fout << endl;
}


tResultsFile::tResultsFile() {
    daily_records.resize (no_of_days_per_year);
}

int tResultsFile::read (const char *filename) {                         // returns day that should be solved next
    ifstream fin (filename);
    if (!fin.is_open())
        return 0;

    int curr_day = 0;
    for (;;) {
        int v; fin >> v;
        if (fin.eof()) break;
        if (v != curr_day) error ("tResultsFile::read", "Expected to read number of day.");
        
        daily_records[curr_day].read (fin);
        
        ++curr_day;
    }
    
    return curr_day;
}

void tResultsFile::write (const char *filename, int curr_day) {
    ofstream fout (filename);
    
    for (int d = 0; d <= curr_day; ++d) {
        fout << d << endl;
        daily_records[d].write (fout);
        fout << endl;
    }
}

#if 0

//const string method_name[3] = {"primal", "dual", "mixed"};
//const string method_name[3] = {"1", "5", "10"};
const string method_name[3] = {"dd", "ds", "sd"};

// {method[*]}_dailyreservoirlevels.csv      -- rows = days, columns = reservoirs
// {method[*]}_dailywatervalues.csv          -- rows = days, columns = reservoirs
// {methods[*]}_hourlybids.csv               -- rows = hours, columns = cumulative bids per market
// {methods[*]}_hourlyrevenues.csv           -- rows = hours, columns = cumulative revenues per market (spot/up/down),
//                                           -- joint cumulative revenues
// {methods[*]}_missedreservecommitments.csv -- rows = hours, columns = cumulative reserve-up commitments
//                                           -- cumulative missed reserve-up commitments
//                                           -- cumulative reserve-down commitments
//                                           -- cumulative missed reserve-down commitments
// {methods[*]}_missedreservecommitments.csv -- rows = hours, columns = cumulative reserve-up commitments
//                                           -- cumulative missed reserve-up commitments
//                                           -- cumulative reserve-down commitments
//                                           -- cumulative missed reserve-down commitments
// {methods[*]}_hydraulics.csv               -- rows = hours, columns = generate[no_of_arcs], pump[no_of_arcs], spill[no_of_arcs]
void tResultsFile::write_tables (int curr_day) {
    /*
     ******************************************************************************************************************************
     {methods[*]}_dailyreservoirlevels.csv     -- rows = days, columns = reservoirs
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_dailyreservoirlevels.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Day #";
        for (int r = 0; r < no_of_reservoirs; ++r)
            fout << ", Initial daily water level of reservoir " << r;
        fout << endl;
            
        for (int d = 0; d <= curr_day; ++d) {
            fout << d;
            for (int r = 0; r < no_of_reservoirs; ++r)
                fout << ", " << daily_records[d][method].initial_reservoir_levels[r];
            fout << endl;
        }
    }
        
    /*
     ******************************************************************************************************************************
     {methods[*]}_dailywatervalues.csv         -- rows = days, columns = reservoirs
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_dailywatervalues.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Day #";
        for (int r = 0; r < no_of_reservoirs; ++r)
            fout << ", Water value of reservoir " << r;
        fout << endl;
            
        for (int d = 0; d <= curr_day; ++d) {
            fout << d;
            for (int r = 0; r < no_of_reservoirs; ++r)
                fout << ", " << daily_records[d][method].water_values[r];
            fout << endl;
        }
    }

    /*
     ******************************************************************************************************************************
     {methods[*]}_hourlybids.csv               -- rows = hours, columns = cumulative bids per market
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_hourlybids.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Hour #, Cumulative bids on spot market, Cumulative bids on reserve-up market, Cumulative bids on reserve-down market" << endl;
        
        double cum_spot = 0.0, cum_up = 0.0, cum_down = 0.0;
        
        for (int d = 0; d <= curr_day; ++d)
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                cum_spot += daily_records[d][method].hourly_records[t].spot_bid;
                cum_up += daily_records[d][method].hourly_records[t].reserve_up_bid;
                cum_down += daily_records[d][method].hourly_records[t].reserve_down_bid;
                
                fout << d * no_of_hours_per_day + t << ", "
                     << cum_spot << ", "
                     << cum_up << ", "
                     << cum_down << endl;
            }
    }
    
    /*
     ******************************************************************************************************************************
     {methods[*]}_hourlyrevenues.csv           -- rows = hours, columns = cumulative revenues per market (spot/up/down),
                                               -- joint cumulative revenues
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_hourlyrevenues.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Hour #, Cumulative revenues on spot market, Cumulative revenues on reserve-up market, Cumulative revenues on reserve-down market, Cumulative revenues on all markets" << endl;
        
        double cum_spot = 0.0, cum_up = 0.0, cum_down = 0.0;
        
        for (int d = 0; d <= curr_day; ++d)
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                cum_spot += daily_records[d][method].hourly_records[t].spot_revenues;
                cum_up += daily_records[d][method].hourly_records[t].reserve_up_revenues;
                cum_down += daily_records[d][method].hourly_records[t].reserve_down_revenues;
                
                fout << d * no_of_hours_per_day + t << ", "
                     << cum_spot << ", "
                     << cum_up << ", "
                     << cum_down << ", "
                     << cum_spot + cum_up + cum_down << endl;
            }
    }
    
    /*
     ******************************************************************************************************************************
     {methods[*]}_missedreservecommitments.csv -- rows = hours, columns = cumulative reserve-up commitments,
                                               -- cumulative missed reserve-up commitments
                                               -- cumulative reserve-down commitments
                                               -- cumulative missed reserve-down commitments
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_missedreservecommitments.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Hour #, Cumulative reserve-up commitments, Cumulative missed reserve-up commitments, Cumulative reserve-down commitments, Cumulative missed reserve-down commitments" << endl;
        
        double cum_reserve_up = 0.0, cum_reserve_up_missed = 0.0;
        double cum_reserve_down = 0.0, cum_reserve_down_missed = 0.0;

        for (int d = 0; d <= curr_day; ++d)
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                cum_reserve_up += daily_records[d][method].hourly_records[t].reserve_up_bid;
                cum_reserve_up_missed += daily_records[d][method].hourly_records[t].missed_reverse_up_commitment;

                cum_reserve_down += daily_records[d][method].hourly_records[t].reserve_down_bid;
                cum_reserve_down_missed += daily_records[d][method].hourly_records[t].missed_reverse_down_commitment;
                
                fout << d * no_of_hours_per_day + t << ", "
                     << cum_reserve_up << ", "
                     << cum_reserve_up_missed << ", "
                     << cum_reserve_down << ", "
                     << cum_reserve_down_missed << endl;
            }
    }
    
    /*
     ******************************************************************************************************************************
     {methods[*]}_hydraulics.csv               -- rows = hours, columns = generate[no_of_arcs], pump[no_of_arcs], spill[no_of_arcs]
     ******************************************************************************************************************************
     */

    for (int method = 0; method < 3; ++method){
        stringstream filename;
        filename << directory << "/Tables/" << method_name[method] << "_hydraulics.csv";
        ofstream fout (filename.str().c_str());
        
        fout << "Hour #, ";
        for (int a = 0; a < no_of_arcs; ++a) fout << "Generation on Arc #" << a << ", ";
        for (int a = 0; a < no_of_arcs; ++a) fout << "Pumping on Arc #" << a << ", ";
        for (int a = 0; a < no_of_arcs; ++a) {
            fout << "Spilling on Arc #" << a;
            if (a < no_of_arcs - 1) fout << ", "; else fout << endl;
        }
        
        for (int d = 0; d <= curr_day; ++d)
            for (int t = 0; t < no_of_hours_per_day; ++t) {
                fout << d * no_of_hours_per_day + t << ", ";
                for (int a = 0; a < no_of_arcs; ++a) fout << daily_records[d][method].hourly_records[t].generate[a] << ", ";
                for (int a = 0; a < no_of_arcs; ++a) fout << daily_records[d][method].hourly_records[t].pump[a] << ", ";
                for (int a = 0; a < no_of_arcs; ++a) {
                    fout << daily_records[d][method].hourly_records[t].spill[a];
                    if (a < no_of_arcs - 1) fout << ", "; else fout << endl;
                }
            }
    }
}

#endif
