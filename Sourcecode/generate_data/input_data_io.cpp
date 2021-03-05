//
//  input_data_io.cpp
//  Hydroscheduling
//
//  Created by Wolfram Wiesemann on 01/12/2020.
//  Copyright © 2020 Wolfram Wiesemann. All rights reserved.
//

#include "input_data_io.hpp"

#include <random>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "auxiliary.hpp"
#include "data_and_parameters.hpp"

using namespace std;

// generates 2(!) years of inflow data, so that for each day in the year we can look forward for an entire year
void generate_inflows (double a_squared, const vector<double> &b_squared, const vector<double> &mu, vector<double> &inflows) {
    inflows.clear();
    inflows.resize (2 * no_of_hours_per_year);

    default_random_engine generator;
    
    normal_distribution<double> distribution_alpha (1.0, sqrt (a_squared));
    double alpha = distribution_alpha (generator);

    for (int d = 0; d < 2 * no_of_days_per_year; ++d) {
        normal_distribution<double> distribution_beta (0.0, sqrt (b_squared[d % no_of_days_per_year]));
        double beta = distribution_beta (generator);
        
        double inflow_for_current_day = max (alpha * mu[d % no_of_days_per_year] + beta, 0.0) / 1000.0;
        
        for (int t = 0; t < no_of_hours_per_day; ++t)
            inflows[d * no_of_hours_per_day + t] = inflow_for_current_day / no_of_hours_per_day;
    }
}

// generate the input data -- only to be run once for a given problem instance
void generate_data (int curr_out_of_sample_no) {
    tic();
    
    /*
     ***************************
     *** DIRECTORY MANAGMENT ***
     ***************************
     */

    tic();
    cout << "Directory management..." << flush;
    
    stringstream command;
    command << "rm -rf \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "\"";
    system (command.str().c_str());

    command.str (string());
    command << "mkdir \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "\"";
    system (command.str().c_str());

    for (int d = 0; d < no_of_days_per_year; ++d) {
        command.str (string());
        command << "mkdir \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "\"";
        system (command.str().c_str());
    }

    cout << "done: " << toc() << " seconds." << endl;

    /*
     ***********************
     *** SPOT PRICE DATA ***
     ***********************
     */
    
    // *** IN-SAMPLE SPOT PRICE DATA ***

    // open the in-sample spot price file corresponding to curr_out_of_sample_no
    stringstream filename;
    filename << directory << "/Raw Input Data/pi_is_";
    filename << setfill ('0') << setw (4) << curr_out_of_sample_no + 1;
    filename << ".txt";
    
    ifstream fin_in (filename.str().c_str());
    if (!fin_in.is_open()) error ("generate_data", "In-sample spot price file not found.");

    // read out and write all in-sample spot prices
    for (int d = 0; d < no_of_days_per_year; ++d) {
        tic();
        cout << "In-sample spot prices for day " << d << "..." << flush;

        for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
            filename.str (string());
            filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_spot_prices.insample_" << s << ".txt";
            ofstream fout (filename.str().c_str());

            for (int t = 0; t < no_of_hours_per_year; ++t) {
                double v;
                fin_in >> v;
                fout << v << endl;
            }
        }
        
        cout << "done: " << toc() << " seconds." << endl;
    }
    
    // *** OUT-OF-SAMPLE SPOT PRICE DATA ***
    
    // open the out-of-sample spot price file corresponding to curr_out_of_sample_no
    filename.str (string());
    filename << directory << "/Raw Input Data/pi_os_";
    filename << setfill ('0') << setw (4) << curr_out_of_sample_no + 1;
    filename << ".txt";
    
    ifstream fin_out (filename.str().c_str());
    if (!fin_out.is_open()) error ("generate_data", "Out-of-sample spot price file not found.");

    // read out and write all out-of-sample spot prices
    for (int d = 0; d < no_of_days_per_year; ++d) {
        tic();
        cout << "Out-of-sample spot prices for day " << d << "..." << flush;

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_spot_prices.outofsample.txt";
        ofstream fout (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t) {
            double v;
            fin_out >> v;
            fout << v << endl;
        }
        
        cout << "done: " << toc() << " seconds." << endl;
    }

    /*
     **************************
     *** RESERVE PRICE DATA ***
     **************************
     */

    const double reserve_up_price__peak   =  737.98, reserve_up_price__offpeak   =  802.40;     // in EUR/MWh
    const double reserve_down_price__peak = 1898.03, reserve_down_price__offpeak = 1983.30;     // in EUR/MWh

    // *** IN-SAMPLE RESERVE PRICE DATA ***
    
    for (int d = 0; d < no_of_days_per_year; ++d) {
        tic();
        cout << "In-sample reserve prices for day " << d << "..." << flush;

        for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
            filename.str (string());
            filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_up_prices.insample_" << s << ".txt";
            ofstream fout_up (filename.str().c_str());

            filename.str (string());
            filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_down_prices.insample_" << s << ".txt";
            ofstream fout_down (filename.str().c_str());
            
            for (int ww = 0; ww < no_of_weeks_per_year; ++ww)
                for (int dd = 0; dd < no_of_days_per_week; ++dd)
                    for (int hh = 0; hh < no_of_hours_per_day; ++hh) {
                        bool weekday = (((d + dd) % no_of_days_per_week) < 5);
                        bool weekday_peak = (hh > 7) && (hh < 20);

                        bool peak = (weekday && weekday_peak);

                        if (peak) {
                            fout_up << reserve_up_price__peak << endl;
                            fout_down << reserve_down_price__peak << endl;
                        } else {
                            fout_up << reserve_up_price__offpeak << endl;
                            fout_down << reserve_down_price__offpeak << endl;
                        }
                    }
        }
        
        cout << "done: " << toc() << " seconds." << endl;
    }

    // *** OUT-OF-SAMPLE RESERVE PRICE DATA ***
    
    for (int d = 0; d < no_of_days_per_year; ++d) {
        tic();
        cout << "Out-of-sample reserve prices for day " << d << "..." << flush;

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_up_prices.outofsample.txt";
        ofstream fout_up (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_down_prices.outofsample.txt";
        ofstream fout_down (filename.str().c_str());
                      
        for (int hh = 0; hh < no_of_hours_per_day; ++hh) {
            bool weekday = ((d % no_of_days_per_week) < 5);
            bool weekday_peak = (hh > 7) && (hh < 20);

            bool peak = (weekday && weekday_peak);

            if (peak) {
                fout_up << reserve_up_price__peak << endl;
                fout_down << reserve_down_price__peak << endl;
            } else {
                fout_up << reserve_up_price__offpeak << endl;
                fout_down << reserve_down_price__offpeak << endl;
            }
        }
    
        cout << "done: " << toc() << " seconds." << endl;
    }

    /*
     *********************
     *** RESERVE CALLS ***
     *********************
     */

    tic();
    cout << "Reserve calls..." << flush;

    for (int d = 0; d < no_of_days_per_year; ++d) {
        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_up_calls.outofsample.txt";
        ofstream fout_up (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_reserve_down_calls.outofsample.txt";
        ofstream fout_down (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t) {
            double v = (double)(rand()) / (double)(RAND_MAX);
            fout_up << (v <= prob_of_reserve_up_call) << endl;
            fout_down << (v >= 1.0 - prob_of_reserve_down_call) << endl;
        }
    }

    cout << "done: " << toc() << " seconds." << endl;

    /*
     *******************
     *** INFLOW DATA ***
     *******************
     */

    // *** READ INFLOW MODEL PARAMETERS ***
    
    // open files
    filename.str (string());
    filename << directory << "/Raw Input Data/Bockhartsee.txt";
    ifstream fin_bock (filename.str().c_str());
    if (!fin_bock.is_open()) error ("generate_data", "Inflow data file not found.");

    filename.str (string());
    filename << directory << "/Raw Input Data/Nassfeld.txt";
    ifstream fin_nass (filename.str().c_str());
    if (!fin_nass.is_open()) error ("generate_data", "Inflow data file not found.");

    // read a^2
    double bock__a_squared; fin_bock >> bock__a_squared;
    double nass__a_squared; fin_nass >> nass__a_squared;

    // read b^2
    vector<double> bock__b_squared, nass__b_squared;
    bock__b_squared.resize (no_of_days_per_year);
    nass__b_squared.resize (no_of_days_per_year);
    for (int d = 0; d < no_of_days_per_year; ++d) {
        fin_bock >> bock__b_squared[d];
        fin_nass >> nass__b_squared[d];
    }
 
    // read mu
    vector<double> bock__mu, nass__mu;
    bock__mu.resize (no_of_days_per_year);
    nass__mu.resize (no_of_days_per_year);
    for (int d = 0; d < no_of_days_per_year; ++d) {
        fin_bock >> bock__mu[d];
        fin_nass >> nass__mu[d];
    }

    // *** GENERATE OUT-OF-SAMPLE INFLOWS ***

    tic();
    cout << "Out-of-sample inflows..." << flush;

    // out-of-sample inflows start at day 0 of the year, looking one year ahead
    vector<double> hourly_inflows__bock__out_of_sample;         // [2 * no_of_hours_per_year]
    vector<double> hourly_inflows__nass__out_of_sample;         // [2 * no_of_hours_per_year]

    // generate data
    generate_inflows (bock__a_squared, bock__b_squared, bock__mu, hourly_inflows__bock__out_of_sample);
    generate_inflows (nass__a_squared, nass__b_squared, nass__mu, hourly_inflows__nass__out_of_sample);

    // write data
    for (int d = 0; d < no_of_days_per_year; ++d) {
        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_inflows.outofsample.txt";
        ofstream fout (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t)
            fout << hourly_inflows__bock__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                 << hourly_inflows__nass__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                 << "0.0" << endl;
    }
    
    cout << "done: " << toc() << " seconds." << endl;

    // *** GENERATE IN-SAMPLE INFLOWS ***

    // in-sample inflows start at the current day of the year, looking one year ahead from the current day
    // the first 24 hours of in-sample inflows [0]...[no_of_hours_per_day - 1] coincide with corresponding out-of-sample
    // inflows [curr_day * no_of_hours_per_day + 0]...[curr_day * no_of_hours_per_day + no_of_hours_per_day - 1]
    vector<double> hourly_inflows__bock__in_sample;              // [2 * no_of_hours_per_year]
    vector<double> hourly_inflows__nass__in_sample;              // [2 * no_of_hours_per_year]

    for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
        tic();
        cout << "In-sample inflows for sample " << s << "..." << flush;

        // generate data
        generate_inflows (bock__a_squared, bock__b_squared, bock__mu, hourly_inflows__bock__in_sample);
        generate_inflows (nass__a_squared, nass__b_squared, nass__mu, hourly_inflows__nass__in_sample);

        // write data
        for (int d = 0; d < no_of_days_per_year; ++d) {
            // open file
            filename.str (string());
            filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "/hourly_inflows.insample_" << s << ".txt";
            ofstream fout (filename.str().c_str());
            
            for (int t = 0; t < no_of_hours_per_year; ++t) {
                if (t < no_of_hours_per_day)
                    fout << hourly_inflows__bock__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                         << hourly_inflows__nass__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                         << "0.0" << endl;
                else
                    fout << hourly_inflows__bock__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                         << hourly_inflows__nass__out_of_sample[d * no_of_hours_per_day + t] << "\t"
                         << "0.0" << endl;
            }
        }
    
        cout << "done: " << toc() << " seconds." << endl;
    }
    
    /*
     *************************
     *** ARCHIVE MANAGMENT ***
     *************************
     */

    tic();
    cout << "Archive management..." << flush;
    
    for (int d = 0; d < no_of_days_per_year; ++d) {
        command.str (string());
        command << "tar -zcf \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day_" << d << ".tar.gz\" "
                << "\"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "\"";
        system (command.str().c_str());
        command.str (string());
        command << "rm -rf \"" << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << d << "\"";
        system (command.str().c_str());
    }

    cout << "done: " << toc() << " seconds." << endl;

    
    // that's it!
    cout << "OVERALL GENERATION TIME: " << toc() << " seconds." << endl;
}

// read the input data
void read_data (int curr_out_of_sample_no, int curr_day) {
    tic();

    // read data: spot_price__in_sample, reserve_up_price__in_sample, reserve_down_price__in_sample
    tic();
    cout << "Reading in-sample spot and reserve prices..." << flush;

    spot_price__in_sample.resize (no_of_insample_scenarios__data);
    reserve_up_price__in_sample.resize (no_of_insample_scenarios__data);
    reserve_down_price__in_sample.resize (no_of_insample_scenarios__data);

    for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
        spot_price__in_sample[s].resize (no_of_hours_per_year);
        reserve_up_price__in_sample[s].resize (no_of_hours_per_year);
        reserve_down_price__in_sample[s].resize (no_of_hours_per_year);

        stringstream filename;
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_spot_prices.insample_" << s << ".txt";
        ifstream fin_spot (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_up_prices.insample_" << s << ".txt";
        ifstream fin_up (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_down_prices.insample_" << s << ".txt";
        ifstream fin_down (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_year; ++t) {
            fin_spot >> spot_price__in_sample[s][t];
            fin_up >> reserve_up_price__in_sample[s][t];
            fin_down >> reserve_down_price__in_sample[s][t];
        }
    }

    cout << "done: " << toc() << " seconds." << endl;

    // read data: spot_price__out_of_sample, reserve_up_price__out_of_sample, reserve_down_price__out_of_sample
    tic();
    cout << "Reading out-of-sample spot and reserve prices..." << flush;

    spot_price__out_of_sample.resize (no_of_hours_per_day);
    reserve_up_price__out_of_sample.resize (no_of_hours_per_day);
    reserve_down_price__out_of_sample.resize (no_of_hours_per_day);

    {
        stringstream filename;
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_spot_prices.outofsample.txt";
        ifstream fin_spot (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_up_prices.outofsample.txt";
        ifstream fin_up (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_down_prices.outofsample.txt";
        ifstream fin_down (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t) {
            fin_spot >> spot_price__out_of_sample[t];
            fin_up >> reserve_up_price__out_of_sample[t];
            fin_down >> reserve_down_price__out_of_sample[t];
        }
    }

    cout << "done: " << toc() << " seconds." << endl;

    // read data: reserve_up_call__out_of_sample, reserve_down_call__out_of_sample
    reserve_up_call__out_of_sample.resize (no_of_hours_per_day);
    reserve_down_call__out_of_sample.resize (no_of_hours_per_day);

    tic();
    cout << "Reading out-of-sample reserve calls..." << flush;

    {
        stringstream filename;
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_up_calls.outofsample.txt";
        ifstream fin_up (filename.str().c_str());

        filename.str (string());
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_reserve_down_calls.outofsample.txt";
        ifstream fin_down (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t) {
            int v;
            fin_up >> v; reserve_up_call__out_of_sample[t] = (v == 1);
            fin_down >> v; reserve_down_call__out_of_sample[t] = (v == 1);
        }
    }

    cout << "done: " << toc() << " seconds." << endl;

    // read data: inflow__in_sample
    tic();
    cout << "Reading in-sample inflows..." << flush;

    inflow__in_sample.resize (no_of_insample_scenarios__data);
    for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
        inflow__in_sample[s].resize (no_of_reservoirs);
        for (int r = 0; r < no_of_reservoirs; ++r)
            inflow__in_sample[s][r].resize (no_of_hours_per_year);
    }
    
    for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
        stringstream filename;
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_inflows.insample_" << s << ".txt";
        ifstream fin (filename.str().c_str());
        
        for (int t = 0; t < no_of_hours_per_year; ++t)
            for (int r = 0; r < no_of_reservoirs; ++r)
                fin >> inflow__in_sample[s][r][t];
    }

    cout << "done: " << toc() << " seconds." << endl;

    // read data: inflow__out_of_sample
    tic();
    cout << "Reading out-of-sample inflows..." << flush;

    inflow__out_of_sample.resize (no_of_reservoirs);
    for (int r = 0; r < no_of_reservoirs; ++r)
        inflow__out_of_sample[r].resize (no_of_hours_per_day);

    {
        stringstream filename;
        filename << directory << "/Processed Input Data/Outer Sample " << curr_out_of_sample_no << "/Day " << curr_day << "/hourly_inflows.outofsample.txt";
        ifstream fin (filename.str().c_str());

        for (int t = 0; t < no_of_hours_per_day; ++t)
            for (int r = 0; r < no_of_reservoirs; ++r)
                fin >> inflow__out_of_sample[r][t];
    }
    
    cout << "done: " << toc() << " seconds." << endl;

    // that's it!
    cout << "OVERALL READING TIME: " << toc() << " seconds." << endl;
}

inline bool approx_equal (double a, double b) { return (fabs (a - b) < 0.000001); }

// conduct some basic sensibility checks on the data
void verify_data() {
    // check that the first no_of_hours_per_day many in-sample prices are the same across all samples
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int s = 0; s < no_of_insample_scenarios__data; ++s) {
            if (!approx_equal (spot_price__in_sample[s][t], spot_price__in_sample[0][t]))
                error ("verify_data", "In-sample spot prices inconsistent.");
            if (!approx_equal (reserve_up_price__in_sample[s][t], reserve_up_price__in_sample[0][t]))
                error ("verify_data", "In-sample reserve-up prices inconsistent.");
            if (!approx_equal (reserve_down_price__in_sample[s][t], reserve_down_price__in_sample[0][t]))
                error ("verify_data", "In-sample reserve-down prices inconsistent.");
        }

    // check that the first no_of_hours_per_day many in-sample prices coincide with the out-of-sample prices
    for (int t = 0; t < no_of_hours_per_day; ++t) {
        if (!approx_equal (spot_price__in_sample[0][t], spot_price__out_of_sample[t]))
            error ("verify_data", "In-sample spot prices do not match out-of-sample prices.");
        if (!approx_equal (reserve_up_price__in_sample[0][t], reserve_up_price__out_of_sample[t]))
            error ("verify_data", "In-sample reserve-up prices do not match out-of-sample prices.");
        if (!approx_equal (reserve_down_price__in_sample[0][t], reserve_down_price__out_of_sample[t]))
            error ("verify_data", "In-sample reserve-down prices do not match out-of-sample prices.");
    }

    // check that the first no_of_hours_per_day many in-sample inflows are the same across all samples
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int r = 0; r < no_of_reservoirs; ++r)
            for (int s = 0; s < no_of_insample_scenarios__data; ++s)
                if (!approx_equal (inflow__in_sample[s][r][t], inflow__in_sample[0][r][t]))
                    error ("verify_data", "In-sample inflows inconsistent.");

    // check that the first no_of_hours_per_day many in-sample inflows coincide with the out-of-sample inflows
    for (int t = 0; t < no_of_hours_per_day; ++t)
        for (int r = 0; r < no_of_reservoirs; ++r)
            if (!approx_equal (inflow__in_sample[0][r][t], inflow__out_of_sample[r][t]))
                error ("verify_data", "In-sample inflows do not math out-of-sample inflows.");
}
