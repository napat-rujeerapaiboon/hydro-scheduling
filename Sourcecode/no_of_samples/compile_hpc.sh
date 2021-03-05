rm ./a.out
g++ *.cpp -O3 -std=c++11 -I /apps/gurobi/9.0.1/linux64/include -L /apps/gurobi/9.0.1/linux64/lib -lgurobi_c++ -lgurobi90
