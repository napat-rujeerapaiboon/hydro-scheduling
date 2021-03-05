rm ./a.out
g++ -O3 -I /Library/gurobi902/mac64/include -L /Library/gurobi902/mac64/lib -lgurobi_c++ -lgurobi90 -Ofast *.cpp
