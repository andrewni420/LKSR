# LKSR TTP solver
A Traveling Thief Problem (TTP) solver for the 2024 GECCO TTP competition based on expanding the space of possible tour changes from simple 2-opt moves to longer moves selected by the Lin-Kernighan TSP algorithm

## To build
Go to the LINKERN folder `cd LINKERN`
Make the executable `make LKSR`

## To run
Run the executable `./LKSR instance_name.ttp random_seed`
The resulting output file will be located at "ttpinstance.ttp.LKSR.systemtime" in the working directory

Example usage: `./LKSR a280_n279_bounded-strongly-corr_01.ttp 24` which will read the ttp instance in the `a280_n279_bounded-strongly-corr_01.ttp` file and output to, for example, `a280_n279_bounded-strongly-corr_01.ttp.LKSR.1582756281`
