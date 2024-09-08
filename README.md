# LKSR TTP solver
A TTP solver based on expanding possible tour changes from simple 2-opt moves to longer moves selected by the lin-kernighan algorithm

## To build
Go to the LINKERN folder `cd LINKERN`
Make the executable `make LKSR`

## To run
Run the executable `./LKSR instance_name.ttp random_seed`
The resulting output file will be located at "ttpinstance.ttp.LKSR.systemtime" in the working directory

Example usage: `./LKSR a280_n279_bounded-strongly-corr_01.ttp 24` which will read the ttp instance in the `a280_n279_bounded-strongly-corr_01.ttp` file and output to, for example, `a280_n279_bounded-strongly-corr_01.ttp.LKSR.1582756281`