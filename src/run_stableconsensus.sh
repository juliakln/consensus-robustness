#!/bin/bash

# examples
python analysis_stableconsensus.py stubborns_parameters \
    --model voter \
    --N 100 \
    --stubborn_type z \
    --plot

python analysis_stableconsensus.py stubborns_groupsizes \
    --model voter \
    --group_sizes 20 50 100 1000 \
    --stubborn_type z \
    --plot

python analysis_stableconsensus.py xy_stubborns \
    --model voter \
    --N 100 \
    --stubborn_type z \
    --ratex 1.05 \
    --ratey 0.95 \
    --plot


## Exploring Consensus Robustness in Swarms with Disruptive Individuals, ISoLA 2024
#Figure 4(a)
python analysis_stableconsensus.py stubborns_parameters \
    --model crossinh \
    --N 100 \
    --stubborn_type z \
    --plot

#Figure 5(a)
python analysis_stableconsensus.py stubborns_parameters \
    --model crossinh \
    --N 100 \
    --stubborn_type c \
    --plot

#Figure 7(a)
python analysis_stableconsensus.py stubborns_parameters \
    --model crossinh \
    --N 500 \
    --stubborn_type z \
    --plot

#Figure 7(b)
python analysis_stableconsensus.py stubborns_parameters \
    --model crossinh \
    --N 500 \
    --stubborn_type c \
    --plot

#Figure 8(a)
python analysis_stableconsensus.py stubborns_groupsizes \
    --model crossinh \
    --group_sizes 10 30 50 70 100 300 500 700 1000 \
    --stubborn_type z \
    --plot

#Figure 9(a)
python analysis_stableconsensus.py stubborns_groupsizes \
    --model crossinh \
    --group_sizes 10 30 50 70 100 300 500 700 1000 \
    --stubborn_type c \
    --plot