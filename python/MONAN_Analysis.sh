#!/bin/bash

#Functions
print_header() {
    echo -e "\n################################################################"
    echo "$1"
    echo -e "################################################################\n"
}

print_status() {
    echo -e "\n----------------------------------------------------------------"
    echo "$1"
    echo -e "----------------------------------------------------------------\n"
}

THRESHOLD=$1

prec="MONAN_Prec.py"
bias="MONAN_Bias.py"
mae="MONAN_Mae.py"
skill="MONAN_Skill.py"



#Load CDO monule for data manipulation
module load cdo

#Date range for the simulations
start_date="2025-11-26" 
#end_date="2025-12-30"  #Full range simulation
end_date="2025-11-27"   #Short range testing simulation

#Simulation span in hours (e.g., 120 for 5 days)
LENGTH=120

#Convert dates to Unix timestamps for looping
start_ts=$(date -d "$start_date" +%s)
end_ts=$(date -d "$end_date" +%s)

print_header "* * *              MONAN Analysis starting...              * * *"

#Looping over dates (the sane way)
#Strategy: Convert dates to Unix timestamps (seconds since 1970-01-01),
#loop numerically, them convert back to human dates.
#This avoids calendar madness.
for (( ts = start_ts; ts <= end_ts; ts += 86400 )); do

    #Get individual date elements
    YEAR=$(date -d "@$ts" +%Y)
    MONTH=$(date -d "@$ts" +%m)
    DAY=$(date -d "@$ts" +%d)
    HOUR=$(date -d "@$ts" +%H)

    simulation_date=$(date -d "@$ts" +%Y%m%d00)
    message="* * * Simulation for date: ${YEAR} ${MONTH} ${DAY} with span: $LENGTH hours * * *"
    print_header "$message"

    #Run the analysis scripts
    print_status "24 hours precipitation accumulation"
    python ${prec} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH}

    print_status "Bias analysis"
    python ${bias} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH}

    print_status "MAE analysis"
    python ${mae} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH}

    print_status "Skill score analysis"
    python ${skill} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${THRESHOLD}

done

print_header "* * *               MONAN Analysis finished.               * * *"