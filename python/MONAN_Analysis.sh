#!/bin/bash

##########################
# Analysis scripts files #
##########################

prec="MONAN_Prec.py"
bias="MONAN_Bias.py"
mae="MONAN_Mae.py"
skill="MONAN_Skill.py"

#####################################################################
# Input and output directories                                      #
# NETCDF_PATH: Where the comparison NetCDF files are located        #
# MONAN_PATH: Where the MONAN simulation output files are located   #
# OUTPUT_PATH: Where the analysis results will be saved             #
#####################################################################

NETCDF_PATH="/home2/eduardo.eras/workspace/NetCDFs/"
MONAN_PATH="/p/projetos/monan_atm/eduardo.eras/MONAN/scripts_CD-CT/dataout/"
OUTPUT_PATH="/home2/eduardo.eras/workspace/python/output/"

#####################################################################
# Erase old output files to avoid confusion (optional, be careful!) #
#####################################################################
#rm -rf ${OUTPUT_PATH}/*

##################################
# Date range for the simulations #
##################################

#Full range simulation
#start_date="2025-11-26" 
#end_date="2025-12-30"

#Short range testing simulation
start_date="2025-12-29" 
end_date="2025-12-30"

###################################################
# Simulation span in hours (e.g., 120 for 5 days) #
###################################################

LENGTH=120

####################################################
# Precipitation threshold for skill score analysis #
####################################################

THRESHOLD=(1 2 5 10 20 50)

##################################
# Functions, Methods and Imports #
##################################

#Load CDO monule for data manipulation
module load cdo

#Function to print headers for better readability
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

#################################################################################
#                                                                               #
#                              Main analysis loop                               #
#                                                                               #
#################################################################################

print_header "* * *              MONAN Analysis starting...              * * *"

echo -e "\n##############################################"
echo "  _____   _______    ___    _____    _______  "
echo " / ____| |__   __|  /   \  |  __ \  |__   __| " 
echo "| (___      | |    /  ^  \ | |__) |    | |    "
echo " \___ \     | |   /  /_\  \|  _  /     | |    "
echo " ____) |    | |  /  _____  \ | \ \     | |    "
echo "|_____/     |_| /__/     \__\|  \_\    |_|    "
echo -e "##############################################\n"



#Convert dates to Unix timestamps for looping
start_ts=$(date -d "$start_date" +%s)
end_ts=$(date -d "$end_date" +%s)

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
    python ${prec} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${MONAN_PATH} ${OUTPUT_PATH}

    print_status "Bias analysis"
    python ${bias} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${NETCDF_PATH} ${OUTPUT_PATH}

    print_status "MAE analysis"
    python ${mae} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${NETCDF_PATH} ${OUTPUT_PATH}

    print_status "Skill score analysis"
    for THR in "${THRESHOLD[@]}"; do
        print_status "Skill score analysis for threshold: ${THR} mm"
        python ${skill} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${THR} ${NETCDF_PATH} ${OUTPUT_PATH}
    done

done



print_header "* * *               MONAN Analysis finished.               * * *"

echo -e "\n############################################"
echo " ______  _____  _   _  _____  _____  _    _ "
echo "|  ____||_   _|| \ | ||_   _|/ ____|| |  | |"
echo "| |__     | |  |  \| |  | | | (___  | |__| |"
echo "|  __|    | |  | .   |  | |  \___ \ |  __  |"
echo "| |      _| |_ | |\  | _| |_ ____) || |  | |"
echo "|_|     |_____||_| \_||_____||_____/|_|  |_|"
echo -e "############################################\n"
