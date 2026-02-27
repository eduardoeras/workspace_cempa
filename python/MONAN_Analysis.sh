#!/bin/bash

##########################
# Analysis scripts files #
##########################

prec="MONAN_Prec.py"                    #OK
bias="MONAN_Bias.py"                    #OK
mae="MONAN_Mae.py"                      #OK
skill="MONAN_Skill.py"                  #OK
monthly_mean="MONAN_monthly_mean.sh"    #OK
mean_bias="MONAN_mean_Bias.py"          #OK
mean_Mae="MONAN_mean_Mae.py"            #OK
mean_RMSE="MONAN_mean_RMSE.py"          #OK
mean_skill="MONAN_mean_Skill.py"        #OK
heatmap="MONAN_Heatmap.py"

#####################################################################
# Input and output directories                                      #
# NETCDF_PATH: Where the comparison NetCDF files are located        #
# MONAN_PATH: Where the MONAN simulation output files are located   #
# OUTPUT_PATH: Where the analysis results will be saved             #
#####################################################################

NETCDF_PATH="/home2/eduardo.eras/workspace/NetCDFs/"
MONAN_PATH="/p/projetos/monan_atm/eduardo.eras/MONAN/scripts_CD-CT/dataout/"
OUTPUT_PATH="/home2/eduardo.eras/workspace/python/output/"
LISTA_PATH="/home2/eduardo.eras/workspace/python/listas_medias"

#####################################################################
# Erase old output files to avoid confusion (optional, be careful!) #
#####################################################################
rm -rf ${OUTPUT_PATH}/*
rm -rf ${LISTA_PATH}/*

##########################################################
# Date range for the simulations (Minimum range: 5 days) #
##########################################################

#Full range simulation
#start_date="2025-11-26" 
#end_date="2025-12-30"

#Short range testing simulation
start_date="2025-12-25" 
end_date="2025-12-30"

###################################################
# Simulation span in hours (e.g., 120 for 5 days) #
###################################################

LENGTH=120

####################################################
# Precipitation threshold for skill score analysis #
####################################################

THRESHOLD=(1 2 5 10 20 50)

######################################
# Analysis name for output directory #
######################################
ANALYSIS_NAME="5_days_mean" # e.g., "full_range", "short_test", etc.

#####################
# Generate figures? #
#####################
GENERATE_FIGURES=true

#################################################################################
#                                                                               #
#                           .,,uod8B8bou,,.                                     #
#                  ..,uod8BBBBBBBBBBBBBBBBRPFT?l!i:.                            #
#             ,=m8BBBBBBBBBBBBBBBRPFT?!||||||||||||||                           #
#             !...:!TVBBBRPFT||||||||||!!^^""'   ||||                           #
#             !.......:!?|||||!!^^""'            ||||                           #
#             !.........||||                     ||||                           #
#             !.........||||  ##                 ||||                           #
#             !.........||||                     ||||                           #
#             !.........||||                     ||||                           #
#             !.........||||                     ||||                           #
#             !.........||||                     ||||                           #
#             `.........||||                    ,||||                           #
#              .;.......||||               _.-!!|||||                           #
#       .,uodWBBBBb.....||||       _.-!!|||||||||!:'                            #
#    !YBBBBBBBBBBBBBBb..!|||:..-!!|||||||!iof68BBBBBb....                       #
#    !..YBBBBBBBBBBBBBBb!!||||||||!iof68BBBBBBRPFT?!::   `.                     #
#    !....YBBBBBBBBBBBBBBbaaitf68BBBBBBRPFT?!:::::::::     `.                   #
#    !......YBBBBBBBBBBBBBBBBBBBRPFT?!::::::;:!^"`;:::       `.                 #
#    !........YBBBBBBBBBBRPFT?!::::::::::^''...::::::;         iBBbo.           #
#    `..........YBRPFT?!::::::::::::::::::::::::;iof68bo.      WBBBBbo.         #
#      `..........:::::::::::::::::::::::;iof688888888888b.     `YBBBP^'        #
#        `........::::::::::::::::;iof688888888888888888888b.     `             #
#          `......:::::::::;iof688888888888888888888888888888b.                 #
#            `....:::;iof688888888888888888888888888888888899fT!                #
#              `..::!8888888888888888888888888888888899fT|!^"'                  #
#                `' !!988888888888888888888888899fT|!^"'                        #
#                    `!!8888888888888888899fT|!^"'                              #
#                      `!988888888899fT|!^"'                                    #
#                        `!9899fT|!^"'                                          #
#                          `!^"'                                                #
#                                                                               #  
#################################################################################

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

base_analysis() {
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

    done
}

monthly_analysis() {
    processes=("bias" "rmse" "skill")
    INI_VALID=$(date -d "${start_date} 00:00" +%Y%m%d%H)
    FIM_VALID=$(date -d "${end_date} 23:00" +%Y%m%d%H)
    for proc in "${processes[@]}"; do
        print_status "Computing monthly mean for process: ${proc}"
        if [ "${proc}" == "skill" ]; then
            for THR in "${THRESHOLD[@]}"; do
                print_status "Computing monthly mean for skill score with threshold: ${THR} mm"
                bash ${monthly_mean} ${INI_VALID} ${FIM_VALID} ${proc} ${OUTPUT_PATH} ${ANALYSIS_NAME} ${THR}
            done
        else
            bash ${monthly_mean} ${INI_VALID} ${FIM_VALID} ${proc} ${OUTPUT_PATH} ${ANALYSIS_NAME} 
        fi
    done

    print_status "Skill score analysis"
    for THR in "${THRESHOLD[@]}"; do
        print_status "Skill score analysis for threshold: ${THR} mm"
        python ${skill} ${YEAR} ${MONTH} ${DAY} ${HOUR} ${LENGTH} ${THR} ${NETCDF_PATH} ${OUTPUT_PATH}
    done
}

mean_analysis() {
    print_status "Mean bias, MAE, RMSE and skill score across all lead times"
    ano=${start_date:0:4}
    mes=${start_date:5:2}
    dominio=(GLB AMS ACC)
    for dom in "${dominio[@]}"; do
        python ${mean_bias} ${ano} ${mes} ${dom} ${ANALYSIS_NAME} ${NETCDF_PATH} ${OUTPUT_PATH}
        python ${mean_Mae} ${ano} ${mes} ${dom} ${ANALYSIS_NAME} ${NETCDF_PATH} ${OUTPUT_PATH}
        python ${mean_RMSE} ${ano} ${mes} ${dom} ${ANALYSIS_NAME} ${NETCDF_PATH} ${OUTPUT_PATH}
    done
}

skill_analysis() {
    print_status "Mean skill score across all lead times for different thresholds"
    ano=${start_date:0:4}
    mes=${start_date:5:2}
    for THR in "${THRESHOLD[@]}"; do
        python ${mean_skill} ${ano} ${mes} ${THR} ${ANALYSIS_NAME} ${NETCDF_PATH} ${OUTPUT_PATH}
    done
}

heatmap_analysis() {
    print_status "Generating heatmaps for skill scores"
    ano=${start_date:0:4}
    mes=${start_date:5:2}
    for reg in "GLB" "AMS" "ACC"; do
        python ${heatmap} ${ano} ${mes} ${reg} ${ANALYSIS_NAME} ${OUTPUT_PATH}
    done
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

base_analysis

monthly_analysis

mean_analysis

skill_analysis

heatmap_analysis

print_header "* * *               MONAN Analysis finished.               * * *"

echo -e "\n############################################"
echo " ______  _____  _   _  _____  _____  _    _ "
echo "|  ____||_   _|| \ | ||_   _|/ ____|| |  | |"
echo "| |__     | |  |  \| |  | | | (___  | |__| |"
echo "|  __|    | |  | .   |  | |  \___ \ |  __  |"
echo "| |      _| |_ | |\  | _| |_ ____) || |  | |"
echo "|_|     |_____||_| \_||_____||_____/|_|  |_|"
echo -e "############################################\n"
