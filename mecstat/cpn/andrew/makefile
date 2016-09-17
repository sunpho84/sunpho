# Compiler is icc
CC=g++
ICC=icc

# Define variables for simulation
# Number of threads
NUM_THREADS?=4
# Number of threads to use for analysis
A_THREADS?=8
# Number of threads per core
SMT=1
# Side length of 2D lattice
LENGTH?=80
# Beta (higher beta, higher correlation length)
BETA?=0.8
# Number of thermalisation steps
THERM_STEPS?=20000
# Number of update steps
UP_STEPS?=100000
# Number of steps to drop in analysis
DROP_STEPS?=0
# Path to home directory
HOME_DIR?=/home/al1g13
# VCL include directory
VCL_INC_DIR?=/usr/include
# dSFMT include directory
dSFMT_INC_DIR?=/usr/include/dSFMT
# NTL include directory
NTL_INC_DIR?=$(HOME_DIR)/Documents/ntl-6.2.1/include
# Restart or not
RESTART?=0
# Save measurements or not
SAVE_MEAS?=1
# Load thermalised state or not
LOAD_THERM=0
# Directory to save executable to
SAVE_DIR?=$(HOME_DIR)/Documents
# Measurements to perform
MEASUREMENTS?=15
# Number of updates between correlator measurements
CORR_N_COR?=1
# Highest value of autocorrelation function to compute for chi_m
MAG_WIDTH?=50000
# Highest value of autocorrelation function to compute for chi_t
TOP_WIDTH?=50000
# Max cutoff to calculate for chi_m autocorrelation array
MS_MAX?=2000
# Max cutoff to calculate for chi_t autocorrelation array
TS_MAX?=20
# If using no boundary conditions, set feducial volume size.
L_GAP?=20

# Set compiler flags
CFLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fvariable-expansion-in-unroller -fprefetch-loop-arrays -frerun-loop-opt -fabi-version=6 -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_THERM=$(THERM_STEPS) -DN_CF=$(UP_STEPS) -DMEASUREMENTS=$(MEASUREMENTS) -DCORR_N_COR=$(CORR_N_COR) -DOMP_NUM_THREADS=$(NUM_THREADS) -DSMT=$(SMT) -DDSFMT_MEXP=1279 -DSAVE_MEASUREMENTS=$(SAVE_MEAS) -DRESTART=$(RESTART) -DLOAD_THERM=$(LOAD_THERM) -I$(dSFMT_INC_DIR) -I$(dSFMT_INC_DIR)/jump -I$(VCL_INC_DIR) -I$(NTL_INC_DIR)
SSE_CFLAGS=-std=c++11 -O3 -fopenmp -msse2 -funroll-loops -fvariable-expansion-in-unroller -fprefetch-loop-arrays -frerun-loop-opt -fabi-version=6 -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_THERM=$(THERM_STEPS) -DN_CF=$(UP_STEPS) -DMEASUREMENTS=$(MEASUREMENTS) -DCORR_N_COR=$(CORR_N_COR) -DOMP_NUM_THREADS=$(NUM_THREADS) -DSMT=$(SMT) -DDSFMT_MEXP=1279 -DSAVE_MEASUREMENTS=$(SAVE_MEAS) -DRESTART=$(RESTART) -DLOAD_THERM=$(LOAD_THERM) -I$(dSFMT_INC_DIR) -I$(dSFMT_INC_DIR)/jump -I$(VCL_INC_DIR) -I$(NTL_INC_DIR)
E_FLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fabi-version=6 -DOMP_NUM_THREADS=$(A_THREADS) -I$(VCL_INC_DIR) -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_CF=$(UP_STEPS) -DDROP=$(DROP_STEPS) -DMEASUREMENT=0 -DERROR=1 -DAUTO=0
MAG_FLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fabi-version=6 -DOMP_NUM_THREADS=$(A_THREADS) -I$(VCL_INC_DIR) -DWIDTH=$(MAG_WIDTH) -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_CF=$(UP_STEPS) -DDROP=$(DROP_STEPS) -DMEASUREMENT=1 -DERROR=1 -DAUTO=1 -DMAX=$(MS_MAX)
G_1_FLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fabi-version=6 -DOMP_NUM_THREADS=$(A_THREADS) -I$(VCL_INC_DIR) -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_CF=$(UP_STEPS) -DDROP=$(DROP_STEPS) -DMEASUREMENT=2 -DERROR=1 -DAUTO=0
TOP_FLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fabi-version=6 -DOMP_NUM_THREADS=$(A_THREADS) -I$(VCL_INC_DIR) -DWIDTH=$(TOP_WIDTH) -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_CF=$(UP_STEPS) -DDROP=$(DROP_STEPS) -DMEASUREMENT=3 -DERROR=1 -DAUTO=1 -DMAX=$(TS_MAX)

MAG_AUTO_FLAGS=-std=c++11 -O3 -fopenmp -march=native -funroll-loops -fabi-version=6 -DOMP_NUM_THREADS=$(A_THREADS) -I$(VCL_INC_DIR) -DWIDTH=$(MAG_WIDTH) -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_CF=$(UP_STEPS) -DDROP=$(DROP_STEPS) -DMEASUREMENT=1 -DJACK_AUTO=1 -DMAX=$(MS_MAX)
#CFLAGS=-std=c++11 -O3 -fopenmp -funroll-loops -fvariable-expansion-in-unroller -fprefetch-loop-arrays -frerun-loop-opt -fabi-version=6 -msse4 -DLEN=$(LENGTH) -DBETA=$(BETA) -DN_THERM=$(THERM_STEPS) -DN_CF=$(UP_STEPS) -DOMP_NUM_THREADS=$(NUM_THREADS) -DSAVE_MEASUREMENTS=$(SAVE_MEAS) -DRESTART=$(RESTART) -I$(VCL_INC_DIR)
	#$(CC) $(CFLAGS) CPN_acceptance_rate.cpp $(dSFMT_INC_DIR)/dSFMT.c $(dSFMT_INC_DIR)/jump/dSFMT-jump.c -lntl -lm -o CPN
	#$(CC) $(CFLAGS) CPN_sitmo.cpp -o CPN
all:
	$(CC) $(CFLAGS) CPN.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_$(LENGTH) -g
sse:
	$(CC) $(CFLAGS) CPN.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_$(LENGTH) -g

heatbath:
	$(CC) $(CFLAGS) -DHEATBATH=1 CPN.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_hb_$(LENGTH) -g

cooling:
	$(CC) $(CFLAGS) CPN_cooling.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_$(LENGTH) -g
rectangle:
	$(CC) $(CFLAGS) -DRECTANGLE=1 -DL_GAP=$(L_GAP) CPN.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_rect_$(LENGTH)_$(L_GAP) -g
no_bc:
	$(CC) $(CFLAGS) -DNO_BC=1 -DRECTANGLE=1 -DL_GAP=$(L_GAP) CPN.cpp $(dSFMT_INC_DIR)/dSFMT.cpp $(dSFMT_INC_DIR)/jump/dSFMT-jump.cpp -lm -o $(SAVE_DIR)/CPN_no_bc_$(LENGTH)_$(L_GAP) -g
	$(CC) $(E_FLAGS) -DNO_BC=1 -DL_GAP=$(L_GAP) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_no_bc_E_$(LENGTH)_$(L_GAP) -g
	$(CC) $(MAG_FLAGS) -DNO_BC=1 -DL_GAP=$(L_GAP) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_no_bc_mag_sus_$(LENGTH)_$(L_GAP) -g
	$(CC) $(G_1_FLAGS) -DNO_BC=1 -DL_GAP=$(L_GAP) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_no_bc_G_1_$(LENGTH)_$(L_GAP) -g
	$(CC) $(TOP_FLAGS) -DNO_BC=1 -DL_GAP=$(L_GAP) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_no_bc_top_sus_$(LENGTH)_$(L_GAP) -g

analysis:
	$(CC) $(E_FLAGS) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_E_$(LENGTH) -g
	$(CC) $(MAG_FLAGS) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_mag_sus_$(LENGTH) -g
	$(CC) $(G_1_FLAGS) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_G_1_$(LENGTH) -g
	$(CC) $(TOP_FLAGS) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_top_sus_$(LENGTH) -g

heatbath_analysis:
	$(CC) $(E_FLAGS) -DHEATBATH=1 CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_hb_analysis_E_$(LENGTH) -g
	$(CC) $(MAG_FLAGS) -DHEATBATH=1 CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_hb_analysis_mag_sus_$(LENGTH) -g
	$(CC) $(G_1_FLAGS) -DHEATBATH=1 CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_hb_analysis_G_1_$(LENGTH) -g
	$(CC) $(TOP_FLAGS) -DHEATBATH=1 CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_hb_analysis_top_sus_$(LENGTH) -g

mag_sus_auto_bin:
	$(CC) $(MAG_AUTO_FLAGS) CPN_perform_full_analysis.cpp -lm -o $(SAVE_DIR)/CPN_analysis_mag_sus_auto_$(LENGTH) -g

