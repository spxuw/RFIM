# RFIM

This is the repository for the "Disease module detection using random-field Ising model (RFIM)".

# Installation
* Specify the path:  (1) open the Makefile in folder "BUILD". (2) Revise the directory of BASEBGL, BASENR, BASE and LIBRARY.
* Build the code: (1): change to working directory to the "Source/BUILD/". (2): type: "make" to compile code.

# code
* Source: C++ codes for the disease module detection using RFIM.
* Prepare_data_RFIM.R: prepare the input data format of RFIM.
* figure_code: code to generate figures.

# data
All relevant data have been processed and are saved in the folder "Source/BUILD/data".

# Run RFIM
* Run RFIM: RFIM accepts several command formats to calculate the M-H curve. See Source/MAIN/main_Network_Real.cpp for details.
* Example command: ./HMS_Network_Real 22 ./data/ Asthma_new_ppirevise_string_v3_f_0.txt 100 -8 8. This will output the module (less than 5% of the total nodes) through calculating the M-H curve ranging from -8 to 8 with 100 H values. The code will stop as long as the module size exceeds 5% of the total nodes.
