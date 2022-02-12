#!/bin/bash

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.txt")

declare -a H_min=(-6 -5 -5)
declare -a H_max=(-4 -3 -3)

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.1.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.1.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.1.txt")

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.2.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.2.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.2.txt")

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.3.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.3.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.3.txt")

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.4.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.4.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.4.txt")

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_string_v3_f_0.5.txt" "Breast_cancer_new_ppirevise_string_v3_f_0.5.txt" "Lung_cancer_new_ppirevise_string_v3_f_0.5.txt")

for ((i=0;i<0;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min[$i]} ${H_max[$i]}
done


declare -a H_min_2=(-7 -6 -6)
declare -a H_max_2=(-5 -4 -4)

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.txt")


for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.1.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.1.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.1.txt")

for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.2.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.2.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.2.txt")

for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.3.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.3.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.3.txt")

for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.4.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.4.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.4.txt")

for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done

declare -a phenotype=("Asthma_new_ppirevise_iRefIndex_v3_f_0.5.txt" "Breast_cancer_new_ppirevise_iRefIndex_v3_f_0.5.txt" "Lung_cancer_new_ppirevise_iRefIndex_v3_f_0.5.txt")

for ((i=0;i<3;i++)); 
do
	./HMS_Network_Real 22 ./data/Dandi/ "${phenotype[$i]}" 100 ${H_min_2[$i]} ${H_max_2[$i]}
done
