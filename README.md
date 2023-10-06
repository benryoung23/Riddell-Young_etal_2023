# Riddell-Young_etal_2023
 The Matlab package used in the manuscript entitled "Tropical sources mainly controleld last glacial maximum and deglacial methane variability". The code calculates the methane interpolar difference and methane source distribution across the last glacial maximum and deglaciation. The code used to produce figures 1 to 3 are also included.


The code provides the scripts and data to calculate the methane interpolar difference across the specified climate interval, then use the output to calculate methane source distributions using a four-box troposphere model, then plot the results. The workflow is as follows: 

Step 1: Open and run GISP2_and_NEEMrIPDa.m.
	Input: WAIS_Chronology_Buizert2015.csv, WAISfinal_NEEM_may22.csv, NEEMfinal_nov22.csv, WAISfinal_GISP2_may22.csv, and GISP2final_gradient_nov22.csv.
	Output (saved to project folder): NEEM_Chronology.txt, GISP2_Chronology.txt, gradient.mat, NEEM_IA.mat, and GISP2_IA.mat

Step 2: Open and run NEEM_IPD_Cont_Source_4box_Model.m
	Input: WAIS_Chronology_Buizert2015.csv, NEEMfinal_nov22.csv, WAISfinal_NEEM_may22,
	Output (saved to project folder): NEEM4box.mat

Step 3: Open and run GISP2_IPD_Cont_Source_4box_Model.m
	Input: WAIS_Chronology_Buizert2015.csv, WAISfinal_GISP2_may22.csv, and GISP2final_gradient_nov22.csv.
	Output (saved to project folder): GISP24box.mat

Step 4: Open and run finaldataplot.m, fourbox_fig.m, gradient_fig.m
	- These scripts plot figures one two and three of the main text of Riddell-Young et al., 2023. 
	- The latter two take the output of steps 1 - 3 and other datasets included in the Matlab package to make the figures. These two scripts also use the package jblill.m to plot error bars around the calculated rIPD and emissions.