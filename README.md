# The following document describes the data and scripts in their order of command as used to analyze the entries in our dataset. They will be explained in order and should be excecuted in the same way.

**N2 Lifespans FINAL.xlsx**
  - This is the finalized dataset containing all entries and their respective details.

**N2_functions.py**
  - This file should be located in the same folder as all other scripts and data files. 
  - This script contains all functions used to create groups of entries, to sort entries by different experimental conditions (Transfers to fresh plates, FUdR, and temperature), and sort entries in different ways for each figure. 

**temperature.py**
	- This script returns statistical information about the dataset, such as percent of information missing for each experimental condition, number of entries, and unique labs and countries. 
	- This script sorts the entries by growth media and temperature. It then uses functions from 'N2_functions.py' to sort and plot the data by % alive on each day and mean lifespan. 
  
**FUdR.py**
  - This script sorts the entries by growth media, temperature, and use of FUdR. It then uses functions from 'N2_functions.py' to sort and plot the data by % alive on each day and mean lifespan. 
  - This script also gives information about the distribution of FUdR concentrations and analyzes the effect on lifespans as the concentrations vary. 

**manipulations.py**
  - This script sorts the entries by growth media, temperature, and manipulations of plates. It then uses functions from 'N2_functions.py' to sort and plot the data by % alive on each day and mean lifespan. 
  - This script tests the effects of different frequencies of plate manipulations on lifespan.
  - This script will also consider the effects of considering both plate manipulations and use of FUdR.

**country_state.py**
  - This script takes the same dataset formed in 'temperature.py' and analyzes the effects on % alive for each day and mean lifespan by country and state of the lab location from which the entry was taken.

**all_mean_lifespans.py**
  - This script retrieves the mean lifespan data from all assortments of data found in 'temperature.py', 'FUdR.py', and 'manipulations.py'. 
  - This script takes all of the mean lifespans of the entire dataset and order them from smallest to largest. It will also group and save mean lifespan data from all groups in the dataset created in 'temperature.py'.

** If any data would like to be saved from these scripts then a folder named 'saved_data' and 'figures' will need to be created in the same folder all of the above mentioned scripts are located. **
