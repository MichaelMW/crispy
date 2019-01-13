#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


## 
           	           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_XR046_XR054(Beaker)(G+M-)" \
            -b "XR043,XR050" \
            -f "XR046,XR054" \
           	-n 0.113 \
           	-c 2

./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_XR045_XR052(G-M+)" \
            -b "XR043,XR050" \
            -f "XR045,XR052" \
           	-n 0.073 \
           	-c 2

./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_cis-" \
            -b "XR043,XR050" \
            -f "XR045,XR052,XR046,XR054" \
           	-n 0.032 \
           	-c 2
           	           	
           	   
## quantile normalization with negative population   
         	           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_XR046_XR054(Beaker)(G+M-)_QWN" \
            -b "XR043,XR050" \
            -f "XR046,XR054" \
            -q "XR043,XR050;XR044,XR051;XR045,XR052,XR046,XR054" \
           	-n 0.068 \
           	-c 2

./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_XR045_XR052(G-M+)_QWN" \
            -b "XR043,XR050" \
            -f "XR045,XR052" \
            -q "XR043,XR050;XR044,XR051;XR045,XR052,XR046,XR054" \
           	-n 0.069 \
           	-c 2
           	
           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_neuron_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_Neuron_CRISPY_v1.4.1_3L_3R_12262018_FDR0.01" \
            -p "sin3a_Neuron_Bin_XR043_XR050(unsorted_2w)_cis-_QWN" \
            -b "XR043,XR050" \
            -f "XR045,XR052,XR046,XR054" \
            -q "XR043,XR050;XR044,XR051;XR045,XR052,XR046,XR054" \
           	-n 0.062 \
           	-c 2
