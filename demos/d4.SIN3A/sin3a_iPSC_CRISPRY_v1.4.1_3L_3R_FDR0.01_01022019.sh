#!/bin/bash


###### notes on setting three key parameters:
## -n is fdr for sgRNA cut off. Yarui Used 0.05 in the nature method paper. Increase this value up to 1 to get more peaks.
## -a is rra peaks. Crispy default is use 1 to keep all (recommended). 
## -c higher is more stringent. 3 means 0.001 peaks filter in macs2. Lower this value down to 0+ to get more peaks. 


##
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_XR032_XR037(G-M+)" \
            -b "XR029,XR034" \
            -f "XR032,XR037" \
           	-n 0.029

           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_XR031_XR036(G+M-)" \
            -b "XR029,XR034" \
            -f "XR031,XR036" \
           	-n 0.028
           	
           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_cis-" \
            -b "XR029,XR034" \
            -f "XR031,XR036,XR032,XR037" \
           	-n 0.029
           	          	
         	

## quantile normalization within negative population   

./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_XR032_XR037(G-M+)_QWN" \
            -b "XR029,XR034" \
            -f "XR032,XR037" \
           	-q "XR029,XR034;XR030,XR035;XR032,XR037,XR031,XR036" \
           	-n 0.037
           	
           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_XR031_XR036(G+M-)_QWN" \
            -b "XR029,XR034" \
            -f "XR031,XR036" \
           	-q "XR029,XR034;XR030,XR035;XR032,XR037,XR031,XR036" \
           	-n 0.025
           	
           	
./crispy.sh -i $HOME/crispy/SIN3A_input/sin3a_ipsc_reads_3L_3R_12262018.tsv \
            -r $HOME/crispy/SIN3A_input/sin3a_regions.bed \
            -s $HOME/crispy/SIN3A_input/sin3a_oligos.tsv \
            -o "results.sin3a_iPSC_CRISPY_1.4.1_3L_3R_12262018_FDR0.01_01022019" \
            -p "sin3a_iPSC_Bin_XR029_XR034(unsorted)_cis-QWN" \
            -b "XR029,XR034" \
            -f "XR032,XR037,XR031,XR036" \
           	-q "XR029,XR034;XR030,XR035;XR032,XR037,XR031,XR036" \
           	-n 0.025
           	   	
           	