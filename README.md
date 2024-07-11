Thank you for using PhagePlasmidQuantifier, the purpose of this tool is to locate Phage-Plasmids, a unique mobile genetic element, to do this, an HMMsearch is performed and compared with HMM profiles gathers that are either plasmid or phage assosiated.
When a HMM is hit for a scaffold, it will be marked and counted, when a final percentage of hit HMMs is given by measuring hits over total genes.
The proteins that hit the given HMMs are quantified for each
This output will give an overall idea of how many plasmid assosiated genes are found within a phage, or how many phage assosiated genes are found within a plasmid
Using this information as a start, further investigation into each scaffold can be performed, this tool is a useful step in location of potenial Phage-Plasmids
If there is a scaffold you wish to investigate further, within this repository is a tool called ScaffoldExtract

USAGE (PhagePlasmidQuantifier)
To start using PhagePlasmidQuantifier, make sure you have the following python dependancies (pandas, sys, subprocess, re)
You will also need HMMER installed onto your enviroment, instructions can be found here https://github.com/EddyRivasLab/hmmer 
Here is the command to run PhagePlasmidQuantifier
"python3 PhagePlasmidQuantifier.py 'HMM database' 'Input fasta .faa' 'Genomad or Vibrant output' 'Final Output Name'"
HMM database depends on type of scaffold, use the 'Profiles of Phages' on plasmids and the 'Profile of Plasmids' on phages
Inputs should be in the .faa format, this output can be obtained using both Genomad and VIBRANT
Currently only Genomad output is supported, this should be the file 'GENOMADNAME_summary/GENOMADNAME_plasmid_summary' or 'GENOMADNAME_summary/GENOMADNAME_phage_summary'


FOR PETER'S REFERENCE, ADJUST THE GENOMAD OR VIBRANT TO WORK FOR A CUSTOM FORMAT

USAGE (ScaffoldExtract)
To use ScaffoldExtract
"python3 ScaffoldExtract.py 'MasterFastaFile' 'Scaffold name' 'Scaffoldname.fasta'






