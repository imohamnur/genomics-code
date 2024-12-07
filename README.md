# OUR FINAL PROJECT
## Haroon Shahzad, Taha Rao, Nafi Khan, Danial Syed

# To Run:

## Trimmomatic Command:
### ```java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 covid_files/sars_spike_protein_raw_reads.fastq output.fastq```

- This trims the amount of sequences down from 5000 to 4780, trimming 220 sequences or about 4.4% of the sequences. The quality expected in a 4 character window for the Trimmomatic command was 30%, so the rest 95.6% sequences met the 30% margin.

## Lighter Command
### First, we need to import our Lighter library by cloning the github repo
### ```git clone https://github.com/mourisl/Lighter.git```
### Now we need to run the cloned repo
### ```make```
### next, cd into the Lighter folder:
### ```cd Lighter```
### Then run this command:
### ```./lighter -r ../output.fastq -K 4 4700000 -t 4 -od corrected_output```

- This corrects the sequences after they are trimmed to make sure that they are properly formatted, the average coverage when running this command was 0.161 and the alpha was 43.576
- The k-mer length to analyze over a 4 character window and made 4-mers.
- 3408 reads out of 5000 are error-free, corrected 2169 bases