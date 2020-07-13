import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Annotate VCF using snpEff")

parser.add_argument('input_vcf', help='Annotated VCF File')
parser.add_argument('output_database', help='Name of output database (name needs to end in .db)')
parser.add_argument('--fam_file', help="If you have family's in the VCF file, include a fam file")
parser.add_argument('--num_cores', help='Loading will go quicker if more cores are available', default=2)

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
databaseName = args.output_database
famFile = args.fam_file
cores = args.num_cores

# Download necessary requirments
os.system("cd / && git clone https://github.com/quinlan-lab/vcf2db && cd vcf2db \
&& conda install -y -c conda-forge snappy && conda install -y -c conda-forge python-snappy \
&& conda install -y -c bioconda cyvcf2 && conda install -y -c bioconda peddy && pip install -r requirements.txt \
&& cd /proj")

# Load annotated file into a GEMINI database
if famFile != None:
    os.system(f"python /vcf2db/vcf2db.py {inputFile} \
    {famFile} {databaseName} --expand gt_alt_freqs")
elif famFile == None:
    os.system(f"python /vcf2db/vcf2db.py {inputFile} \
    {databaseName} --expand gt_alt_freqs")

# Print output information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')