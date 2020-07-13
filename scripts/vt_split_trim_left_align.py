# Import necessary modules
import os
import time
import argparse
import concurrent.futures
import glob

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Trim and normalize VCF file using vt tools.")

parser.add_argument('input_vcf', help='Name of input file')
parser.add_argument('output_file', help='Name of output file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file.rstrip(".gz")

# Download reference files if needed
if not os.path.exists("/references/Homo_sapiens_assembly38.fasta"):
    os.system("wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9ddd0ba7bc73000ce87e38/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references \
    && rm /tmp/references.zip \
    && gzip -d /tmp/references/Homo_sapiens_assembly38.fasta.gz \
    && gzip -d /tmp/references/Homo_sapiens_assembly38.dict.gz \
    && gzip -d /tmp/references/Homo_sapiens_assembly38.fasta.fai.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}") and fileName != "readme" and "Homo_sapiens_assembly38" in fileName:
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# Use VT to split, trim and left align the phased samples.
os.system(f"/root/miniconda2/bin/vt decompose -s {inputFile} \
| /root/miniconda2/bin/vt normalize -n -r /references/Homo_sapiens_assembly38.fasta - > \
{outputFile} && gzip {outputFile}")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')