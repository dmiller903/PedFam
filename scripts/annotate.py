import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Annotate VCF using snpEff. GRCH37.75 is used as the reference genome.")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_file', help='Name of output VCF file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file

# Download annotation files if necessary
if not os.path.exists("/snpEff/./data/GRCh38.86/sequence.1.bin"):
    os.system("java -jar /snpEff/snpEff.jar download -v GRCh38.86")

# Annotate the vt trimmed file
os.system(f"java -Xmx40g -jar /snpEff/snpEff.jar GRCh38.86 -v \
{inputFile} > \
{outputFile}")

# Print output information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')