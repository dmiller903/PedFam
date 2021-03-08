import gzip
import time
import os
import concurrent.futures
import argparse
import re

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Creates a file that has clinVar, gnomAD AF, 1K AF, and CADD\
    values for all variant positions.")

parser.add_argument('output_file', help='Name of output file')

args = parser.parse_args()

#Create variables of each argument from argparse
outputFile = args.output_file

def create_cadd_dict(currentChr):
    tempDict = {}
    with gzip.open(f"/tmp/cadd_chr{currentChr}.tsv.gz", "rt") as caddFile:
        for line in caddFile:
            if line.startswith("#Chrom"):
                lineList = line.rstrip("\n").split("\t")
                posIndex = lineList.index("Pos")
                refIndex = lineList.index("Ref")
                altIndex = lineList.index("Alt")
                phredIndex = lineList.index("PHRED")
            else:
                lineList = line.rstrip("\n").split("\t")
                pos = lineList[posIndex]
                ref = lineList[refIndex]
                alt = lineList[altIndex]
                phred = lineList[phredIndex]
                if pos not in tempDict:
                    tempDict[pos] = [[ref, alt, phred]]
                elif pos in tempDict:
                    tempDict[pos].append([ref, alt, phred])
        return(tempDict)

def create_clinVar_dict(currentChr):
    tempDict = {}
    with gzip.open(f"/tmp/clinVar_chr{currentChr}.tsv.gz", "rt") as clinVarFile:
        for line in clinVarFile:
            if line.startswith("rsID"):
                lineList = line.rstrip("\n").split("\t")
                rsIndex = lineList.index("rsID")
                sigIndex = lineList.index("clinicalSig")
            else:
                lineList = line.rstrip("\n").split("\t")
                rs = lineList[rsIndex]
                sig = lineList[sigIndex]
                if rs not in tempDict:
                    tempDict[rs] = {sig}
                else:
                    tempDict[rs].add(sig)
        return(tempDict)

def create_maf_dict(currentChr):
    tempDict = {}
    fileName = f"/tmp/1K_chr{currentChr}.tsv.gz"
    with gzip.open(fileName, "rt") as mafFile:
        for line in mafFile:
            if line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                posIndex = lineList.index("POS")
                refIndex = lineList.index("REF")
                altIndex = lineList.index("ALT")
                afIndex = lineList.index("AF")
            else:
                lineList = line.rstrip("\n").split("\t")
                pos = lineList[posIndex]
                ref = lineList[refIndex]
                alt = lineList[altIndex]
                af = lineList[afIndex] 
                if pos not in tempDict:
                    tempDict[pos] = [[ref, alt, af]]
                elif pos in tempDict:
                    tempDict[pos].append([ref, alt, af])  
        return(tempDict)

def get_value_from_dict(aDict, pos, ref, alt):
    if pos in aDict:
        for posList in aDict[pos]:
            if posList[0] == ref and posList[1] == alt:
                aValue = posList[-1]
                return(aValue)
    return("None")

def parse_cadd(file):
    currentChrom = re.findall(r"_chr(\d+)\.", file)[0]
    with gzip.open("/tmp/whole_genome_SNVs.tsv.gz", 'rt') as cadd, gzip.open(file, "wb") as output:
        meta = cadd.readline()
        header = cadd.readline()
        headerList = header.rstrip("\n").split("\t")
        chromIndex = headerList.index("#Chrom")
        output.write(header.encode())
        for line in cadd:
            lineList = line.rstrip("\n").split("\t")
            chrom = lineList[chromIndex]
            if chrom == currentChrom:
                output.write(line.encode())

def parse_clinVar(file):
    currentChrom = re.findall(r"_chr(\d+)\.", file)[0]
    with gzip.open("/tmp/variant_summary.txt.gz", 'rt') as clinVar, gzip.open(file, "wb") as output:
        header = clinVar.readline()
        headerList = header.rstrip("\n").split("\t")
        chromIndex = headerList.index("Chromosome")
        rsIndex = headerList.index("RS# (dbSNP)")
        sigIndex = headerList.index("ClinicalSignificance")
        output.write(f"rsID\tclinicalSig\n".encode())
        for line in clinVar:
            lineList = line.rstrip("\n").split("\t")
            chrom = lineList[chromIndex]
            if chrom == currentChrom:
                if lineList[rsIndex] != "-":
                    rs = "rs" + lineList[rsIndex]
                    sig = lineList[sigIndex]
                    output.write(f"{rs}\t{sig}\n".encode())

def getFiles(file):
        os.system(f"wget --no-check-certificate {file}")

def parse_1K(file):
    currentChrom = re.findall(r"_chr(\d+)\.", file)[0]
    with gzip.open(f"/tmp/ALL.chr{currentChrom}_GRCh38_sites.20170504.vcf.gz", 'rt') as mafFile, gzip.open(file, "wb") as output:
        for line in mafFile:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                chromIndex = lineList.index("#CHROM")
                posIndex = lineList.index("POS")
                refIndex = lineList.index("REF")
                altIndex = lineList.index("ALT")
                infoIndex = lineList.index("INFO")
                rsIndex = lineList.index("ID")
                output.write(f"#CHROM\tPOS\tREF\tALT\tAF\n".encode())
            else:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex]
                if currentChrom == chrom:
                    pos = lineList[posIndex]
                    ref = lineList[refIndex]
                    alt = lineList[altIndex]
                    info = lineList[infoIndex]
                    rs = lineList[rsIndex]
                    rs = rs.split(";")
                    alt = alt.split(",")
                    infoList = info.split(";")
                    maf = infoList[1].replace("AF=", "")
                    maf = maf.split(",")
                    for j, value in enumerate(alt):
                        output.write(f"{chrom}\t{pos}\t{ref}\t{value}\t{maf[j]}\n".encode())

def parse_gnomAD(file):
    currentChrom = re.findall(r"_chr(\d+)_", file)[0]
    with gzip.open("/tmp/gnomad.genomes.r3.0.sites.vcf.bgz", 'rt') as gnomAD, gzip.open(file, "wb") as output:
        output.write(f"#CHROM\tPOS\tREF\tALT\trsID\tAF\n".encode())
        for line in gnomAD:
            if not line.startswith("#"):
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[0].lstrip("chr")
                if chrom == currentChrom:
                    pos = lineList[1]
                    rsID = lineList[2]
                    ref = lineList[3]
                    alt = lineList[4]
                    # info is the last item in the list and the elements of info are separated by :
                    info = lineList[-1]
                    infoList = info.split(";")
                    af = infoList[2]
                    #remove "AF=" from af and convert to float instead of scientific notation.
                    if "AF=" in af:
                        af = float(af[3:])
                        af = f"{af:.10f}"
                        output.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rsID}\t{af}\n".encode())

# Download the CADD file and separate it by chromosome
if not os.path.exists("/tmp/whole_genome_SNVs.tsv.gz"):
    os.system("wget --no-check-certificate https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz -P /tmp")

fileList = []
for i in range(1, 23):
    fileList.append(f"/tmp/cadd_chr{i}.tsv.gz")

with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(parse_cadd, fileList)

# Download the clinVar file and separate by chromosome
if not os.path.exists("/tmp/variant_summary.txt.gz"):
    os.system("wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz -P /tmp")

fileList = []
for i in range(1, 23):
    fileList.append(f"/tmp/clinVar_chr{i}.tsv.gz")

with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(parse_clinVar, fileList)
    
# Download 1000 genome maf file and separate by chromosome
if not os.path.exists(f"/tmp/ALL.chr1_GRCh38_sites.20170504.vcf.gz"):
    fileList = []
    for i in range(1, 23):
        fileList.append(f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{i}_GRCh38_sites.20170504.vcf.gz -P /tmp")

    with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
        executor.map(getFiles, fileList)

fileList = []
for i in range(1, 23):
    fileList.append(f"/tmp/1K_chr{i}.tsv.gz")
with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(parse_1K, fileList)

# Download gnomAD maf file and separate by chromosome
if not os.path.exists("/tmp/gnomad.genomes.r3.0.sites.vcf.bgz"):
    os.system("wget --no-check-certificate https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz -P /tmp")

fileList = []
for i in range(1, 23):
    fileList.append(f"/tmp/gnomAD_chr{i}_AF.tsv.gz")
                        
with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(parse_gnomAD, fileList)

# Create a dictionaries containing clinvar, maf and cadd information and out to file
for i in range(1, 23):
    #Current integer is current chromosome
    currentChr = str(i)

    #Create empty dictionaries
    clinVarDict = {}
    mafDict = {}
    caddDict = {}

    #Create a dictionary with ClinVar information
    clinVarDict = create_clinVar_dict(currentChr)
    timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
    timeElapsedHours = round(timeElapsedMinutes / 60, 2)
    print(f'{char}ClinVar dictionary created for chromosome {currentChr}. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')

    #Create a dictionary with 1K information
    mafDict = create_maf_dict(currentChr)
    timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
    timeElapsedHours = round(timeElapsedMinutes / 60, 2)
    print(f'{char}1K dictionary created for chromosome {currentChr}. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')

    #Create a dictionary with CADD information
    caddDict = create_cadd_dict(currentChr)
    timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
    timeElapsedHours = round(timeElapsedMinutes / 60, 2)
    print(f'{char}CADD dictionary created for chromosome {currentChr}. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')

    #Iterate through gnomAD file, outputting a new file with gnomAD, CADD, 1K, and Clinvar information
    with gzip.open(f"/tmp/gnomAD_chr{currentChr}_AF.tsv.gz", "rt") as gnomAD:
        if currentChr == "1":
            with gzip.open(outputFile, "wb") as output:
                for line in gnomAD:
                    if line.startswith("#CHROM"):
                        lineList = line.rstrip("\n").split("\t")
                        chromIndex = lineList.index("#CHROM")
                        posIndex = lineList.index("POS")
                        refIndex = lineList.index("REF")
                        altIndex = lineList.index("ALT")
                        afIndex = lineList.index("AF")
                        rsIndex = lineList.index("rsID")
                        output.write(f"#CHROM\tPOS\trsID\tREF\tALT\tcadd\tclinVar\t1K_AF\tgnomAD_AF\n".encode())
                    else:
                        lineList = line.rstrip("\n").split("\t")
                        chrom = lineList[chromIndex]
                        if chrom == currentChr:
                            pos = lineList[posIndex]
                            ref = lineList[refIndex]
                            alt = lineList[altIndex]
                            af = lineList[afIndex]
                            rs = lineList[rsIndex]
                            cadd = get_value_from_dict(caddDict, pos, ref, alt)
                            maf = get_value_from_dict(mafDict, pos, ref, alt)
                            if rs != "." and rs in clinVarDict:
                                sigSet = clinVarDict[rs]
                                sig = ", ".join(sigSet)
                            else:
                                sig = "None" 
                            output.write(f"{currentChr}\t{pos}\t{rs}\t{ref}\t{alt}\t{cadd}\t{sig}\t{maf}\t{af}\n".encode())
        else:
            with gzip.open(outputFile, "ab") as output:
                for line in gnomAD:
                    if line.startswith("#CHROM"):
                        lineList = line.rstrip("\n").split("\t")
                        chromIndex = lineList.index("#CHROM")
                        posIndex = lineList.index("POS")
                        refIndex = lineList.index("REF")
                        altIndex = lineList.index("ALT")
                        afIndex = lineList.index("AF")
                        rsIndex = lineList.index("rsID")
                    else:
                        lineList = line.rstrip("\n").split("\t")
                        chrom = lineList[chromIndex]
                        if chrom == currentChr:
                            pos = lineList[posIndex]
                            ref = lineList[refIndex]
                            alt = lineList[altIndex]
                            af = lineList[afIndex]
                            rs = lineList[rsIndex]
                            cadd = get_value_from_dict(caddDict, pos, ref, alt)
                            maf = get_value_from_dict(mafDict, pos, ref, alt)
                            if rs != "." and rs in clinVarDict:
                                sigSet = clinVarDict[rs]
                                sig = ", ".join(sigSet)
                            else:
                                sig = "None" 
                            output.write(f"{currentChr}\t{pos}\t{rs}\t{ref}\t{alt}\t{cadd}\t{sig}\t{maf}\t{af}\n".encode())
    
    timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
    timeElapsedHours = round(timeElapsedMinutes / 60, 2)
    print(f'{char}Chromosome {currentChr} written to file. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')