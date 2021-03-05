import os
import time
import argparse
import concurrent.futures
import gzip

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Uses a GEMINI database as input to identify CH variants. If a \
--fam_file' is not used, false positives may result since parental haplotypes not being takin in to consideration.")

parser.add_argument('input_file', help='GEMINI database')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('cadd_maf_file', help='Name of file with CADD and MAF values')
parser.add_argument('--cadd', help='Indicate cadd cut-off value.', default='15')
parser.add_argument('--af', help='Indicate allele frequency cut-off value.', default='0.01')
parser.add_argument('--fam_file', help='If family relationships are known among the samples, use a fam file to help \
with the CH identification process.')
parser.add_argument('--impact_filter_only', help='"y" will only filter based on "HIGH" impact without regard to allele \
frequency or cadd scores', default="n")

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_file
outputFile = args.output_file
caddMafFile = args.cadd_maf_file
inputCadd = float(args.cadd)
inputAF = args.af
if inputAF != "None":
    inputAF = float(args.af)
familyFile = args.fam_file
impactFilterOnly = args.impact_filter_only

#Function to get convert sample genotype from alpha to numeric
def getNumericGenotype(genotype, ref, alt):
    if "." not in genotype:
        if "|" in genotype:
            genotypeList = genotype.split("|")
            genotypeSymbol = "|"
        elif "/" in genotype:
            genotypeList = genotype.split("/")
            genotypeSymbol = "/"
        firstAllele = ""
        secondAllele = ""
        if "," not in alt:
            if genotypeList[0] == ref:
                firstAllele = "0"
            elif genotypeList[0] == alt:
                firstAllele = "1"
            else:
                firstAllele = "."
            if genotypeList[1] == ref:
                secondAllele = "0"
            elif genotypeList[1] == alt:
                secondAllele = "1"
            else:
                secondAllele = "."
            newGenotype = f"{firstAllele}{genotypeSymbol}{secondAllele}"
            return(newGenotype)
        else:
            altList = alt.split(",")
            for i, alt in enumerate(altList):
                if genotypeList[0] == ref:
                    firstAllele = "0"
                elif genotypeList[0] == alt:
                    firstAllele = str(i + 1)
                else:
                    firstAllele = "."
                if genotypeList[1] == ref:
                    secondAllele = "0"
                elif genotypeList[1] == alt:
                    secondAllele = str(i + 1)
                else:
                    secondAllele = "."
                newGenotype = f"{firstAllele}{genotypeSymbol}{secondAllele}"
                return(newGenotype)
    else:
        return(".|.")

#Function to grab information from header of input file
def getHeaderInfo(headerList):
    startIndex = headerList.index("start")
    geneIndex = headerList.index("gene")
    refIndex = headerList.index("ref")
    altIndex = headerList.index("alt")
    impactIndex = headerList.index("impact_severity")
    caddIndex = headerList.index("cadd")
    af_1k_index = headerList.index("af_1k")
    af_gnomAD_index = headerList.index("af_gnomAD")
    lofIndex = headerList.index("is_lof")
    exonicIndex = headerList.index("is_exonic")
    samples = headerList[15:]
    return(startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, af_1k_index, af_gnomAD_index, lofIndex, exonicIndex, samples)

#Function to grab information from line of input file
def getLineInfo(lineList):
    start = lineList[startIndex]
    gene = lineList[geneIndex]
    ref = lineList[refIndex]
    alt = lineList[altIndex]
    impact = lineList[impactIndex]
    cadd = lineList[caddIndex]
    af_1k = lineList[af_1k_index]
    af_gnomAD = lineList[af_gnomAD_index]
    lof = lineList[lofIndex]
    exonic = lineList[exonicIndex]
    return(start, gene, ref, alt, impact, cadd, af_1k, af_gnomAD, lof, exonic)

def iterateThroughSamples():
    for sampleIndex in sampleIndexes:
        sample = headerList[sampleIndex]
        genotype = lineList[sampleIndex]
        newGenotype = getNumericGenotype(genotype, ref, alt)
        if gene not in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene] = [newGenotype]
            samplePositions[sample][gene] = [start]
        elif gene in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene].append(newGenotype)
            samplePositions[sample][gene].append(start)

# Create a .tsv that has all pertinent information for compound heterozygous identification
impactSeverity = "'LOW'"
tempTsv = f"/tmp/gemini.tsv"
geminiTsv = f"{inputFile.replace('.db', '_gemini.tsv')}"
if not os.path.exists(geminiTsv):
    os.system(f'gemini query --header -q "select chrom, start, vcf_id, ref, alt, gene, is_exonic, impact_severity, \
        is_lof, impact, (gts).(*) from variants where impact_severity != {impactSeverity}" \
        {inputFile} \
        > {tempTsv}')
    
    # Create a dictionary of all the chromosome positions in the gemini tsv file
    posDict = {}
    with open(tempTsv) as inputFile:
        for line in inputFile:
            if line.startswith("chrom"):
                lineList = line.rstrip("\n").split("\t")
                chromIndex = lineList.index("chrom")
                posIndex = lineList.index("start")
            else:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex].replace("chr", "")
                pos = str(int(lineList[posIndex]) + 1)
                if chrom not in posDict:
                    posDict[chrom] = {pos}
                elif chrom in posDict:
                    posDict[chrom].add(pos)
    # Use the posDict to determine the cadd, af_1k, and gnomAD values for each chromosomal position
    caddMafDict = {}
    with gzip.open(caddMafFile, 'rt') as inputFile:
        for line in inputFile:
            if line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                chromIndex = lineList.index("#CHROM")
                posIndex = lineList.index("POS")
                refIndex = lineList.index("REF")
                altIndex = lineList.index("ALT")
                af_1k_index = lineList.index("1K_AF")
                af_gnomAD_index = lineList.index("gnomAD_AF")
                caddIndex = lineList.index("cadd")
                rsIndex = lineList.index("rsID")
                clinVarIndex = lineList.index("clinVar")
            else:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex]
                pos = lineList[posIndex]
                ref = lineList[refIndex]
                alt = lineList[altIndex]
                af_1k = lineList[af_1k_index]
                af_gnomAD = lineList[af_gnomAD_index]
                cadd = lineList[caddIndex]
                rs = lineList[rsIndex]
                clinVar = lineList[clinVarIndex]
                if chrom in posDict and pos in posDict[chrom]:
                    if chrom not in caddMafDict:
                        caddMafDict[chrom] = {pos: [[ref, alt, af_1k, af_gnomAD, cadd, rs, clinVar]]}
                    elif chrom in caddMafDict and pos not in caddMafDict[chrom]:
                        caddMafDict[chrom][pos] = [[ref, alt, af_1k, af_gnomAD, cadd, rs, clinVar]]
                    elif chrom in caddMafDict and pos in caddMafDict[chrom]:
                        caddMafDict[chrom][pos].append([ref, alt, af_1k, af_gnomAD, cadd, rs, clinVar])

    #Create a new gemini file that includes the af_1k, af_gnomAD and cadd values present in the caddMafDict
    with open(tempTsv) as inputFile, open(geminiTsv, "w") as geminiFile:
        for line in inputFile:
            if line.startswith("chrom"):
                lineList = line.rstrip("\n").split("\t")
                chromIndex = lineList.index("chrom")
                posIndex = lineList.index("start")
                refIndex = lineList.index("ref")
                altIndex = lineList.index("alt")
                newLine = lineList[0:5] + ["af_1k", "af_gnomAD", "cadd", "rs", "clinvar"] + lineList[5:]
                newLine = "\t".join(newLine)
                geminiFile.write(f"{newLine}\n")
            else:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex].replace("chr", "")
                pos = str(int(lineList[posIndex]) + 1)
                ref = lineList[refIndex]
                alt = lineList[altIndex]
                if chrom in caddMafDict and pos in caddMafDict[chrom]:
                    for posList in caddMafDict[chrom][pos]:
                        matchFound = False
                        if posList[0] == ref and posList[1] == alt:
                            clinVar = posList[-1]
                            rs = posList[-2]
                            cadd = posList[-3]
                            af_gnomAD = posList[-4]
                            af_1k = posList[-5]
                            matchFound = True
                            newLine = lineList[0:1] + [pos] + lineList[2:5] + [af_1k, af_gnomAD, cadd, rs, clinVar] + lineList[5:]
                            newLine = "\t".join(newLine)
                            geminiFile.write(f"{newLine}\n")
                            break
                    if matchFound == False:
                        newLine = lineList[0:1] + [pos] + lineList[2:5] + ["None", "None", "None", "None", "None"] + lineList[5:]
                        newLine = "\t".join(newLine)
                        geminiFile.write(f"{newLine}\n")

# Use fam file to create a list of samples, list of parents, and a parent dictionary where each key is a parent ID and value is sample ID
if familyFile is not None:
    parentDict = {}
    parentList = []
    patientList = []
    familyDict = {}
    with open(familyFile) as familyFile:
        for line in familyFile:
            line = line.replace("-", "_")
            lineList = line.rstrip("\n").split()
            if lineList[-1] is "2":
                parentDict[f"gts.{lineList[2]}"] = f"gts.{lineList[1]}"
                parentDict[f"gts.{lineList[3]}"] = f"gts.{lineList[1]}"
                patientList.append(f"gts.{lineList[1]}")
                familyDict[f"gts.{lineList[1]}"] = [f"gts.{lineList[2]}", f"gts.{lineList[3]}"]
            else:
                parentList.append(f"gts.{lineList[1]}")

"""
Iterate through inputFile in order to create two dictionaries: sampleGenotype and samplePositions. The key of 
sampleGenotype is the sample ID and value is a dictionary where the key is a gene and the value is a list of all 
genotypes ("0|1", "1|0", or "0|0") for that gene that meet specific CADD score, minor allele frequency, and impact 
severity criteria. The samplePositions has the same information, except the list for each gene is genotype positions.
"""

sampleGenotype = {}
samplePositions = {}
sampleIndexes = []
with open(geminiTsv) as geminiFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, af_1k_index, af_gnomAD_index, lofIndex, exonicIndex, samples = getHeaderInfo(headerList)
    for sample in samples:
        sampleIndexes.append(headerList.index(sample))
        sampleGenotype[sample] = {}
        samplePositions[sample] = {}    
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, af_1k, af_gnomAD, lof, exonic = getLineInfo(lineList)
        if af_gnomAD != "None":
            af = af_gnomAD
        elif af_gnomAD == "None" and af_1k != "None":
            af = af_1k
        else:
            af = "None"
        if cadd != "None" and af != "None" and impactFilterOnly == "n" and exonic == "1":
            if float(cadd) >= inputCadd and float(af) <= inputAF and impact == "HIGH":
                iterateThroughSamples()
        elif impactFilterOnly == "y" and impact == "HIGH" and gene != "None" and exonic == "1":
            iterateThroughSamples()
print("Sample Dictionaries Created.")

"""
Use sampleGenotype and samplePositions to generate a new dictionaries where the key is the sample ID and the value is 
a dictionary where the key is a gene and the value is a list of genotypes (or positions) where CH variant(s) are found.
"""

homAltPositionDict = {}
homAltGenotypeDict = {}

for patient in patientList:
    homAltPositionDict[patient] = {}
    homAltGenotypeDict[patient] = {}
    parent1 = familyDict[patient][0]
    parent2 = familyDict[patient][1]
    for gene, genotypes in sampleGenotype[patient].items():
        for i, genotype in enumerate(genotypes):
            positionList = samplePositions[patient][gene]
            position = positionList[i]
            #This part helps eliminate genotypes being added to the CH list where either parent is homozygous recessive
            parentGenotype1 = ""
            parentGenotype2 = ""
            if gene in samplePositions[parent1] and position in samplePositions[parent1][gene]:
                parentPosIndex = samplePositions[parent1][gene].index(position)
                parentGenotype1 = sampleGenotype[parent1][gene][parentPosIndex]
            if gene in samplePositions[parent2] and position in samplePositions[parent2][gene]:
                parentPosIndex = samplePositions[parent2][gene].index(position)
                parentGenotype2 = sampleGenotype[parent2][gene][parentPosIndex]
            if genotype in ["1|1", "1/1", "2|2", "2/2", "3|3", "3/3"] and parentGenotype1 not in ["1|1", "1/1", "2|2", "2/2", "3|3", "3/3"] and parentGenotype2 not in ["1|1", "1/1", "2|2", "2/2", "3|3", "3/3"]:
                if gene not in homAltPositionDict[patient]:
                    homAltPositionDict[patient][gene] = [position]
                    homAltGenotypeDict[patient][gene] = [genotype]
                elif gene in homAltPositionDict[patient]:
                    homAltPositionDict[patient][gene].append(position)
                    homAltGenotypeDict[patient][gene].append(genotype)
print("CH variant dictionaries created.")

#Iterate through the input file and use the homAltPositionDict in order to output CH variant data for each sample
with open(geminiTsv) as geminiFile, open(outputFile, "w") as outputFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, af_1k_index, af_gnomAD_index, lofIndex, exonicIndex, samples = getHeaderInfo(headerList)
    columnInfo = headerList[0:15]
    newHeader = "\t".join(columnInfo) + "\tgenotype\tsample\n"
    outputFile.write(newHeader)
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, af_1k, af_gnomAD, lof, exonic = getLineInfo(lineList)
        for sampleIndex in sampleIndexes:
            sample = headerList[sampleIndex]
            if sample in homAltPositionDict and gene in homAltPositionDict[sample] and start in homAltPositionDict[sample][gene]:
                genotype = lineList[sampleIndex]
                numericGenotype = getNumericGenotype(genotype, ref, alt)
                if "." not in numericGenotype and numericGenotype in ["1|1", "1/1", "2|2", "2/2", "3|3", "3/3"]:
                    columnInfo = lineList[0:15]
                    columnStr = "\t".join(columnInfo)
                    newLine = f"{columnStr}\t{numericGenotype}\t{sample.replace('gts.', '')}\n"
                    outputFile.write(newLine)

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')