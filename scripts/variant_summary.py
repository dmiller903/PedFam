import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Summarizes to indicate which genes have variants that are unique to each sample.")

parser.add_argument('input_file_variants', help='Variant Gemini tsv file')
parser.add_argument('input_file_dels', help='Dels Gemini tsv file')
parser.add_argument('input_file_svs', help='SVS Gemini tsv file')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('fam_file', help='A .fam file is used to determine relationships among samples')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile_variants = args.input_file_variants
inputFile_dels = args.input_file_dels
inputFile_svs = args.input_file_svs
outputFile = args.output_file
familyFile = args.fam_file
fileList = [inputFile_variants, inputFile_dels, inputFile_svs]

# Function to find and output unique genes
def getUniqueVariants(file):
    sampleVariants = {}
    sampleGenes = {}
    with open(file) as dataFile:
        headerList = dataFile.readline().rstrip("\n").split("\t")
        dataList = headerList[0:15]
        sampleList = headerList[15:]
        sampleIndexes = []
        for sample in sampleList:
            sampleIndexes.append(headerList.index(sample))
        geneIndex = dataList.index("gene")
        for line in dataFile:
            lineList = line.rstrip("\n").split("\t")
            for i, sample in enumerate(sampleList):
                genotype = lineList[sampleIndexes[i]]
                sampleName = sampleList[i]
                gene = lineList[geneIndex]
                if "." not in genotype:
                    geneList = gene.split("&")
                    variantAdded = False
                    for gene in geneList:
                        if sampleName not in sampleGenes and gene != "None":
                            sampleGenes[sampleName] = {gene}
                            sampleVariants[sampleName] = 1
                            variantAdded = True
                        elif sampleName in sampleGenes and gene != "None" and variantAdded == False:
                            sampleGenes[sampleName].add(gene)
                            sampleVariants[sampleName] += 1
                            variantAdded = True
    return(sampleVariants, sampleGenes)

def summarize_data(variantDict, geneDict):
    for sample, value in variantDict.items():
        outputFile.write(f"{sample}\t{value}\tvariants\n")
    for sample, geneList in geneDict.items():
        outputFile.write(f"{sample}\t{len(geneList)}\tgenes\n")

# Create a dictionary of family members
familyDict = {}
sampleDisease = {}
with open(familyFile) as familyFile:
    for line in familyFile:
        line = line.replace("-", "_")
        lineList = line.rstrip("\n").split()
        if lineList[-1] is "2" and lineList[0] not in familyDict:
            familyDict[f"{lineList[0]}"] = [f"{lineList[1]}"]
            sampleDisease[lineList[1]] = lineList[0]
        elif lineList[-1] is "2" and lineList[0] in familyDict:
            familyDict[f"{lineList[0]}"].append(f"{lineList[1]}")
            sampleDisease[lineList[1]] = lineList[0]


# Use getUniqueVariants to write unique genes to output file
with open(outputFile, "w") as outputFile:
    sampleVariants, sampleGenes = getUniqueVariants(inputFile_variants)
    outputFile.write(f"# Variant Data\n")
    print(sampleVariants)
    print(sampleGenes)
    summarize_data(sampleVariants, sampleGenes)
    sampleVariants, sampleGenes = getUniqueVariants(inputFile_dels)
    outputFile.write(f"# Deletion Data\n")
    summarize_data(sampleVariants, sampleGenes)
    sampleVariants, sampleGenes = getUniqueVariants(inputFile_svs)
    outputFile.write(f"# Structural Variant Data\n")
    summarize_data(sampleVariants, sampleGenes)