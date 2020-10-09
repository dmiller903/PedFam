import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Summarizes to indicate which genes have variants that are unique to each sample.")

parser.add_argument('input_file_CH', help='CH tsv file that is the output of "add_GDI_and_gene_lengths.py"')
parser.add_argument('input_file_homAlt', help='homAlt tsv file that is the output of "add_GDI_and_gene_lengths.py"')
parser.add_argument('input_file_deNovo', help='deNovo tsv file that is the output of "add_GDI_and_gene_lengths.py"')
parser.add_argument('input_file_het', help='het tsv file that is the output of "add_GDI_and_gene_lengths.py"')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('fam_file', help='A .fam file is used to determine relationships among samples')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile_CH = args.input_file_CH
inputFile_homAlt = args.input_file_homAlt
inputFile_deNovo = args.input_file_deNovo
inputFile_het = args.input_file_het
outputFile = args.output_file
familyFile = args.fam_file
fileList = [inputFile_CH, inputFile_homAlt, inputFile_deNovo, inputFile_het]

# Function to find and output unique genes
def getUniqueVariants(file):
    variantType = file.rstrip("\n").split("/")[-1]
    variantType = variantType.rstrip("\n").split("_")[1]

    sampleGenes = {}
    with open(file) as dataFile:
        headerList = dataFile.readline().rstrip("\n").split("\t")
        sampleIndex = headerList.index("sample")
        geneIndex = headerList.index("gene")
        for line in dataFile:
            lineList = line.rstrip("\n").split("\t")
            sample = lineList[sampleIndex]
            gene = lineList[geneIndex]
            disease = sampleDisease[sample]
            if sample not in sampleGenes and gene != "None":
                sampleGenes[sample] = {gene}
            elif sample in sampleGenes and gene != "None":
                sampleGenes[sample].add(gene)

    diseaseGenes = {}
    for sample, genes in sampleGenes.items():
        disease = sampleDisease[sample]
        if disease not in diseaseGenes:
            diseaseGenes[disease] = []
        for gene in genes:
            diseaseGenes[disease].append(gene)

    uniqueGenes = {}
    for sample, genes in sampleGenes.items():
        disease = sampleDisease[sample]
        uniqueGenes[sample] = []
        for gene in genes:
            if diseaseGenes[disease].count(gene) == 1:
                uniqueGenes[sample].append(gene)
    with open(outputFile, 'a') as output:
        for sample, genes in sorted(uniqueGenes.items()):
            disease = sampleDisease[sample]
            if len(genes) > 0:
                for gene in genes:
                    output.write(f"{sample}\t{disease}\t{gene}\t{variantType}\n")

# Create a dictionary of family members
familyDict = {}
sampleDisease = {}
with open(familyFile) as familyFile:
    for line in familyFile:
        line = line.replace("-", "_")
        lineList = line.rstrip("\n").split()
        sampleDisease[lineList[1]] = lineList[0]
        if lineList[-1] is "2" and lineList[0] not in familyDict:
            familyDict[f"{lineList[0]}"] = [f"{lineList[1]}"]
        elif lineList[-1] is "2" and lineList[0] in familyDict:
            familyDict[f"{lineList[0]}"].append(f"{lineList[1]}")

# Write header to output file
with open(outputFile, 'w') as output:
    output.write(f"sample\tdisease\tunique_genes\tvariant_type\n")

# Use getUniqueVariants to write unique genes to output file
for file in fileList:
    print(file)
    getUniqueVariants(file)