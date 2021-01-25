import glob
import gzip
import concurrent.futures
import argparse
from sys import argv

#Create variables of each argument from argv
fileList = argv[1:-1]
outputPath = argv[-1].rstrip("/")

def getStats(file):
    lineCount = 0
    passed = 0
    phased = 0
    numMulti = 0
    fileName = file.split("/")[-1]
    sample = file.split("/")[-3]
    outputName = f"{outputPath}/{sample}_{fileName}"
    with gzip.open(file, 'rt') as inputFile, gzip.open(outputName, 'wb') as outputFile:
        for line in inputFile:
            if line.startswith("##"):
                outputFile.write(line.encode())
            elif line.startswith("#CHROM"):
                headerList = line.rstrip("\n").split("\t")
                filterIndex = headerList.index("FILTER")
                altIndex = headerList.index("ALT")
                qualIndex = headerList.index("QUAL")
                outputFile.write(line.encode())
            else:
                lineList = line.rstrip("\n").split("\t")
                filterCol = lineList[filterIndex]
                altCol = lineList[altIndex]
                qualCol = lineList[qualIndex]
                genotype = lineList[-1]
                genotype = genotype.split(":")[0]
                if filterCol == "PASS" and "," not in altCol and float(qualCol) >= 20.0:
                    outputFile.write(line.encode())
                    lineCount += 1
                    passed += 1
                    if "|" in genotype:
                        phased += 1
                elif filterCol == "PASS" and "," in altCol and float(qualCol) < 20.0:
                    lineCount += 1
                    numMulti += 1
                else:
                    lineCount += 1
    print(f"{passed} out of {lineCount} variants passed ({(passed / lineCount) * 100}%)\n\
    for {file}. Of those that passed, {phased} were phased ({(phased / passed) * 100}%).\n\
    There were {numMulti} that passed but were not included in outputfile and were not\n\
    included in the number of passed variants.")

with concurrent.futures.ProcessPoolExecutor(max_workers=44) as executor:
        executor.map(getStats, fileList)