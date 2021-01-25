import gzip
from sys import argv

#Create variables of each argument from argv
fileList = argv[1:-2]
outputPath = argv[-2].rstrip("/")
summaryFile = argv[-1]

summaryDict = {}
def getStats(file):
    lineCount = 0
    passed = 0
    phased = 0
    numMulti = 0
    fileName = file.split("/")[-1]
    sample = file.split("/")[-3]
    outputName = f"{outputPath}/{sample}_{fileName}"
    with gzip.open(file, 'rt') as inputFile, gzip.open(outputName, 'wb') as outputFile, open(summaryFile, "a") as summary:
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
                if qualCol == ".":
                    lineCount += 1
                    continue
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

        summary.write(f"{sample}\t{passed}\t{phased}\t{lineCount}\n")                
        print(f"{passed} out of {lineCount} variants passed ({(passed / lineCount) * 100}%)\n\
        for {file}. Of those that passed, {phased} were phased ({(phased / passed) * 100}%).\n\
        There were {numMulti} that passed but were not included in outputfile and were not\n\
        included in the number of passed variants.")

with open(summaryFile, "w") as summary:
    summary.write(f"sample\tnum_passed\tnum_passed_phased\ttotal_variants\n")

for file in fileList:
    getStats(file)