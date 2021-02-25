import gzip
import time

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

def create_cadd_dict(currentChr):
    tempDict = {}
    with gzip.open(f"cadd/cadd_chr{currentChr}.tsv.gz", "rt") as caddFile:
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
    with gzip.open(f"clinVar/clinVar_chr{currentChr}.tsv.gz", "rt") as clinVarFile:
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
    fileName = f"maf/1K_chr{currentChr}.tsv.gz"
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
    with gzip.open(f"gnomAD/gnomAD_chr{currentChr}_AF.tsv.gz", "rt") as gnomAD:
        if currentChr == "1":
            with gzip.open("gnomAD_1K_cadd_GRCh38.tsv.gz", "wb") as output:
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
            with gzip.open("gnomAD_1K_cadd_GRCh38.tsv.gz", "ab") as output:
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