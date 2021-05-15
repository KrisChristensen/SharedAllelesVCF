##########################################################
### Import Necessary Modules

import argparse                        #provides options at the command line
import sys                             #take command line arguments and uses it in the script
import gzip                            #allows gzipped files to be read
import re                              #allows regular expressions to be used

##########################################################
### Command-line Arguments

parser = argparse.ArgumentParser(description="A script to identify the number of shared alleles between all individuals in vcf file (not written for phased genotypes)")
parser.add_argument("-vcf", help = "A vcf file that has multiple individuals and bi-allelic SNPs, assumes all filtering has already been completed", default=sys.stdin, required=True)
parser.add_argument("-pop", help = "A population map, format sampleID(matches vcf)<tab>population.  This file is used to help partition the output file.", default=sys.stdin, required=True)
args = parser.parse_args()

#########################################################
### Variables
class Variables():
    popMap = {}


class OpenFile():
    def __init__ (self, f, typ, fnum):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r')
        if typ == "vcf":
            sys.stderr.write("\n\tOpened vcf file: {}\n\n".format(f))
            OpenVCF(self.filename,fnum)
        elif typ == "pop":
            sys.stderr.write("\n\tOpened popMap file: {}\n\n".format(f))
            OpenPop(self.filename,fnum)        
            
class OpenPop():
    def __init__ (self, f, fnum):
        self.open_pop = f
        for line in self.open_pop:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')   
            if not re.search("^#", line):
                self.individual, self.pop = line.split()[0:2]
                if self.individual in Variables.popMap:
                    sys.stderr.write("\tError: {} found more than once in population map file.  Please fix before resubmitting.\n".format(self.individual))
                    exit()
                else:
                    Variables.popMap[self.individual] = self.pop 
        self.open_pop.close()    

class OpenVCF():
    def __init__ (self, f, fnum):
        self.open_vcf = f
        self.lineCount = 0
        self.allIndividuals = []
        self.sharedAlleles = {}
        for line in self.open_vcf:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')  
            if re.search("^#CHROM", line):
                self.allIndividuals = line.split()[9:]
                for self.ind1 in self.allIndividuals:
                    for self.ind2 in self.allIndividuals:
                        self.sharedAlleles["{}\t{}".format(self.ind1, self.ind2)] = 0   
            elif not re.search("^#", line):  
                self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format = line.split()[0:9]
                self.genotypes = line.split()[9:]
                for self.index1, self.indGeno1 in enumerate(self.genotypes):
                    self.geno1 = self.indGeno1.split(":")[0]
                    for self.index2, self.indGeno2 in enumerate(self.genotypes):
                        self.geno2 = self.indGeno2.split(":")[0]
                        self.match = 0
                        if self.geno1 != "./." and self.geno2 != "./.": 
                            if self.geno1 == self.geno2:
                                self.match = 2 
                            elif self.geno1 == "0/1" or self.geno1 == "1/0" or self.geno2 == "0/1" or self.geno2 == "1/0":
                                self.match = 1
                            self.sharedAlleles["{}\t{}".format(self.allIndividuals[self.index1], self.allIndividuals[self.index2])] += int(self.match)
                self.lineCount += 1        
                if self.lineCount % 10000 == 0:
                    sys.stderr.write("\tProcessed {} variants\n".format(self.lineCount))
        sys.stderr.write("\tProcessed {} total variants\n\n".format(self.lineCount))
        for self.combo in self.sharedAlleles:
            print("{}\t{}".format(self.combo, self.sharedAlleles[self.combo]))
        self.open_vcf.close()

if __name__ == '__main__':
    Variables()
    open_pop = OpenFile(args.pop, "pop", 1)
    open_vcf = OpenFile(args.vcf, "vcf", 1)
