##########################################################
### Import Necessary Modules

import argparse                        #provides options at the command line
import sys                             #take command line arguments and uses it in the script
import gzip                            #allows gzipped files to be read
import re                              #allows regular expressions to be used

##########################################################
### Command-line Arguments

parser = argparse.ArgumentParser(description="A script to find the average and standard deviation of shared alleles using the output from VCFsharedAlleles.v1.0.py.  Remove individual from population map (or make a copy with it removed) if you would like to remove from the analysis.")
parser.add_argument("-file", help = "The file output be VCFsharedAlleles.v1.0.py", default=sys.stdin, required=True)
parser.add_argument("-pop", help = "A population map, format sampleID(matches vcf)<tab>population.  This file is used to help partition the output file.", default=sys.stdin, required=True)
args = parser.parse_args()

#########################################################
class Variables():
    popMap = {}

class OpenFile():
    def __init__ (self, f, typ, fnum):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r')
        if typ == "file":
            sys.stderr.write("\n\tOpened shared allele file: {}\n\n".format(f))
            OpenShared(self.filename,fnum)   
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
                    #sys.stderr.write("{}\t{}\n".format(self.individual, self.pop))
        self.open_pop.close()    
            
class OpenShared():
    def __init__ (self, f, fnum):
        self.open_file = f
        self.sharedAlleles = {}
        self.populations = {}
        self.printWarning = 0
        for line in self.open_file:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')   
            if not re.search("^#", line):
                self.individual1, self.individual2, self.sharedAlleleCount = line.split()[0:3]
                if self.individual1 != self.individual2:
                    if "{}\t{}".format(self.individual1, self.individual2) not in self.sharedAlleles and "{}\t{}".format(self.individual2, self.individual1) not in self.sharedAlleles:
                        self.sharedAlleles["{}\t{}".format(self.individual1, self.individual2)] = 1
                        self.sharedAlleles["{}\t{}".format(self.individual2, self.individual1)] = 1
                        try:
                            self.pop1 = Variables.popMap[self.individual1]
                            self.pop2 = Variables.popMap[self.individual2]
                            if "{}\t{}".format(self.pop2, self.pop1) in self.populations: ### Makes it so that it doesn't have multipe combinations of the same comparison
                                self.pop1 = Variables.popMap[self.individual2]
                                self.pop2 = Variables.popMap[self.individual1]                            
                            if "{}\t{}".format(self.pop1, self.pop2) in self.populations:
                    	         self.populations["{}\t{}".format(self.pop1, self.pop2)] += ",{}".format(self.sharedAlleleCount)
                            else:
                                self.populations["{}\t{}".format(self.pop1, self.pop2)] = "{}".format(self.sharedAlleleCount)
                        except:
                            self.printWarning = "{}\t{}".format(self.individual1, self.individual2)
                            
        if self.printWarning != 0:
            sys.stderr.write("\tWarning: at least one individual not found in the population map file\n\tThe last instance found was with this pair of individuals {}\n\n".format(self.printWarning))
        self.printed = {}                        
        for self.pop in self.populations:
            self.average = 0
            self.stdDev = 0
            self.total = 0
            for self.count in self.populations[self.pop].split(","):
                self.total += 1
                self.average += int(self.count)
            if int(self.total) > 0:
                self.average = float(self.average)/int(self.total)
                self.sqDist = 0
                for self.count in self.populations[self.pop].split(","):
                    self.sqDist += (int(self.count) - float(self.average))**2
                self.stdDev = (self.sqDist/len(self.populations[self.pop].split(",")))**0.5
                print("{}\n\tTotal Comparisons: {}\n\tAverage: {}\n\tStandard Deviation: {}\n".format(self.pop, self.total, self.average, self.stdDev))
            else:
            	sys.stderr.write("Error, didn't find any comparisons\n")
        self.open_file.close()    



if __name__ == '__main__':
    Variables()
    open_pop = OpenFile(args.pop, "pop", 1)
    open_pop = OpenFile(args.file, "file", 1)
