# SharedAllelesVCF
A script to identify shared alleles between all the individuals in a vcf file.  For example, if one individual is a heterozygote and the other individuals is a homozygote, they will share one out of two alleles (assuming they are both diploid) at a locus.  This script will count the combined shared alleles for each all loci in the vcf file.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

This script has been tested with Python 2.7 and 3 and should work with either.
The script requires a vcf file.  The vcf file can be compressed with gzip.  This script also requires a population map (individual<tab>population).  The individual names must match the names in the vcf files

<!-- usage -->
## Usage

Find Pairwise:
python VCFsharedAlleles.v1.0.py -vcf file.vcf -pop populationMap.txt > pairwiseComparison.txt

To see the usage and get futher information: python VCFsharedAlleles.v1.0.py -h

See population values:
python AnalyzeSharedAlleles.v1.0.py -file pairwiseComparison.txt -pop populationMap.txt

To see the usage and get further information: python AnalyzeSharedAlleles.v1.0.py -h

<!-- license -->
## License 

Distributed under the MIT License.
