# Introduction

HLA Allele To Serotype (HATS) was developed with Perl 5.10 on CentOS 7. The program was developed using the reference sequences based on IPD-IMGT/HLA Database release version 3.38.0.

## Installation
Download HATS directory in a Linux computer


## Usage

```perl

# Create an output directory
mkdir output

# Input file
# Download hla_prot.fasta file from IPD-IMGT/HLA Database. IPD-IMGT/HLA Database release version 3.38.0 must be added at the end of hla_prot.fasta file and saved in input directory
# input/hla_prot.fasta.3.44.0
# An empty file is created to show database version when one of the following commands are executed

hla_prot.fasta.3.44.0

# Generate HLA-A allele to serotype table
./runHlaA.pl

# Generate HLA-B allele to serotype table
./runHlaB.pl

# Generate HLA-C allele to serotype table
./runHlaC.pl

# Generate HLA-DRB1 allele to serotype table
./runDRB1.pl

# Generate HLA-DRB3 allele to serotype table
./runDRB3.pl

# Generate HLA-DRB4 allele to serotype table
./runDRB4.pl

# Generate HLA-DRB5 allele to serotype table
./runDRB5.pl

# Generate HLA-DQB1 allele to serotype table
./runDQB1.pl

# Generate HLA-DQA1 allele to serotype table
./runDQA1.pl

# Generate HLA-DPB1 allele to serotype table
./runDPB1.pl

# Generate HLA-DPA1 allele to serotype table
./runDPA1.pl

# Untruncated HLA allele to serotype table is saved in RESULTS directory
# Truncated two-field HLA allele to serotype table is saved in TWORESULTS directory

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## REFERENCES
#Publication: https://pubmed.ncbi.nlm.nih.gov/35538616/

#input/hla_prot.fasta file was downloaded from https://github.com/ANHIG/IMGTHLA/ and saved in input directory for test.

#IPD-IMGT/HLA Database:
https://pubmed.ncbi.nlm.nih.gov/31667505/

#CWD2.0, CIWD3.0 and European CWD catalogues were downloaded from the following literatures, formatted and saved in CWD2, CIWD and EURCWD directories, respectively.

#CWD 2.0 catalogue:
https://pubmed.ncbi.nlm.nih.gov/23510415/

#CIWD 3.0 catalogue:
https://pubmed.ncbi.nlm.nih.gov/31970929/

#European CWD catalogue:
https://pubmed.ncbi.nlm.nih.gov/28102034/

## License
© 2022 Stanford Blood Center L.L.C.

Please cite this work as:
Osoegawa, Kazutoyo et al. “A new strategy for systematically classifying HLA alleles into serological specificities.” HLA vol. 100,3 (2022): 193-231. doi:10.1111/tan.14662

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
