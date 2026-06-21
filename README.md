# HLA Allele To Serotype (HATS)

[![Perl Version](https://img.shields.io/badge/perl-v5.26.3-blue.svg)](https://www.perl.org/)
[![Platform](https://img.shields.io/badge/platform-RHEL%208.6-red.svg)](https://www.redhat.com/)
[![License](https://img.shields.io/badge/License-SBC%20BSD--3--Clause-green.svg)](#license)

## 📌 Introduction
The **HLA Allele To Serotype (HATS)** tool generates various HLA allele-to-serotype conversion mapping files using the reference sequences provided by the **IPD-IMGT/HLA Database**. 

HATS was developed using **Perl 5 (v5.26.3)** on **Red Hat Enterprise Linux 8.6**.

---

## 📂 Output Files & Directories
If you only wish to access the generated HATS datasets without running the full pipeline, you can simply download the standalone `output/` directory.

### Directory Structure Overview
Inside the `output/` directory, files are structured into the following subdirectories:

| Subdirectory | Description |
| :--- | :--- |
| **`RESULTS`** | Untruncated HLA allele-to-serotype equivalency tables for 11 HLA loci. |
| **`TWORESULTS`** | Truncated two-field HLA allele-to-serotype equivalency tables for 11 HLA loci. |
| **`PRACTICAL`** | HLA allele-to-serotype equivalency tables where antigens are reported at the **highest-resolution** level. |
| **`TWOPRACTICAL`** | Two-field HLA allele-to-serotype equivalency tables where antigens are reported at the **highest-resolution** level. |
| **`COMBINED`** | Consolidated files combining all 11-loci tables from the `(TWO)RESULTS` and `(TWO)PRACTICAL` directories (see details below). |
| **`RESIDUES`** | Epitope mapping datasets containing sequence structural identifiers (`_DEP_`, `target_gene_`, and `_LAX_`). |
| **`LEGACY`** | Historical HLA allele-to-serotype equivalency tables matching nomenclature prior to the **2026 HLA Nomenclature updates**. |
| **`COUNT`** | Statistical summaries tracking total alleles/proteins assigned to **FULL**, **SEROTYPE**, **INCOMPLETE**, and **InSilico** categories. |

---

### 🔬 Deep-Dive: Key Combined Output Files

#### 1. `Antigen_Bw4Bw6` File
This file details specific Bw4- and Bw6-associated structural reactivities:
* **HLA-B Bw4 reactivity** is explicitly defined by residues **82L** and **83R** in the HLA-B protein. This identical motif exists within certain HLA-A antigens, which are subsequently cataloged as having Bw4-"associated" reactivity.
* **HLA-B Bw6 reactivity** is defined by residues **76E**, **82R**, and **83G** in the HLA-B protein.
* **HLA-B Bw4/Bw6 "negative" reactivity** is defined by the residues **76V**, **82R**, and **83G**.

#### 2. `combined_A_B_C_DPA1_DPB1_DQA1_DQB1_DRB1_DRB3_DRB4_DRB5_Protein_Antigen_Table_IMGT_HLA_VERSION_DATE.csv`
* **Protein:** Shows the two-field allele name (excludes "Null" and "Questionable" alleles).
* **Qualifier:** Assignment Status markers:
  * **`F` (FULL):** Indicates identical residues at all determining epitope positions (DEPs). Uses the `target_gene_` file definitions. It is viable to use the prototype as a surrogate in virtually all instances. *(e.g., A\*24:07 is serologically equivalent to the associated antigen A2402 with a FULL assignment, belonging to the A24 split / A9 broad antigens).*
  * **`S` (SEROTYPE):** Assigned when proteins do not achieve FULL status but match the 5–10 critical residues historically used to define pre-2026 Nomenclature antigens (found in the `_LAX_` file). Prototype surrogates work in *most* instances, but rigorous donor-specific antibody (DSA) assessment requires closer serum-antibody pattern epitope analysis.
  * **`I` (INCOMPLETE):** Assigned when a protein contains a **single amino acid (AA) mismatch** within the predefined 5–10 `_LAX_` residue positions. Mismatches follow the format `Antigen_position` *(e.g., A\*01:08 is cataloged as `A1_163`, indicating an A1 INCOMPLETE status due to a single mismatch at position 163).*
  * **`U` (UNASSIGNED)**
* **Associated / Split / Broad:** Represents official WHO-assigned antigen equivalencies.
* **CIWD3.0 / CWD2.0 / EURCWD:** Tracks status designations: **C** (Common), **I** (Intermediate), and **WD** (Well-Documented).

#### 3. `combined_A_B_C_DPA1_DPB1_DQA1_DQB1_DRB1_DRB3_DRB4_DRB5_Allele_Antigen_Table_IMGT_HLA_VERSION_DATE.csv`
* **Allele:** Displays all alleles captured in the target IPD-IMGT/HLA Database release version, **including** "Null" and "Questionable" alleles.
* All remaining data columns map identically to the *Protein* file layout.

---

## 💻 Installation & Setup

### Prerequisites
HATS requires a Linux environment (tested natively on Red Hat Enterprise Linux 8.6) with **Perl 5** installed.

### Input Data Preparation
Before running the tool, you must source reference datasets from the [ANHIG/IMGTHLA GitHub Repository](https://github.com/ANHIG/IMGTHLA/) and the IPD-IMGT/HLA Database. 

1. Download the following required files:
   `hla_prot.fasta`, `A_prot.msf`, `B_prot.msf`, `C_prot.msf`, `DRB_prot.msf`, `DPA1_prot.msf`, `DPB1_prot.msf`, `DQA1_prot.msf`, and `DQB1_prot.msf`.
2. Move all downloaded files into `input/` folder in the root directory.
3. *(Optional)* Append the database release version identifier to the primary fasta filename (e.g., saving it explicitly as `input/hla_prot.fasta.3.64.0`).

---

## 🚀 Usage

### Run the Complete Pipeline
To compute serotype mappings for all 11 HLA loci globally in a single invocation:

```bash
./run.sh
# When completed, output files are generated in output directory: RESULTS, TWORESULTS, PRACTICAL, TWOPRACTICAL, COMBINED, LEGACY, RESIDUES, COUNT directories.

### Run Loci Separately
If you prefer to generate files for each locus separately, use the following execution scripts:
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

ℹ️ Locus-Specific Execution Outputs: When HATS is executed for an individual locus, an output/gene_Allele_Residues_DB_Version_Date.csv file is generated. It contains all expressed two-field alleles, assigned serotypes, comments, and specific residues for each key residue position.

    Untruncated HLA allele-to-serotype tables are saved in the RESULTS/ directory.

    Truncated two-field HLA allele-to-serotype tables are saved in the TWORESULTS/ directory.

```

🤝 Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## REFERENCES
Osoegawa K, Marsh SGE, Holdsworth R, Heidt S, Fischer G, Murphey C, Maiers M, Fernández Viňa MA. A new strategy for systematically classifying HLA alleles into serological specificities. HLA. 2022 Sep;100(3):193-231. doi: 10.1111/tan.14662. Epub 2022 Jun 22. PMID: 35538616.
Osoegawa K, Yim K, Jeracki M, Nguyen TN, Wang L, Cho A, David R, Son J, Mankey A, Marsh SGE, Gendzekhadze K, Murphey C, Fernández Viňa MA. A new strategy for systematically classifying HLA alleles into serological specificities: Update and refinement. HLA. 2024 Oct;104(4):e15702. doi: 10.1111/tan.15702. PMID: 39435845.
Osoegawa K, Son J, Yim K, Marsh SGE, Fernández Viňa MA. Replacements at Structural or Functional Dimorphisms 103, 109 and 167 Distinguish HLA Class I Serologically Defined Antigens. HLA. 2025 Sep;106(3):e70387. doi: 10.1111/tan.70387. PMID: 40944463; PMCID: PMC12432678.
Marsh SGE, Osoegawa K, Bodmer WF, Bontrop RE, Carrington MN, Erlich HA, Heidt S, Holdsworth R, Mayr WR, Maiers M, Parham P, Petersdorf EW, Robinson J, Trowsdale J, Fernández-Viña M. Nomenclature for Factors of the HLA System, 2026. HLA. 2026 Mar;107(3):e70595. doi: 10.1111/tan.70595. PMID: 41742599; PMCID: PMC12936402.

## input directory:
# The following files: hla_prot.fasta, A_prot.msf, B_prot.msf, C_prot.msf, DRB_prot.msf, DPA1_prot.msf, DPB1_prot.msf, DQA1_prot.msf and DQB1_prot.msf are downloaded from https://github.com/ANHIG/IMGTHLA/ and saved in input directory for test.

#IPD-IMGT/HLA Database reference:
Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE. IPD-IMGT/HLA Database. Nucleic Acids Res. 2020 Jan 8;48(D1):D948-D955. doi: 10.1093/nar/gkz950. PMID: 31667505; PMCID: PMC7145640.

#CWD2.0, CIWD3.0 and European CWD catalogues were downloaded from the following literatures, formatted and saved in CWD2, CIWD and EURCWD directories, respectively.

#CWD 2.0 catalogue:
Mack SJ, Cano P, Hollenbach JA, He J, Hurley CK, Middleton D, Moraes ME, Pereira SE, Kempenich JH, Reed EF, Setterholm M, Smith AG, Tilanus MG, Torres M, Varney MD, Voorter CE, Fischer GF, Fleischhauer K, Goodridge D, Klitz W, Little AM, Maiers M, Marsh SG, Müller CR, Noreen H, Rozemuller EH, Sanchez-Mazas A, Senitzer D, Trachtenberg E, Fernandez-Vina M. Common and well-documented HLA alleles: 2012 update to the CWD catalogue. Tissue Antigens. 2013 Apr;81(4):194-203. doi: 10.1111/tan.12093. PMID: 23510415; PMCID: PMC3634360.

#CIWD 3.0 catalogue:
Hurley CK, Kempenich J, Wadsworth K, Sauter J, Hofmann JA, Schefzyk D, Schmidt AH, Galarza P, Cardozo MBR, Dudkiewicz M, Houdova L, Jindra P, Sorensen BS, Jagannathan L, Mathur A, Linjama T, Torosian T, Freudenberger R, Manolis A, Mavrommatis J, Cereb N, Manor S, Shriki N, Sacchi N, Ameen R, Fisher R, Dunckley H, Andersen I, Alaskar A, Alzahrani M, Hajeer A, Jawdat D, Nicoloso G, Kupatawintu P, Cho L, Kaur A, Bengtsson M, Dehn J. Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0. HLA. 2020 Jun;95(6):516-531. doi: 10.1111/tan.13811. Epub 2020 Jan 31. PMID: 31970929; PMCID: PMC7317522.

#European CWD catalogue:
Sanchez-Mazas A, Nunes JM, Middleton D, Sauter J, Buhler S, McCabe A, Hofmann J, Baier DM, Schmidt AH, Nicoloso G, Andreani M, Grubic Z, Tiercy JM, Fleischhauer K. Common and well-documented HLA alleles over all of Europe and within European sub-regions: A catalogue from the European Federation for Immunogenetics. HLA. 2017 Feb;89(2):104-113. doi: 10.1111/tan.12956. PMID: 28102034.

## License
© 2022 Stanford Blood Center L.L.C.

Please cite this work as:
Osoegawa, Kazutoyo et al. “A new strategy for systematically classifying HLA alleles into serological specificities.” HLA vol. 100,3 (2022): 193-231. doi:10.1111/tan.14662

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
