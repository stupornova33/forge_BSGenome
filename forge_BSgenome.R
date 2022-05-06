if (!require("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!require("BSgenome", quietly = TRUE)) BiocManager::install("BSgenome")
if (!require("GenomeInfoDb", quietly = TRUE)) BiocManager::install("GenomeInfoDb")
if (!require("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
library(GenomeInfoDb)
library(Biostrings)
library(BSgenome)
library(rtracklayer)
### This code allows the user to forge a BSGenome out of any genome available through NCBI
### Download GCF_003640425.2_ASM364042v2_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/640/425/GCF_003640425.2_ASM364042v2/
### modify path
dna <- readDNAStringSet("C:/Users/tjarva/Desktop/RTestFiles/GCF_016746395.2_Prin_Dsim_3.1_genomic.fa")

### extract genome info
seqinfo <- getChromInfoFromNCBI("GCF_016746395.2",
                                assembled.molecules.only=FALSE,
                                assembly.units=NULL,
                                recache=FALSE,
                                as.Seqinfo=TRUE) 


### Check seq names
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))


chrominfo <- getChromInfoFromNCBI("GCF_016746395.2")
expected_RefSeqAccn <- chrominfo[ ,"RefSeqAccn"]

### check to see if the names match between the two objects
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### In this case the expected RefSeq names did not match the actual names. Any mismatched names should be resolved. Since this was a
### non-existent sequence, I simply removed the name from chrominfo
### This step is not necessary if the names match
chrominfo <- chrominfo[-c(7),]

### check again
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences. An alternative would be to rename them to
### chrominfo[ , "SequenceName"] but these names are VERY ugly (e.g.
### "ScRZk8e_1;HRSCAF=1").
names(dna) <- expected_RefSeqAccn

### Export as 2bit.
export.2bit(dna, "GCF_016746395.2.sorted.2bit")

file <- "GCF_016746395.2_Prin_Dsim_3.1_genomic.fa"
### Define path to seed file
### See example_seed_file.txt for format of seed file 
seed_files <- "C:/Users/tjarva/Desktop/RTestFiles/BSgenome.Dsimulans.NCBI.v3.1-seed"

### forge the genome package
### This will create a directory containing the genome package wherever your R library is
forgeBSgenomeDataPkg("BSgenome.Dsimulans.NCBI.v3.1-seed")

### You can now load the library whenever you wish by using
### library(BSgenome.Dsimulans.NCBI.v3.1-seed)




