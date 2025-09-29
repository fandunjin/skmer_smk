# Skmer snakemake workflow
<span style="font-size:12px">&nbsp;&nbsp;&nbsp;&nbsp;This Snakemake-based workflow provides an automated pipeline for phylogenetic analysis starting from raw sequencing data, integrating multiple bioinformatics tools including Bowtie2 for nuclear genome filtering, BBTools for sequence repair and merging, Skmer for k-mer-based distance matrix calculation, and FastME/RAxML for robust phylogenetic tree construction. The workflow efficiently processes multiple samples through parallel computing, ensures reproducibility with containerized environments, and delivers comprehensive outputs including distance matrices and consensus trees for evolutionary studies.
## Test Procedure (Using Test Data)</span>
```bash
git@github.com:fandunjin/skmer_smk.git
conda create -n skmer snakemake Jellyfish Mash seqtk pandas=1.5.2 scipy biopython
conda activate skmer
git clone https://github.com/shahab-sarmashghi/Skmer.git
cd Skmer
python setup.py install
cd ../
snakemake --core 4
```

# How to analyse your own data with Skmer workflow?
First, clone this repository
```bash
git@github.com:fandunjin/skmer_smk.git
```
Next, create the skmer working environment
```bash
conda create -n skmer snakemake Jellyfish Mash seqtk pandas=1.5.2 scipy biopython
# In skmer correct, the pandas version must be <2; otherwise, a rep91 ValueError will occur, so an older pandas version needs to be installed separately.
conda activate skmer #get into environment
pip show pandas # confirm the version of pandas
```
Since the version of Skmer available in Conda is outdated, it needs to be manually downloaded from GitHub: https://github.com/shahab-sarmashghi/Skmer
```bash
git clone https://github.com/shahab-sarmashghi/Skmer.git
cd Skmer
python setup.py install
```
Finally, simply run Snakemake to start the analysis
```bash
#place the reference genome file (ref.fna) into the raw_data/ directory, and put the paired-end sequencing reads into raw_data/raw_data/.
snakemake --core 48 --latency-wait 120 #Return to the directory containing the Snakefile and run the following command to complete the analysis
```
### NOTICE
1.Naming must be consistent
Reference genome:ref_fna
Paired-end reads: sample.R1.fq.gz and sample.R2.fq.gz (the "sample" prefix can be changed, but the suffix must be exactly ".R1.fq.gz" and ".R2.fq.gz")
2.Place the data in the correct locations: Place ref.fna in the raw_data/ directory, and put sample.R1.fq.gz (along with sample.R2.fq.gz) inside raw_data/raw_data/
```text
skmer_smk/
├── raw_data/
│   ├── ref.fna
│   └── raw_data/
│       ├── sample1.R1.fq.gz
│       ├── sample1.R2.fq.gz
│       ├── sample2.R1.fq.gz
│       ├── sample2.R2.fq.gz
│       └── ...
├── Snakefile
├── tsv_to_phymat.sh
└── merge_consensus.py
```

## Output results include
Main output files
```text
results/skmer/dimtrx_main_cor_OK.txt.tre：Branch tree
results/skmer/RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS_fixed.tre: Bootstrap tree
results/skmer/integration.tre：Integrate branch and bootstrap tree
```
Other output files
```text
results/skmer/dimtrx_main.txt：Main distance matrix
results/{sample}/nDNA/：Per-sample DNA filtering results
results/{sample}/nDNAOK/：Sequences after head 25,000,000
results/skmer/subsample/：bootstrap Bootstrap replicate results
ref_dir/：Collection directory of sequences for Skmer runs
```
