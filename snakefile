rep_n = 100
SAMPLES, = glob_wildcards("raw_data/raw_data/{sample}.R1.fq.gz")
def get_all_fq(wildcards):
    return expand("results/{sample}/nDNAOK/{sample}.fq", sample=SAMPLES)

rule all: 
    input:
        expand("results/{sample}/nDNAOK/{sample}.fq", sample=SAMPLES),
        "results/skmer/dimtrx_main.txt", # 距离矩阵
        expand("results/skmer/subsample/rep{i}/dimtrx_rep.txt", i=range(rep_n)),
        "results/skmer/dimtrx_main_cor_OK.txt.tre", # 最终树
        "results/skmer/RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS_fixed.tre", # skmer树
	    "results/skmer/integration.tre"

rule bowtie2_build: 
    input: ref="raw_data/ref.fna" 
    output: expand("raw_data/ref.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]) 
    shell: "bowtie2-build {input.ref} raw_data/ref" 

rule bowtie2_filter:
    input:
        ref=expand("raw_data/ref.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
        r1="raw_data/raw_data/{sample}.R1.fq.gz",
        r2="raw_data/raw_data/{sample}.R2.fq.gz"
    output:
        un_conc_gz_r1="results/{sample}/nDNA/un-conc-mate.1",
        un_conc_gz_r2="results/{sample}/nDNA/un-conc-mate.2",
        sam="results/{sample}/nDNA/{sample}.sam"
    threads: 48
    shell:
        "mkdir -p results &&"
        "bowtie2 -x raw_data/ref -1 {input.r1} -2 {input.r2} --un-conc-gz results/{wildcards.sample}/nDNA -S {output.sam} -p {threads}"

rule rename_unconc:
    input:
        r1="results/{sample}/nDNA/un-conc-mate.1",
        r2="results/{sample}/nDNA/un-conc-mate.2"
    output:
        r1="results/{sample}/nDNA/un-conc-mate.1.fq.gz",
        r2="results/{sample}/nDNA/un-conc-mate.2.fq.gz"
    shell:
        "mv {input.r1} {output.r1} && mv {input.r2} {output.r2}"

rule repair_pairs:
    input:
        r1="results/{sample}/nDNA/un-conc-mate.1.fq.gz",
        r2="results/{sample}/nDNA/un-conc-mate.2.fq.gz"
    output:
        r1="results/{sample}/nDNA/un-conc-mate.1.fixed.fq.gz",
        r2="results/{sample}/nDNA/un-conc-mate.2.fixed.fq.gz",
        single="results/{sample}/nDNA/singletons.fq.gz"
    threads: 4
    shell:
        "repair.sh in1={input.r1} in2={input.r2} "
        "out1={output.r1} out2={output.r2} outs={output.single} -Xmx4g threads={threads}"

rule bbmerge: 
    input:
        r1="results/{sample}/nDNA/un-conc-mate.1.fixed.fq.gz",
        r2="results/{sample}/nDNA/un-conc-mate.2.fixed.fq.gz" 
    output: merged="results/{sample}/sample_merged.fq", unmerged1="results/{sample}/sample_unmerged1.fq", unmerged2="results/{sample}/sample_unmerged2.fq" 
    shell: "bbmerge.sh t=2 in1={input.r1} in2={input.r2} out={output.merged} outu1={output.unmerged1} outu2={output.unmerged2}" 

rule concatenate_reads: 
    input: "results/{sample}/sample_merged.fq", "results/{sample}/sample_unmerged1.fq", "results/{sample}/sample_unmerged2.fq" 
    output: "results/{sample}/sample.fq" 
    shell: "cat {input} > {output}" 

rule head_reads: 
    input: "results/{sample}/sample.fq" 
    output: "results/{sample}/nDNAOK/{sample}.fq" 
    shell: "head -25000000 {input} > {output}" 
    
import glob
import os
rule link_to_ref_dir:
    input: get_all_fq
    output: touch("ref_dir/.moved") # 全局标记 
    params: ref_dir="ref_dir" 
    shell: """ 
        mkdir -p {params.ref_dir} 
        for fq in {input}; do sample=$(basename $fq .fq) 
        cp $(realpath $fq) {params.ref_dir}/${{sample}}.fq 
        done 
        """ 

rule skmer_reference:
    input: 
        moved="ref_dir/.moved"
    output: 
        temp("results/skmer/dimtrx_main.txt.txt")
    params:
        ref_dir = "ref_dir",
        prefix  = "results/skmer/dimtrx_main.txt"
    threads: 4
    shell:
        "mkdir -p results/skmer && "
        "skmer reference {params.ref_dir} -s 100000 -S 42 -p {threads} -t -o {params.prefix}"

rule rename_main: 
    input: "results/skmer/dimtrx_main.txt.txt" 
    output: "results/skmer/dimtrx_main.txt" 
    shell: "mv {input} {output}" 

rule skmer_subsample:
    input:
        main="results/skmer/dimtrx_main.txt"
    output:
        expand("results/skmer/subsample/rep{i}/dimtrx_rep.txt", i=range(rep_n))
    shell:
        """
        skmer subsample ref_dir \
            -sub results/skmer/subsample \
            -b 100 -i 0 \
            -s 100000 -S 42 -t
        """

rule skmer_correct:
    input:
        main="results/skmer/dimtrx_main.txt",
        sub=expand("results/skmer/subsample/rep{i}/dimtrx_rep.txt", i=range(rep_n))
    output:
        "results/skmer/dimtrx_main_cor_.txt"
    shell:
        """
        skmer correct -main {input.main} \
                      -sub results/skmer/subsample
        """

rule tsv_to_phymat: 
    input: "results/skmer/dimtrx_main_cor_.txt" 
    output: "results/skmer/dimtrx_main_cor_OK.txt" 
    shell: "bash tsv_to_phymat.sh {input} {output}" 

rule fastme: 
    input: "results/skmer/dimtrx_main_cor_OK.txt" 
    output: "results/skmer/dimtrx_main_cor_OK.txt.tre" 
    shell: "fastme -i {input} -o {output}" 

rule tsv_to_phymat_sub:
    input:  "results/skmer/subsample/rep{idx}/dimtrx_rep.txt"
    output: "results/skmer/subsample/rep{idx}/dimtrx_rep_cor0.txt"
    shell:
        "bash tsv_to_phymat.sh {input} {output}"

rule fastme_per_sub:
    input:  "results/skmer/subsample/rep{idx}/dimtrx_rep_cor0.txt"
    output: "results/skmer/subsample/rep{idx}/dimtrx_rep_cor.txt.tre"
    shell:  "fastme -i {input} -o {output}"

rule concatenate_trees:
    input:  expand("results/skmer/subsample/rep{i}/dimtrx_rep_cor.txt.tre", i=range(rep_n))
    output: "results/skmer/bootstrap.All_consensus"
    shell:  "cat {input} > {output}"
    
rule raxml:
    input:
        "results/skmer/bootstrap.All_consensus"
    output:
        "RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS"
    shell:
        "raxmlHPC -J MRE -z {input} -p 4424 -m GTRCAT -n BS_TREE_CONS"

rule tree_fixed:
    input:  "RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS"
    output: "results/skmer/RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS_fixed.tre"
    run:
        with open(input[0], 'r', encoding='utf-8') as f_in, \
             open(output[0], 'w', encoding='utf-8') as f_out:
            import re
            content = f_in.read()
            pattern = r'(:[0-9]+\.[0-9]+)\[([0-9]+)\]'
            fixed_content = re.sub(pattern, r'\2\1', content)
            f_out.write(fixed_content)

rule integrate:
    input: "results/skmer/RAxML_MajorityRuleExtendedConsensusTree.BS_TREE_CONS_fixed.tre"
    output: "results/skmer/integration.tre"
    shell: "python merge_consensus.py {input} results/skmer/dimtrx_main_cor_OK.txt.tre {output}"
