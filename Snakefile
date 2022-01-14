import datetime
configfile: "config.yaml"
wildcard_constraints:
    species = "\w+"

rule all:
    input:
        "data/panTro5_hg38_figs/"

def getLiftoverFile(wildcards):
    return config["liftover"][wildcards.species]

def getCorrectBed(wildcards):
    if wildcards.species == config["speciesA"]:
        return "data/intersected.bed"
    else:
        return "data/common" + config["speciesA"] + "Lift{species}.bed"

def getCorrectFasta(wildcards):
    if wildcards.species == config["speciesA"]:
        return "data/commonSeqs_" + config["speciesB"] + ".fa"
    else:
        return "data/commonSeqs_" + config["speciesA"] + ".fa"

def correctWildcard(wildcards):
    return wildcards.species[0].upper() + wildcards.species[1:]

def getBarrierFile(wildcards):
    return config["barriers"][wildcards.ref]

rule getNonOverlapInt:
    shadow: "shallow"
    output:
        bedHO = "data/{ref}.lift.{species}.{ref}.overlap.bed",
        bedHNO = "data/{ref}.lift.{species}.{ref}.nonOverlap.bed",
        bedSO = "data/{ref}.lift.{species}.{species}.overlap.bed",
        bedSNO = "data/{ref}.lift.{species}.{species}.nonOverlap.bed"
    params:
        speciesA = config["speciesA"],
        currentSp = correctWildcard
    shell:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.speciesA}/liftOver/{params.speciesA}To{params.currentSp}.over.chain.gz
        python3 scripts/convertToBed.py {params.speciesA}To{params.currentSp}.over.chain.gz tmp.{params.speciesA}.bed tmp.{params.currentSp}.bed
        
        echo Checking overlaps for ref
        bedtools intersect -c -a tmp.{params.speciesA}.bed -b tmp.{params.speciesA}.bed > tmp.overlaps.bed
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.speciesA}.bed tmp.{params.speciesA}.NOH.bed tmp.{params.speciesA}.OH.bed
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.currentSp}.bed tmp.{params.currentSp}.NOH.bed tmp.{params.currentSp}.OH.bed
        
        echo Checking overlaps for species
        bedtools intersect -c -a tmp.{params.currentSp}.NOH.bed -b tmp.{params.currentSp}.NOH.bed > tmp.overlaps.bed
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.speciesA}.NOH.bed {output.bedHNO} {output.bedHO}
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.currentSp}.NOH.bed {output.bedSNO} {output.bedSO}
        """

rule sort:
    input:
        "data/{ref}.lift.{species}.{ref}.nonOverlap.bed"
    output:
        "data/{ref}Lift{species}.sorted.bed"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule intersect:
    input:
        expand("data/{ref}Lift{species}.sorted.bed", ref = config["speciesA"], species = [config["speciesB"],] + config["outgroups"])
    output:
        "data/intersected.bed"
    shadow: "shallow"
    shell:
        """
        touch {output}
        fileList="{input}"
        read -ra fileArray<<<$fileList
        cat ${{fileArray[0]}} > {output}
        for bed in {input}
        do
        export tempFile=$(date +%N)'_intersectTemp.bed'
        bedtools intersect -sorted -a {output} -b $bed > $tempFile
        sort -k1,1 -k2,2n $tempFile > {output}
        done
        """

rule updateInterval:
    input:
        bed = "data/intersected.bed",
        originalH = "data/{ref}.lift.{species}.{ref}.nonOverlap.bed",
        originalS = "data/{ref}.lift.{species}.{species}.nonOverlap.bed"
    output:
        "data/common{ref}Lift{species}.bed"
    shell:
        "python3 scripts/updateInterval.py {input.originalH} {input.originalS} {input.bed} {output}"

rule getFastas:
    output:
        "data/{species}.fa"
    wildcard_constraints:
        species="[A-Za-z\d]+"
    shadow: "shallow"
    shell:
        """
        wget --retry-connrefused --waitretry=5 -t 10 'https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.fa.gz' -P data/
        gunzip data/{wildcards.species}.fa.gz
        """

rule getSeqsFromInt:
    input:
        fa = "data/{species}.fa",
        bed = getCorrectBed
    output:
        fasta = "data/commonSeqs_{species}.fa"
    shell:
        """
        bedtools getfasta -s -fi {input.fa} -bed {input.bed} > {output.fasta}
        """

rule checkLength:
    input:
        fasta = expand("data/commonSeqs_{species}.fa", species = [config["speciesB"],] + config["outgroups"] + [config["speciesA"],])
    output: touch(".checkCompleted")
    shell:
        "python3 scripts/checkLength.py {input}"

rule getAncestralState:
    input:
        fasta = expand("data/commonSeqs_{outgroups}.fa", outgroups = config["outgroups"]),
        ref = "data/commonSeqs_{species}.fa",
        ref2 = getCorrectFasta,
        check = ".checkCompleted"
    output:
        "data/{species}_ancestralBases.vcf"
    shell:
        "./bin/getAncestralBase.bin {input.ref} {input.ref2} {input.fasta} {output}"

rule getAncestralGenome:
    input:
        vcf = "data/{ref}_ancestralBases.vcf",
        bed = "data/common" + config["speciesA"] + "Lift{ref}.bed",
        fa = "data/{ref}.fa"
    output:
        "data/{ref}_ancestralGenome.fasta"
    shell:
        "./bin/makeAncestralGenome.bin {input} {output}"

rule getAncestralGenomeRef:
    input:
        vcf = "data/" + config["speciesA"] + "_ancestralBases.vcf",
        bed = "data/intersected.bed",
        fa = "data/" + config["speciesA"] + ".fa"
    output:
        "data/" + config["speciesA"] + "_ancestralGenome.fasta"
    shell:
        "./bin/makeAncestralGenome.bin {input} {output}"

rule getInt:
    input:
        "data/{species}_ancestralGenome.fasta",
    output:
        "data/{species}_{type}_ints.bed",
        "data/{species}_n{type}_ints.bed"
    params:
        pattern = "{type}"
    shell:
        "./bin/getPattern.bin {params} {input} {output}"

rule filterBarriers:
    input:
        getBarrierFile
    output:
        "data/barriersAOE_{ref}.tsv"
    shell:
        "Rscript scripts/filterInterNIEBs.R {input} {output}"
        
rule filterVCF:
    input:
        "data/{species}_ancestralBases.vcf"
    output:
        "data/{species}.filtered_ancestralBases.vcf"
    shell:
        'grep -v -P "\tN\t" {input} > {output}'

rule countMuts:
    input:
        "data/barriersAOE_{species}.tsv",
        "data/{species}_{type}_ints.bed",
        "data/{species}.filtered_ancestralBases.vcf"
    output:
        "data/{species}_{type}_muts.tsv"
    resources:
        mem_cons = 50
    shell:
        "./bin/countMuts.bin {input} {output}"

rule countBases:
    input:
        "data/{species}.fa",
        "data/barriersAOE_{species}.tsv",
        "data/{species}_{type}_ints.bed"
    output:
        "data/{species}_{type}_bases.tsv"
    shell:
        "./bin/countBases.bin {input} {output}"

rule archivate:
    input:
        "data/{species}_{type}_bases.tsv",
        "data/{species}_{type}_muts.tsv"
    conda:
        "envs/R.yaml"
    output:
        directory("data/" + datetime.datetime.now().strftime("%Y%m%d%H%M") + "{species}_{type}_formatted/"),
        touch(".archived_{species}_{type}")
    shell:
        "Rscript scripts/archiveData.R {input} {output}"

rule prettyFigures:
    input:
        "data/{species}_CG_bases.tsv",
        "data/{species}_nCG_bases.tsv",
        "data/{species}_CG_muts.tsv",
        "data/{species}_nCG_muts.tsv"
    conda:
        "envs/R.yaml"
    output:
        directory("data/{species}_figures"),
        "data/{species}_figures/{species}.rda"
    shell:
        "Rscript scripts/prettyConformFigures.R {input} {output}"

rule commonFigs:
    input:
        "data/{species1}_figures/{species1}.rda",
        "data/{species2}_figures/{species2}.rda"
    output:
        directory("data/{species1}_{species2}_figs/")
    conda:
        "envs/R.yaml"
    shell:
        "Rscript scripts/commonFigs.R {input} {output}"
