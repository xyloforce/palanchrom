import datetime
configfile: "config.yaml"
wildcard_constraints:
    species = "\w+"

rule all:
    input:
        "." + datetime.datetime.now().strftime("%Y_%m_%d") + ".panTro5.hg38.CG.nCG"

def getLiftoverFile(wildcards):
    return config["liftover"][wildcards.species]

def getCorrectBed(wildcards):
    if wildcards.species == config["speciesA"]:
        return config["result_folder"] + "/intersected.bed"
    else:
        return config["result_folder"] + "/common" + config["speciesA"] + "Lift{species}.bed"

def getCorrectFasta(wildcards):
    if wildcards.species == config["speciesA"]:
        return config["result_folder"] + "/commonSeqs_" + config["speciesB"] + ".fa"
    else:
        return config["result_folder"] + "/commonSeqs_" + config["speciesA"] + ".fa"

def correctWildcard(wildcards):
    return wildcards.species[0].upper() + wildcards.species[1:]

def getBarrierFile(wildcards):
    return config["barriers"][wildcards.ref]

rule getNonOverlapInt:
    shadow: "shallow"
    output:
        bedHO = config["result_folder"] + "/{ref}.lift.{species}.{ref}.overlap.bed",
        bedHNO = config["result_folder"] + "/{ref}.lift.{species}.{ref}.nonOverlap.bed",
        bedSO = config["result_folder"] + "/{ref}.lift.{species}.{species}.overlap.bed",
        bedSNO = config["result_folder"] + "/{ref}.lift.{species}.{species}.nonOverlap.bed"
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
        config["result_folder"] + "/{ref}.lift.{species}.{ref}.nonOverlap.bed"
    output:
        config["result_folder"] + "/{ref}Lift{species}.sorted.bed"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule intersect:
    input:
        expand(config["result_folder"] + "/{ref}Lift{species}.sorted.bed", ref = config["speciesA"], species = [config["speciesB"],] + config["outgroups"])
    output:
        config["result_folder"] + "/intersected.bed"
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
        bed = config["result_folder"] + "/intersected.bed",
        originalH = config["result_folder"] + "/{ref}.lift.{species}.{ref}.nonOverlap.bed",
        originalS = config["result_folder"] + "/{ref}.lift.{species}.{species}.nonOverlap.bed"
    output:
        config["result_folder"] + "/common{ref}Lift{species}.bed"
    shell:
        "python3 scripts/updateInterval.py {input.originalH} {input.originalS} {input.bed} {output}"

rule getFastas:
    output:
        config["result_folder"] + "/{species}.fa"
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
        fa = config["result_folder"] + "/{species}.fa",
        bed = getCorrectBed
    output:
        fasta = config["result_folder"] + "/commonSeqs_{species}.fa"
    shell:
        """
        bedtools getfasta -s -fi {input.fa} -bed {input.bed} > {output.fasta}
        """

rule checkLength:
    input:
        fasta = expand(config["result_folder"] + "/commonSeqs_{species}.fa", species = [config["speciesB"],] + config["outgroups"] + [config["speciesA"],])
    output: touch(".checkCompleted")
    shell:
        "python3 scripts/checkLength.py {input}"

rule getAncestralState:
    input:
        fasta = expand(config["result_folder"] + "/commonSeqs_{outgroups}.fa", outgroups = config["outgroups"]),
        ref = config["result_folder"] + "/commonSeqs_{species}.fa",
        ref2 = getCorrectFasta,
        check = ".checkCompleted"
    output:
        config["result_folder"] + "/{species}_ancestralBases.vcf"
    shell:
        "./bin/getAncestralBase {input.ref} {input.ref2} {input.fasta} {output}"

rule getAncestralGenome:
    input:
        vcf = config["result_folder"] + "/{ref}_ancestralBases.vcf",
        bed = config["result_folder"] + "/common" + config["speciesA"] + "Lift{ref}.bed",
        fa = config["result_folder"] + "/{ref}.fa"
    output:
        config["result_folder"] + "/{ref}_ancestralGenome.fasta"
    shell:
        "./bin/makeAncestralGenome {input} {output}"

rule getAncestralGenomeRef:
    input:
        vcf = config["result_folder"] + "/" + config["speciesA"] + "_ancestralBases.vcf",
        bed = config["result_folder"] + "/intersected.bed",
        fa = config["result_folder"] + "/" + config["speciesA"] + ".fa"
    output:
        config["result_folder"] + "/" + config["speciesA"] + "_ancestralGenome.fasta"
    shell:
        "./bin/makeAncestralGenome {input} {output}"

rule getInt:
    input:
        config["result_folder"] + "/{species}_ancestralGenome.fasta",
    output:
        config["result_folder"] + "/{species}_{type}_ints.bed",
        config["result_folder"] + "/{species}_n{type}_ints.bed"
    params:
        pattern = "{type}"
    shell:
        "./bin/getPattern {params} {input} {output}"

rule filterBarriers:
    input:
        getBarrierFile
    output:
        config["result_folder"] + "/barriersAOE_{ref}.tsv"
    shell:
        "Rscript scripts/filterInterNIEBs.R {input} {output}"
        
rule filterVCF:
    input:
        config["result_folder"] + "/{species}_ancestralBases.vcf"
    output:
        config["result_folder"] + "/{species}.filtered_ancestralBases.vcf"
    shell:
        'grep -v -P "\tN\t" {input} > {output}'

rule countMuts:
    input:
        config["result_folder"] + "/barriersAOE_{species}.tsv",
        config["result_folder"] + "/{species}_{type}_ints.bed",
        config["result_folder"] + "/{species}.filtered_ancestralBases.vcf"
    output:
        config["result_folder"] + "/{species}_{type}_muts.tsv"
    shadow: "shallow"
    resources:
        mem_cons = 25
    shell:
        "./bin/countMuts {input} {output}"

rule countBases:
    input:
        config["result_folder"] + "/{species}_ancestralGenome.fasta",
        config["result_folder"] + "/barriersAOE_{species}.tsv",
        config["result_folder"] + "/{species}_{type}_ints.bed"
    output:
        config["result_folder"] + "/{species}_{type}_bases.tsv"
    shadow: "shallow"
    resources:
        mem_cons = 50
    shell:
        "./bin/countBases {input} {output}"

rule archivate:
    input:
        config["result_folder"] + "/{species}_{type}_bases.tsv",
        config["result_folder"] + "/{species}_{type}_muts.tsv"
    conda:
        "envs/R.yaml"
    output:
        folder = directory(config["result_folder"] + "/{datetime}_{species}_{type}_formatted/"),
        touch = touch(".archived_{datetime}_{species}_{type}")
    shell:
        "Rscript scripts/archiveData.R {input} {output.folder}"

rule setupPackages:
    output:
        touch(".R.setup")
    conda:
        "envs/R.yaml"
    shell:
        "Rscript scripts/setupPackages.R"

rule prettyFigures:
    input:
        folder = config["result_folder"] + "/{datetime}_{species}_{type}_formatted/",
        setup = ".R.setup"
    conda:
        "envs/R.yaml"
    output:
        folder = directory(config["result_folder"] + "/{datetime}_{species}_{type}_figures"),
        touch = touch(".figures_{datetime}_{species}_{type}")
    shell:
        "Rscript scripts/prettyConformFigures.R {input.folder} {output.folder}"

rule commonFigs:
    input:
        folder1 = config["result_folder"] + "/{datetime}_{species1}_{type}_formatted/",
        folder2 = config["result_folder"] + "/{datetime}_{species2}_{type}_formatted/",
        setup = ".R.setup"
    output:
        directory(config["result_folder"] + "/{datetime}_{species1}_{species2}_{type}_figs/")
    conda:
        "envs/R.yaml"
    shell:
        "Rscript scripts/commonFigs.R {input.folder1} {input.folder2} {output}"

rule sendToPhystorage:
    input:
        config["result_folder"] + "/{datetime}_{species1}_{type1}_formatted/",
        config["result_folder"] + "/{datetime}_{species1}_{type2}_formatted/",
        config["result_folder"] + "/{datetime}_{species2}_{type1}_formatted/",
        config["result_folder"] + "/{datetime}_{species2}_{type2}_formatted/",
        config["result_folder"] + "/{datetime}_{species1}_{species2}_{type1}_figs/",
        config["result_folder"] + "/{datetime}_{species1}_{species2}_{type2}_figs/",
        config["result_folder"] + "/{datetime}_{species1}_{type1}_figures",
        config["result_folder"] + "/{datetime}_{species1}_{type2}_figures",
        config["result_folder"] + "/{datetime}_{species2}_{type1}_figures",
        config["result_folder"] + "/{datetime}_{species2}_{type2}_figures"
    output:
        touch(".{datetime}.{species1}.{species2}.{type1}.{type2}")
    params:
        folder = config["result_folder"],
        datetime = "{datetime}"
    shell:
        """
        FOLDER={params.folder}/{params.datetime}_results
        if [ -d "$FOLDER" ]; then
        rm -r $FOLDER
        fi
        mkdir $FOLDER
        mv {input} {params.folder}/{params.datetime}_results
        rsync -av --exclude='.*' {params.folder}/{params.datetime}_results fsassola@phystorage.physique.ens-lyon.fr:/partages/Bioinfo/shared/users/fsassola
        """
