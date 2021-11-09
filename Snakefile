configfile: "config.yaml"
wildcard_constraints:
    species = "\w+"

rule all:
    input:
        "data/hg38_ancestralGenome.fasta",
        "data/panTro5_ancestralGenome.fasta"

def getLiftoverFile(wildcards):
    return config["liftover"][wildcards.species]

def getCorrectFile(wildcards):
    if wildcards.species == config["speciesA"]:
        return "intersected.bed"
    else:
        return "data/common" + config["speciesA"] + "Lift{species}.bed"

def returnPath(wildcards):
    return wildcards.species[0].toupper() + wildcards.species[1:]

rule getDatFiles:
    output:
        path = directory("data/{species}"),
        interesting_stuff = directory("data/{species}/Non_Overlapping_regions/")
    params:
        speciesA = config["speciesA"],
        currentSp = "{species}"
    shell:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.speciesA}/liftOver/{params.speciesA}to{params.currentSp}.over.chain.gz
        python3 scripts/getDatFiles.py -i {params.speciesA}to{params.currentSp}.over.chain.gz -o {output.path}
        """

rule convert:
    input:
        folder = "data/{species}/Non_Overlapping_regions/"
    output:
        "data/" + config["speciesA"] + "Lift{species}.bed"
    shell:
        "python3 scripts/datToBed.py {input} {output}"

rule sort:
    input:
        "data/{ref}Lift{species}.bed"
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
        folder = "data/{species}/Non_Overlapping_regions/"
    output:
        "data/common" + config["speciesA"] + "Lift{species}.bed"
    shell:
        "python3 scripts/updateInterval.py {input.folder} {input.bed} {output}"

rule getFastas:
    output:
        "data/{species}.fa",
        "data/{species}.chrom.sizes"
    wildcard_constraints:
        species="[A-Za-z\d]+"
    shell:
        """
        wget --retry-connrefused --waitretry=5 -t 10 'https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.fa.gz' -P data/
        wget --retry-connrefused --waitretry=5 -t 10 'https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.chrom.sizes' -P data/
        gunzip data/{wildcards.species}.fa.gz
        """

rule getSeqsFromInt:
    input:
        fa = "data/{species}.fa",
        bed = getCorrectFile
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
# will output file with headers and pos and seq set for the ref2 !!!
    input:
        fasta = expand("data/commonSeqs_{outgroups}.fa", outgroups = config["outgroups"]),
        ref = "data/commonSeqs_{species}.fa",
        check = ".checkCompleted",
        chromSizes = "data/{species}.chrom.sizes"
    output:
        "data/{species}_ancestralBases.vcf"
    shell:
        "./bin/getAncestralBase.bin {input.ref} {input.fasta} {output}"

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

#rule getSortCPG:
    #input:
        #"data/ancestralBases.fmna"
    #output:
        #cpg = "data/ancestralBasesCPG.diff",
        #ncpg = "data/ancestralBasesNCPG.diff"
    #shell:
        #"./bin/sortCPG.bin {input} {output.cpg} {output.ncpg}"

#rule filterBarriers:
    #input:
        #config["barriers"]
    #output:
        #"data/barriersAOE.tsv"
    #shell:
        #"Rscript scripts/filterInterNIEBs.R {input} {output}"
        
#rule countMuts:
    #input:
        #"data/barriersAOE.tsv",
        #"data/ancestralBases{type}.diff"
    #output:
        #"data/countBases{type}.tsv"
    #shell:
        #"./bin/countMuts.bin {input} {output}"

#rule intersectChimpBarriers:
    #input:
        #"filteredBarriersFile",
        #"commonLiftSpeciesChimp" # ATTENTION need to allow for a change in reference
    #output:
        #"barriersInCommonLifts"
