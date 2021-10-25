configfile: "config.yaml"
ruleorder: getSeqsRef > getSeqsFromInt

rule all:
    input:
        "data/{ref}ancestralBases.vcf"

def getSpeciesFolder(wildcards):
    return config["liftover"][wildcards.species]

rule convert:
    input:
        folder = getSpeciesFolder
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
        folder = getSpeciesFolder
    output:
        "data/common" + config["speciesA"] + "Lift{species}.bed"
    shell:
        "python3 scripts/updateInterval.py {input.folder} {input.bed} {output}"

rule getFastas:
    output:
        "data/{species}.fa"
    wildcard_constraints:
        species="[A-Za-z\d]+"
    shell:
        """
        wget --retry-connrefused --waitretry=5 -t 10 'https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.fa.gz' -P data/
        gunzip data/{wildcards.species}.fa.gz
        """

rule getSeqsRef:
    input:
        fa = "data/" + config["speciesA"] + ".fa",
        bed = "data/intersected.bed"
    output:
        "data/commonSeqs_" + config["speciesA"] + ".fa"
    shell:
        """
        bedtools getfasta -s -fi {input.fa} -bed {input.bed} > {output}
        """

rule getSeqsFromInt:
    input:
        fa = "data/{species}.fa",
        bed = "data/common" + config["speciesA"] + "Lift{species}.bed"
    output:
        fasta = "data/commonSeqs_{species}.fa"
    #wildcard_constraints:
        #species="^[8]"
    shell:
        """
        bedtools getfasta -s -fi {input.fa} -bed {input.bed} > {output.fasta}
        """
        
rule checkLength:
    input:
        fasta = expand("data/commonSeqs_{species}.fa", species = [config["speciesB"],] + config["outgroups"]),
        ref = rules.getSeqsRef.output
    output: touch(".checkCompleted")
    shell:
        "python3 scripts/checkLength.py {input}"

rule getAncestralState:
# will output file with headers and pos and seq set for the ref2 !!!
    input:
        fasta = expand("data/commonSeqs_{species}.fa", species = config["outgroups"]),
        ref = "data/commonSeqs_{ref}.fa",
        check = ".checkCompleted"
    output:
        "data/{ref}ancestralBases.vcf"
    shell:
        "./bin/getAncestralBase.bin {input.ref} {input.fasta} {output}"

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
