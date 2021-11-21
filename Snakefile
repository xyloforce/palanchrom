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
        return "data/intersected.bed"
    else:
        return "data/common" + config["speciesA"] + "Lift{species}.bed"

def correctWildcard(wildcards):
    return wildcards.species[0].upper() + wildcards.species[1:]

# rule getDatFiles:
#     shadow: "shallow"
#     output:
#         path = directory("data/{species}"),
#         interesting_stuff = directory("data/{species}/Non_Overlapping_regions/")
#     params:
#         speciesA = config["speciesA"],
#         currentSp = correctWildcard
#     shell:
#         """
#         wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.speciesA}/liftOver/{params.speciesA}To{params.currentSp}.over.chain.gz
#         python3 scripts/getDatFiles.py -i {params.speciesA}To{params.currentSp}.over.chain.gz -o {output.path}
#         """

# rule convert:
#     input:
#         folder = "data/{species}/Non_Overlapping_regions/"
#     output:
#         "data/" + config["speciesA"] + "Lift{species}.bed"
#     shell:
#         "python3 scripts/datToBed.py {input} {output}"

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
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.speciesA}.NOH.bed {output.bedHO} {output.bedHNO}
        Rscript scripts/sortDuplicate.R tmp.overlaps.bed tmp.{params.currentSp}.NOH.bed {output.bedSO} {output.bedSNO}
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
        check = ".checkCompleted"
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

rule getSortCPG:
    input:
        "data/{species}_ancestralGenome.fasta",
    output:
        cpg = "data/{species}_Ancestral_CPG.fa",
        ncpg = "data/{species}_Ancestral_CPG.fa"
    shell:
        "./bin/sortCPG.bin {input} {output.cpg} {output.ncpg}"

rule filterBarriers:
    input:
        config["barriers"]
    output:
        "data/barriersAOE.tsv"
    shell:
        "Rscript scripts/filterInterNIEBs.R {input} {output}"
        
rule countMuts:
    input:
        "data/barriersAOE.tsv",
        "data/ancestralBases{type}.diff"
    output:
        "data/countBases{type}.tsv"
    shell:
        "./bin/countMuts.bin {input} {output}"

#rule intersectChimpBarriers:
    #input:
        #"filteredBarriersFile",
        #"commonLiftSpeciesChimp" # ATTENTION need to allow for a change in reference
    #output:
        #"barriersInCommonLifts"
