add_library(bed_tools bed_tools.cpp vcf_tools.cpp)
add_library(fasta_tools fasta_tools.cpp bio_tools.cpp bed_tools.cpp vcf_tools.cpp)
add_library(bio_tools bio_tools.cpp fasta_tools.cpp bed_tools.cpp vcf_tools.cpp)
add_library(vcf_tools vcf_tools.cpp bed_tools.cpp)

add_executable(makeAncestral makeAncestralGenome.cpp)
add_executable(getAncestralBase getAncestralBase.cpp)
add_executable(getPattern getPattern.cpp)
add_executable(countMuts countMuts.cpp)
add_executable(countBases countBases.cpp)

target_link_libraries(makeAncestral PUBLIC bed_tools fasta_tools vcf_tools)
target_link_libraries(getAncestralBase PUBLIC bed_tools fasta_tools vcf_tools)
target_link_libraries(getPattern PUBLIC bed_tools fasta_tools bio_tools)
target_link_libraries(countMuts PUBLIC bed_tools fasta_tools bio_tools)
target_link_libraries(countBases PUBLIC bed_tools fasta_tools bio_tools)

install(TARGETS makeAncestral getAncestralBase getPattern countMuts countBases RUNTIME DESTINATION ~/Documents/Projets/Reprise_Sara/bin/)
