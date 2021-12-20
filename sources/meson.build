bed_tools = static_library('bed_tools',
		'bed_tools.cpp',
		'vcf_tools.cpp')

fasta_tools = static_library('fasta_tools',
		'fasta_tools.cpp',
		'bio_tools.cpp',
		'bed_tools.cpp',
		'vcf_tools.cpp')

bio_tools = static_library('bio_tools',
		'bio_tools.cpp',
		'fasta_tools.cpp',
		'bed_tools.cpp',
		'vcf_tools.cpp')

vcf_tools = static_library('vcf_tools',
		'vcf_tools.cpp',
		'bed_tools.cpp')

makeAncestral_sources = ['makeAncestralGenome.cpp']

executable('makeAncestralGenome.bin',
	makeAncestral_sources,
	install: true,
	link_with : [bed_tools, fasta_tools, vcf_tools])

getAncestral_sources = ['getAncestralBase.cpp']

executable('getAncestralBase.bin',
	getAncestral_sources,
	install: true,
	link_with : [bed_tools, fasta_tools, vcf_tools])

getCPGInt_sources = ['getCPGInt.cpp']

executable('getCPGInt.bin',
	getCPGInt_sources,
	install: true,
	link_with : [bed_tools, fasta_tools, bio_tools])

countsMuts_sources = ['countMuts.cpp']

executable('countMuts.bin',
	countsMuts_sources,
	install: true,
	link_with : [bed_tools, bio_tools, vcf_tools])

countBases_sources = ['countBases.cpp']

executable('countBases.bin',
	countBases_sources,
	install: true,
	link_with: [bed_tools, fasta_tools, bio_tools])

getNonPattern_sources = ['getNonPattern.cpp']

executable('testNonPattern.bin',
	getNonPattern_sources,
	install: true,
	link_with: [bio_tools])
