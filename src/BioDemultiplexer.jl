__precompile__()

module BioDemultiplexer

export

	execute_demultiplexing,
	single_demltplex,
	mlt_demltplex,
	divide_fastq,
	preprocess_bc_file,
	process_fastq_file,
	write_fastq_entry,
	determine_filemname,
	classify_sequence,
	semiglobal_alignment


using Distributed
@everywhere using DataFrames, CSV


include("classification.jl")
include("demultiplexing.jl")


end#module
