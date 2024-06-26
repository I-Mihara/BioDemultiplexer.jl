"""
Preprocesses the barcode file by modifying sequences based on specific criteria.
"""
function preprocess_bc_file(bc_file_path::String, rev::Bool = true)
	bc_df = CSV.read(bc_file_path, DataFrame, delim = "\t")
	for i in 1:nrow(bc_df)
		prefix_region = 1:findfirst('B', bc_df.Full_annotation[i])-1
		suffix_region = findlast('B', bc_df.Full_annotation[i])+1:length(bc_df.Full_annotation[i])
		prefix = SubString(bc_df.Full_seq[i], prefix_region)
		suffix = SubString(bc_df.Full_seq[i], suffix_region)
		bc_df.Full_seq[i] = replace(bc_df.Full_seq[i], prefix => "")
		bc_df.Full_seq[i] = replace(bc_df.Full_seq[i], suffix => "")
	end
	bc_df.Full_seq = uppercase.(bc_df.Full_seq)
	bc_df.Full_seq = replace.(bc_df.Full_seq, "U" => "T")
	if rev == true
		bc_df.Full_seq = reverse.(bc_df.Full_seq)
		bc_df.Full_seq = replace.(bc_df.Full_seq, "A" => "T", "T" => "A", "G" => "C", "C" => "G")
	end
	return bc_df
end

"""
Divides a pair of FASTQ files into smaller parts for parallel processing.
It calculates the number of reads per worker and uses the split command to divide the files.
"""
function divide_fastq(file_R1::String, file_R2::String, output_dir::String, workers::Int)
	divided_dir = joinpath(output_dir, "divided_fastq")
	mkdir(divided_dir)
	total_lines = countlines(file_R1)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	run(`split -l $lines_per_worker -a 5 -d $file_R1 $divided_dir/R1_ --additional-suffix=.fastq`)
	run(`split -l $lines_per_worker -a 5 -d $file_R2 $divided_dir/R2_ --additional-suffix=.fastq`)
end

"""
Divides a single FASTQ file for parallel processing.
"""
function divide_fastq(file_R1::String, output_dir::String, workers::Int)
	divided_dir = joinpath(output_dir, "divided_fastq")
	mkdir(divided_dir)
	total_lines = countlines(file_R1)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	run(`split -l $lines_per_worker -a 5 -d $file_R1 $divided_dir/R1_ --additional-suffix=.fastq`)
end

function mlt_demltplex(thread_num::Int, bc_df::DataFrame, output_dir::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, classify_both::Bool)
	fastq_R1 = output_dir * "/divided_fastq" * "/R1_" * lpad((thread_num - 1), 5, "0") * ".fastq"
	fastq_R2 = output_dir * "/divided_fastq" * "/R2_" * lpad((thread_num - 1), 5, "0") * ".fastq"
	mkdir(output_dir * "/thread" * string(thread_num))
	output_dir = output_dir * "/thread" * string(thread_num)
	classify_sequences(fastq_R1, fastq_R2, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel, classify_both)
end

function mlt_demltplex(thread_num::Int, bc_df::DataFrame, output_dir::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	fastq_R1 = output_dir * "/divided_fastq" * "/R1_" * lpad((thread_num - 1), 5, "0") * ".fastq"
	mkdir(output_dir * "/thread" * string(thread_num))
	output_dir = output_dir * "/thread" * string(thread_num)
	classify_sequences(fastq_R1, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel)
end


function merge_fastq_files(paths::Vector, bc_df::DataFrame, output_dir::String)
	paths_unknown = filter(x -> occursin(r"thread.*/unknown.fastq", x), paths)
	if paths_unknown != []
		run(pipeline(`cat $paths_unknown`, stdout = "$output_dir/unknown.fastq"))
	end

	paths_ambiguous_classification = filter(x -> occursin(r"thread.*/ambiguous_classification.fastq", x), paths)
	if paths_ambiguous_classification != []
		run(pipeline(`cat $paths_ambiguous_classification`, stdout = "$output_dir/ambiguous_classification.fastq"))
	end

	for i in 1:nrow(bc_df)
		regex = r"thread.*/" * string(bc_df.ID[i]) * ".fastq"
		paths_matched = filter(x -> occursin(regex, x), paths)
		if paths_matched != []
			paths_matched = joinpath.(output_dir, paths_matched)
			run(pipeline(`cat $paths_matched`, stdout = "$output_dir/$(bc_df.ID[i]).fastq"))
		end
	end
end

"""
Orchestrates the entire demultiplexing process for FASTQ files.
Handles the preprocessing, dividing, demultiplexing, and merging of files.
"""
function execute_demultiplexing(file_R1::String, file_R2::String, bc_file::String, output_dir::String; max_error_rate::Float64 = 0.22, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, classify_both::Bool = false, bc_rev::Bool = true)
	if isdir(output_dir)
		error("Output directory already exists")
	end
	mkdir(output_dir)
	workers = nworkers()
	bc_df = preprocess_bc_file(bc_file, bc_rev)
	if workers == 1
		classify_sequences(file_R1, file_R2, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel, classify_both)
	else
		divide_fastq(file_R1, file_R2, output_dir, workers)
		pmap(x -> mlt_demltplex(x, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel, classify_both), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(output_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		if classify_both
			paths_R1 = filter(x -> occursin(r"R1", x), paths)
			paths_R2 = filter(x -> occursin(r"R2", x), paths)
			mkdir(output_dir * "/R1")
			mkdir(output_dir * "/R2")
			merge_fastq_files(paths_R1, bc_df, output_dir * "/R1")
			merge_fastq_files(paths_R2, bc_df, output_dir * "/R2")
		else
			merge_fastq_files(paths, bc_df, output_dir)
		end

		rm(joinpath(output_dir, "divided_fastq"), recursive = true)
		for i in 1:workers
			rm(joinpath(output_dir, "thread" * string(i)), recursive = true)
		end
	end
end

function execute_demultiplexing(file_R1::String, bc_file::String, output_dir::String; max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, bc_rev::Bool = true)
	if isdir(output_dir)
		error("Output directory already exists")
	end
	mkdir(output_dir)

	workers = nworkers()
	bc_df = preprocess_bc_file(bc_file, bc_rev)
	if workers == 1
		classify_sequences(file_R1, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel)
	else
		divide_fastq(file_R1, output_dir, workers)
		pmap(x -> mlt_demltplex(x, bc_df, output_dir, max_error_rate, min_delta, mismatch, indel), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(output_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		merge_fastq_files(paths, bc_df, output_dir)

		rm(joinpath(output_dir, "divided_fastq"), recursive = true)
		for i in 1:workers
			rm(joinpath(output_dir, "thread" * string(i)), recursive = true)
		end
	end
end
