@everywhere begin
	function preprocess_bc_file(bc_file_path::String; rev = true)
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

	function execute_demultiplexing(file_R1, file_R2, bc_file, output_dir; max_error_rate = 0.2, min_delta = 0.1, classify = "R2", bc_rev = true)
		if isdir(output_dir)
			error("Output directory already exists")
		end
		mkdir(output_dir)
		workers = nworkers()
		bc_df = preprocess_bc_file(bc_file)
		if workers == 1
			classify_seqences(file_R1, file_R2, bc_df, output_dir, max_error_rate, min_delta, classify = classify)
		else
			divide_fastq(file_R1, file_R2, output_dir, workers)
			pmap(x -> mlt_demltplex(x, bc_df, output_dir, max_error_rate, min_delta, classify = classify), 1:workers)
			paths = []
			for (root, dirs, files) in walkdir(output_dir)
				for file in files
					push!(paths, joinpath(root, file))
				end
			end
			paths_unknown = filter(x -> occursin(r"thread.*/unknown.fastq", x), paths)
			paths_ambiguous_classification = filter(x -> occursin(r"thread.*/ambiguous_classification.fastq", x), paths)
			# merge files
			run(pipeline(`cat $paths_unknown`, stdout = "$output_dir/unknown.fastq"))
			if paths_ambiguous_classification != []
				run(pipeline(`cat $paths_ambiguous_classification`, stdout = "$output_dir/ambiguous_classification.fastq"))
			end
			for i in 1:nrow(bc_df)
				regex1 = r"thread.*/" * string(bc_df.ID[i]) * ".fastq"
				paths_matched = filter(x -> occursin(regex1, x), paths)
				if paths_matched != []
					paths_matched = joinpath.(output_dir, paths_matched)
					run(pipeline(`cat $paths_matched`, stdout = "$output_dir/$(bc_df.ID[i]).fastq"))
				end
			end
			rm(joinpath(output_dir, "divided_fastq"), recursive = true)
			for i in 1:workers
				rm(joinpath(output_dir, "thread" * string(i)), recursive = true)
			end
		end
	end

	function divide_fastq(file_R1, file_R2, output_dir, workers)
		div_dir = joinpath(output_dir, "divided_fastq")
		mkdir(div_dir)
		num_lines = countlines(file_R1)
		num_reads = num_lines ÷ 4
		reads_per_thread = cld(num_reads, workers)
		lines_per_thread = reads_per_thread * 4
		run(`split -l $lines_per_thread -a 5 -d $file_R1 $div_dir/R1_ --additional-suffix=.fastq`)
		run(`split -l $lines_per_thread -a 5 -d $file_R2 $div_dir/R2_ --additional-suffix=.fastq`)
	end

	function mlt_demltplex(thread_num, bc_df, output_dir, max_error_rate, min_delta, classify)
		fastq_R1 = output_dir * "/divided_fastq" * "/R1_" * lpad((thread_num - 1), 5, "0") * ".fastq"
		fastq_R2 = output_dir * "/divided_fastq" * "/R2_" * lpad((thread_num - 1), 5, "0") * ".fastq"
		mkdir(output_dir * "/thread" * string(thread_num))
		output_dir = output_dir * "/thread" * string(thread_num)
		classify_seqences(fastq_R1, fastq_R2, bc_df, output_dir, max_error_rate, min_delta, classify = classify)
	end
end
