@everywhere begin
	"""
	This function aligns `query` and `seq` strings, using semiglobal alignment algorithm. 
	# Returns
	A similarity score as a float, where higher values indicate better alignment.(0<=similarity_score<=1)
	"""
	function semiglobal_alignment(query::String, seq::String; m::Int = Base.length(query), n::Int = Base.length(seq), match::Int = 0, mismatch::Int = -1, indel::Int = -1)
		if m == 0 || n == 0
			return 0
		end
		DP = [indel * i for i in 1:m] # Initialize the DP vector.
		# Run DP column by column.
		for j in 1:n
			previous_score = 0
			for i in 1:m
				insertion_score = DP[i] + (i == m ? 0 : indel)#→
				deletion_score = previous_score + indel#↓
				substitution_score = (i == 1 ? 0 : DP[i-1]) + (query[i] == seq[j] ? match : mismatch)#↘︎
				if i != 1
					DP[i-1] = previous_score
				end
				previous_score = max(insertion_score, deletion_score, substitution_score)
			end
			DP[m] = previous_score
		end

		return 1 + DP[m] / m
	end

	"""
	Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
	# Returns
	A tuple `(max_score_bc, delta)`, where `max_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the highest and second-highest scores.
	"""
	function fine_best_matching_seq(seq::String, bc_df::DataFrame, max_error_rate::Float64)
		max_score = -1.0
		sub_max_score = -1.0
		max_score_bc = 0
		lenseq = length(seq)

		for (i, row) in enumerate(eachrow(bc_df))
			similarity_score = semiglobal_alignment(row.Full_seq, seq, n = lenseq)

			if similarity_score > 1.0 - max_error_rate
				if similarity_score > max_score
					sub_max_score = max_score
					max_score = similarity_score
					max_score_bc = i
				elseif similarity_score > sub_max_score
					sub_max_score = similarity_score
				end
			end
		end

		delta = max_score - sub_max_score
		return max_score_bc, delta
	end

	function determine_filename(seq::String, bc_df::DataFrame, max_error_rate::Float64, min_delta::Float64)
		max_score_bc, delta = find_best_matching_seq(seq, bc_df, max_error_rate)

		if max_score_bc == 0
			return "/unknown.fastq"
		elseif delta < min_delta
			return "/ambiguous_classification.fastq"
		else
			return "/" * string(bc_df.ID[max_score_bc]) * ".fastq"
		end
	end

	function write_fastq_entry(filepath, header, seq, plus, quality)
		open(filepath, "a") do outputfile
			write(outputfile, header * "\n" * seq * "\n" * plus * "\n" * quality * "\n")
		end
	end

	"""
	Compare each sequence in the fastq_R1 file with the sequences in bc_df, and classify the sequences of the specified file based on that comparison.
	"""
	function classify_seqences(fastq_R1::String, fastq_R2::String, bc_df::DataFrame, output_dir::String, max_error_rate::Float64, min_delta::Float64; classify = "R2")
		if classify == "both"
			mkdir(output_dir * "/R1")
			mkdir(output_dir * "/R2")
			open(fastq_R1, "r") do primary_file
				open(fastq_R2, "r") do secondary_file
					header, seq, plus, quality_score = "", "", "", ""
					header2, seq2, plus2, quality_score2 = "", "", "", ""
					mode = "header"
					filename = ""
					for line1 in eachline(primary_file)
						line2 = readline(secondary_file)
						if line1[1] == '@' && mode == "header"
							header = line1
							header2 = line2
							mode = "seq"
						elseif mode == "seq"
							seq = line1
							seq2 = line2
							mode = "plus"
						elseif mode == "plus"
							plus = line1
							plus2 = line2
							mode = "quality_score"
						elseif mode == "quality_score"
							quality_score = line1
							quality_score2 = line2
							filename = determine_filename(seq, bc_df, max_error_rate, min_delta)
							write_fastq_entry(output_dir * "/R1" * filename, header, seq, plus, quality_score)
							write_fastq_entry(output_dir * "/R2" * filename, header2, seq2, plus2, quality_score2)
							mode = "header"
						end
					end
				end
			end
		elseif classify == "R2"
			open(fastq_R1, "r") do primary_file
				open(fastq_R2, "r") do secondary_file
					header2, seq2, plus2, quality_score2 = "", "", "", ""
					mode = "header"
					filename = ""
					for line1 in eachline(primary_file)
						line2 = readline(secondary_file)
						if line1[1] == '@' && mode == "header"
							header2 = line2
							mode = "seq"
						elseif mode == "seq"
							filename = determine_filemname(line1, bc_df, max_error_rate, min_delta)
							seq2 = line2
							mode = "plus"
						elseif mode == "plus"
							plus2 = line2
							mode = "quality_score"
						elseif mode == "quality_score"
							quality_score2 = line2
							write_fastq_entry(output_dir * filename, header2, seq2, plus2, quality_score2)
							mode = "header"
						end
					end
				end
			end
		elseif classify == "R1"
			open(fastq_R1, "r") do file
				header, seq, plus, quality_score = "", "", "", ""
				mode = "header"
				filename = ""
				for line in eachline(file)
					if line[1] == '@' && mode == "header"
						header = line
						mode = "seq"
					elseif mode == "seq"
						seq = line
						mode = "plus"
					elseif mode == "plus"
						plus = line
						mode = "quality_score"
					elseif mode == "quality_score"
						quality_score = line
						filename = determine_filemname(seq, bc_df, max_error_rate, min_delta)
						write_fastq_entry(output_dir * filename, header, seq, plus, quality_score)
						mode = "header"
					end
				end
			end
		end
	end
end
