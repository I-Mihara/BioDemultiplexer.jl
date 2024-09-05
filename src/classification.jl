"""
This function aligns `query` and `ref` strings, using semiglobal alignment algorithm. 
# Returns
A similarity score as a float, where higher values indicate better alignment.(0<=similarity_score<=1)
"""
# function semiglobal_alignment(query::String, ref::String; m::Int = Base.length(query), n::Int = Base.length(ref), match::Int = 0, mismatch::Int = -1, indel::Int = -1)
# 	if m == 0 || n == 0
# 		return 0
# 	end
# 	DP = [indel * i for i in 1:m] # Initialize the DP vector.
# 	# Run DP column by column.
# 	for j in 1:n
# 		previous_score = 0
# 		for i in 1:m
# 			insertion_score = DP[i] + (i == m ? 0 : indel)#→
# 			deletion_score = previous_score + indel#↓
# 			substitution_score = (i == 1 ? 0 : DP[i-1]) + (query[i] == ref[j] ? match : mismatch)#↘︎
# 			if i != 1
# 				DP[i-1] = previous_score
# 			end
# 			previous_score = max(insertion_score, deletion_score, substitution_score)
# 		end
# 		DP[m] = previous_score
# 	end

# 	return 1 + DP[m] / m
# end

function semiglobal_alignment(query::String, ref::String, max_error::Float64; match::Int = 0, mismatch::Int = 1, indel::Int = 1)
	m = Base.length(query)
	n = Base.length(ref)
	if m == 0 || n == 0
		return Inf
	end
	score = Inf
	allowed_error = floor(max_error * m) |> Int
	DP = [indel * i for i in 1:m] # Initialize the DP vector.
	# Run DP column by column.
	lact = min(allowed_error + 1, m)
	for j in 1:n
		if m - div(allowed_error, indel) - n + j >= 1#
			fact = m - div(allowed_error, indel) - n + j#
			previous_score = allowed_error
		else
			fact = 1
			previous_score = 0
		end
		if fact > lact
			return score / m
		end
		for i in fact:lact
			insertion_score = (i == m ? Inf : DP[i]+indel)#→#→
			deletion_score = previous_score + indel#↓
			substitution_score = (i == 1 ? 0 : DP[i-1]) + (query[i] == ref[j] ? match : mismatch)#↘︎
			if i != 1
				DP[i-1] = previous_score
			end
			previous_score = min(insertion_score, deletion_score, substitution_score)
		end
		DP[lact] = previous_score
		while lact > 0 && DP[lact] > allowed_error
			lact -= 1
		end
		if lact == m
			score=min(score,previous_score)
		else
			lact += 1
		end
	end
	return score / m
end


# function semiglobal_alignment(query::String, ref::String, max_error::Float64; match::Int = 0, mismatch::Int = 1, indel::Int = 1)
# 	m = Base.length(query)
# 	n = Base.length(ref)
# 	if m == 0 || n == 0
# 		return -Inf
# 	end
# 	allowed_error = floor(max_error * m) |> Int
# 	DP = [indel * i for i in 1:m] # Initialize the DP vector.
# 	score = Inf
# 	# Run DP column by column.
# 	fact = 1
# 	lact = allowed_error
# 	for j in 1:n
# 		lact = min(lact, m - 1)
# 		if m-div(allowed_error,indel) - n + j >= 2
# 		# if indel * (m-n+j-1) > allowed_error#m-div(allowed_error,indel) - n + j >= 2
# 			fact = max(fact,m - div(allowed_error, indel) - n + j)
# 			previous_score = allowed_error
# 			while (indel*(m-n+j-fact)+ DP[fact-1]) > max_error
# 				fact+=1
# 				if fact > lact + 1
# 					break
# 				end
# 			end
# 		elseif m-div(allowed_error,indel) - n + j == 1
# 			fact = 1
# 			previous_score = allowed_error
# 			if indel*(m-n+j-fact) > max_error
# 				fact+=1
# 				while (indel*(m-n+j-fact)+ DP[fact-1]) > max_error
# 					fact+=1
# 					if fact > lact + 1
# 						break
# 					end
# 				end
# 			end
# 		else
# 			fact = 1
# 			previous_score = 0
# 		end
# 		if fact > lact + 1
# 			return -Inf
# 		end
# 		for i in fact:lact+1
# 			insertion_score = (i == m ? Inf : DP[i] + indel)#→
# 			deletion_score = previous_score + indel#↓
# 			substitution_score = (i == 1 ? 0 : DP[i-1]) + (query[i] == ref[j] ? match : mismatch)#↘︎
# 			if i != 1
# 				DP[i-1] = previous_score
# 			end
# 			previous_score = min(insertion_score, deletion_score, substitution_score)
# 		end
# 		if lact == m-1
# 			score=min(score,previous_score)
# 		else
# 			DP[lact+1] = previous_score
# 		end
# 		while lact > -1 && DP[lact+1] > allowed_error
# 			lact -= 1
# 		end
# 		lact += 1
# 	end
# 	return 1 - score / m
# end

"""
#変えた。
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(max_score_bc, delta)`, where `max_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the highest and second-highest scores.
"""
function find_best_matching_bc(seq::String, bc_df::DataFrame, max_error_rate::Float64, mismatch::Int, indel::Int)
	min_score = Inf
	sub_min_score = Inf
	min_score_bc = 0

	for (i, row) in enumerate(eachrow(bc_df))
		alignment_score = semiglobal_alignment(row.Full_seq, seq, max_error_rate, mismatch = mismatch, indel = indel)

		if alignment_score <= max_error_rate
			if alignment_score < min_score
				sub_min_score = min_score
				min_score = alignment_score
				min_score_bc = i
			elseif alignment_score < sub_min_score
				sub_min_score = alignment_score
			end
		end
	end

	delta = sub_min_score - min_score
	return min_score_bc, delta
end

function determine_filename(seq::String, bc_df::DataFrame, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	min_score_bc, delta = find_best_matching_bc(seq, bc_df, max_error_rate, mismatch, indel)

	if min_score_bc == 0
		return "/unknown.fastq"
	elseif delta < min_delta
		return "/ambiguous_classification.fastq"
	else
		return "/" * string(bc_df.ID[min_score_bc]) * ".fastq"
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
function classify_sequences(fastq_R1::String, fastq_R2::String, bc_df::DataFrame, output_dir::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int = 1, indel::Int = 1, classify_both = false)
	if classify_both
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
	else
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
						filename = determine_filename(line1, bc_df, max_error_rate, min_delta)
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
	end
end

function classify_sequences(fastq_R1::String, bc_df::DataFrame, output_dir::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int = 1, indel::Int = 1)
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
				filename = determine_filename(seq, bc_df, max_error_rate, min_delta, mismatch, indel)
				write_fastq_entry(output_dir * filename, header, seq, plus, quality_score)
				mode = "header"
			end
		end
	end
end
