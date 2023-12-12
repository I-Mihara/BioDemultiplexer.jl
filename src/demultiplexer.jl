@everywhere begin
   """
   semiglobal_alignment is a function that aligns two strings.
   The first string is a query string and the second string is a sequence string.
   The function returns the score of the alignment and the error rate of the alignment.
   """
   #search a string from b string
   function semiglobal_alignment(query::String, seq::String; m::Int=Base.length(query), n::Int=Base.length(seq), match::Int=0, mismatch::Int=-1, indel::Int=-1)
      if m == 0 || n == 0
         return 0
      end
      DP = [indel * i for i in 1:m]
      # run dynamic programming column by column (except the last column)
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
      return 1+DP[m]/m
   end
   
   """
   classify_sequence is a function that classifies a query string into a sequence string.
   The function returns the number of the sequence string and whether the query string is assigned to multiple sequence strings.
   """
   function classify_sequence(seq::String, bc_df::DataFrame,match::Int=0,mismatch::Int=-1,indel::Int=-2,max_error_rate::Float64=0.22)
      maximum_score = -Inf
      maximum_score_number = 0
      multi_assignment = false
      lenseq=length(seq)
      for (i, row) in enumerate(eachrow(bc_df))
         score, errorrate = semiglobal_alignment(row.Full_seq, seq, match=match, mismatch=mismatch, indel=indel, m=length(row.Full_seq), n=lenseq)
         if errorrate < max_error_rate
            if maximum_score < score
               multi_assignment = false
               maximum_score = score
               maximum_score_number = i
            elseif maximum_score == score
               multi_assignment = true
            end
         end
      end
      return maximum_score_number, multi_assignment
   end
   """
   determine_filemname is a function that determines the file name of the output file.
   """
   function determine_filemname(seq::String, bc_df::DataFrame; match::Int=0, mismatch::Int=-1, indel=-2, max_error_rate=0.22)
      class_barcode, multi_assignment = classify_sequence(seq, bc_df, match, mismatch, indel, max_error_rate)
      if multi_assignment
         return "/multi_assignment.fastq"
      elseif class_barcode == 0
         return "/unknown.fastq"
      else
         return "/" * string(bc_df.ID[class_barcode]) * ".fastq"
      end
   end
      

   """
   write_fastq_entry is a function that writes a fastq entry to a file.
   """  
   function write_fastq_entry(filepath, header, seq, plus, quality)
      open(filepath, "a") do outputfile
         write(outputfile, header * "\n" * seq * "\n" * plus * "\n" * quality * "\n")
      end
   end

   """
   process_fastq_file is a function that processes a fastq file.
   """
   function process_fastq_file(fastq_R1::String, fastq_R2::String, bc_df::DataFrame, output_dir::String; match=0, mismatch=-1, indel=-2, classify_r1=false, classify_r2=true)
      if classify_r2
         if classify_r1
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
                        header = line
                        header2 = line2
                        mode = "seq"
                     elseif mode == "seq"
                        seq = line
                        seq2 = line2
                        mode = "plus"
                     elseif mode == "plus"
                        plus = line
                        plus2 = line2
                        mode = "quality_score"
                     elseif mode == "quality_score"
                        quality_score = line
                        quality_score2 = line2
                        filename = determine_filemname(seq, bc_df, match=match, mismatch=mismatch, indel=indel)
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
                  count = 0

                  for line1 in eachline(primary_file)
                     line2 = readline(secondary_file)
                        if line1[1] == '@' && mode == "header"
                           header2 = line2
                           mode = "seq"
                        elseif mode == "seq"
                           filename = determine_filemname(line1, bc_df, match=match, mismatch=mismatch, indel=indel)
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
                     count += 1
                  end
               end
            end
         end
      elseif classify_r1
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
                  filename = determine_filemname(seq, bc_df, match=match, mismatch=mismatch, indel=indel)
                  write_fastq_entry(output_dir * filename, header, seq, plus, quality_score)
                  mode = "header"
               end
            end
         end
      end
   end
end