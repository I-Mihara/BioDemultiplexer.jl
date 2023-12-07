@everywhere begin
   function preprocess_bc_file(bc_file_path::String; rev=true)
      bc_df = CSV.read(bc_file_path, DataFrame, delim="\t")
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
  
  function execute_demultiplexing(file_R1,file_R2,bc_file,output_dir)
      if isdir(output_dir)
          rm(output_dir, recursive=true)
      end
      mkdir(output_dir)
      workers = nworkers()
      bc_df = preprocess_bc_file(bc_file)
      if workers == 1
          process_fastq_file(file_R1, file_R2, bc_df, output_dir)
      else
           divide_fastq(file_R1,file_R2,output_dir, workers)
           pmap(x -> mlt_demltplex(x, bc_df, output_dir), 1:workers)
           f = []
           for (root, dirs ,files) in walkdir(output_dir)
                 for file in files
                    push!(f, joinpath(root, file))
                 end
           end   
           dr_unknown = filter(x -> occursin(r"thread.*/trimmed-unknown.fastq", x), f)
           dr_multi_assignment = filter(x -> occursin(r"thread.*/trimmed-multi_assignment.fastq", x), f)
           run(pipeline(`cat $dr_unknown` ,stdout= "$output_dir/trimmed-unknown.fastq"))
           if dr_multi_assignment != []
              run(pipeline(`cat $dr_multi_assignment` ,stdout= "$output_dir/trimmed-multi_assignment.fastq"))
           end
           for i in 1:nrow(bc_df)
                 regex1 = r"thread.*/trimmed-"* string(bc_df.ID[i]) * ".fastq"
                 dr_trimmed = filter(x -> occursin(regex1, x), f)
                 if dr_trimmed != []
                    dr_trimmed = joinpath.(output_dir, dr_trimmed)
                    run(pipeline(`cat $dr_trimmed`, stdout= "$output_dir/trimmed-$(bc_df.ID[i]).fastq"))
                 end    
           end
           rm(joinpath(output_dir, "divided_fastq"), recursive=true)
           for i in 1:workers
              rm(joinpath(output_dir, "thread"*string(i)), recursive=true)
           end
      end
  end
      
  function divide_fastq(file_R1,file_R2,output_dir,workers)
      div_dir = joinpath(output_dir, "divided_fastq")
      mkdir(div_dir)
      num_lines = countlines(file_R1)
      num_reads = num_lines ÷ 4
      reads_per_thread = cld(num_reads, workers)
      lines_per_thread = reads_per_thread*4
      run(`split -l $lines_per_thread -a 5 -d $file_R1 $div_dir/R1_ --additional-suffix=.fastq`)
      run(`split -l $lines_per_thread -a 5 -d $file_R2 $div_dir/R2_ --additional-suffix=.fastq`)
  end
  
  function mlt_demltplex(thread_num, bc_df, output_dir)
      fastq_R1 = output_dir * "/divided_fastq" * "/R1_" * lpad((thread_num-1),5,"0") * ".fastq"
      fastq_R2 = output_dir * "/divided_fastq" * "/R2_" * lpad((thread_num-1),5,"0") * ".fastq"
      mkdir(output_dir * "/thread" * string(thread_num))
      output_dir = output_dir * "/thread" * string(thread_num)
      process_fastq_file(fastq_R1, fastq_R2, bc_df, output_dir)
  end      
end
