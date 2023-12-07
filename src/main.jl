@everywhere begin
   """
   preprocess_bc_file is a function to preprocess barcode file.
   This function is used in execute_demultiplexing function.
   """
   function preprocess_bc_file(bc_file_path::String;rev=true)
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
   
   """
   divide_fastq is a function to divide fastq files into the number of workers.
   This function is used in execute_demultiplexing function.
   """
   function divide_fastq(fastq_dir::String, output_dir,workers)
   
      # FASTQファイルのリストを取得
      files = readdir(fastq_dir)
      R1_files = filter(x -> occursin("_R1_", x), files)
      R2_files = filter(x -> occursin("_R2_", x), files)
   
      # ディレクトリを作成
      div_dir = joinpath(output_dir, "divided_fastq")
      if isdir(div_dir)
         rm(div_dir, recursive=true)
      end
      mkdir(div_dir)
   
      #各ファイルを分割する
      for file in R1_files
         num_lines = countlines(joinpath(fastq_dir, file))
         num_reads = num_lines ÷ 4
         reads_per_thread = cld(num_reads, workers)
         open(joinpath(fastq_dir,file), "r") do input
            while !eof(input)
               for thread in 1:workers
                  open(joinpath(div_dir,"thread$(thread)_R1.fastq"),"a") do output
                     for _ in 1:reads_per_thread*4
                        line = readline(input)
                        write(output, line * "\n")
                     end
                  end
               end
            end
         end
      end
   
      #R2も同様に
      for file in R2_files
         num_lines = countlines(joinpath(fastq_dir, file))
         num_reads = num_lines ÷ 4
         reads_per_thread = cld(num_reads, workers)
         open(joinpath(fastq_dir,file), "r") do input
            while !eof(input)
               for thread in 1:workers
                  open(joinpath(div_dir,"thread$(thread)_R2.fastq"),"a") do output
                     for _ in 1:reads_per_thread*4
                        line = readline(input)
                        write(output, line * "\n")
                     end
                  end
               end
            end
         end
      end
   end
   
   
   """
   execute_demultiplexing is a function to execute demultiplexing.
   """
   function execute_demultiplexing(fastq_dir, bc_file, output_dir)
       if isdir(output_dir)
          rm(output_dir, recursive=true)
       end
       mkdir(output_dir)
       workers = nworkers()
       if workers == 1
          single_demltplex(fastq_dir, bc_file, output_dir)
       else
          divide_fastq(fastq_dir, output_dir, workers)
          pmap(x -> mlt_demltplex(x, bc_file, output_dir), 1:workers)
          rm(joinpath(output_dir, "divided_fastq"), recursive=true)
       end
   end
   
   
   """
   mlt_demltplex is a function to demultiplex fastq files.
   """ 
   function mlt_demltplex(thread_num, bc_file, output_dir)
       fastq_R1 = output_dir * "/divided_fastq" * "/thread$(thread_num)_R1.fastq"
       fastq_R2 = output_dir * "/divided_fastq" * "/thread$(thread_num)_R2.fastq"
       bc_df = preprocess_bc_file(bc_file)
       process_fastq_file(fastq_R1, fastq_R2, bc_df, output_dir)   
   end
   
   """
   single_demltplex is a function to demultiplex fastq files.
   """
   function single_demltplex(fastq_dir, bc_file, output_dir)
       files = readdir(fastq_dir)
       R1_files = filter(x -> occursin("_R1_", x), files)
       bc_df = preprocess_bc_file(bc_file)
       for file in R1_files
         # m = match(r"^(.+?)_R1_", filename)
         # mkdir(joinpath(output_dir, m.captures[1]))
         process_fastq_file(joinpath(fastq_dir, file), joinpath(fastq_dir, replace(file, "_R1_" => "_R2_")), bc_df , output_dir)
       end
   end      
end
