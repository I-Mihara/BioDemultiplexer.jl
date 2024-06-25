# BioDemultiplexer.jl
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://I-Mihara.github.io/BioDemultiplexer.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://I-Mihara.github.io/BioDemultiplexer.jl/dev/)
[![Build Status](https://github.com/I-Mihara/BioDemultiplexer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/I-Mihara/BioDemultiplexer.jl/actions/workflows/CI.yml?query=branch%3Amain)

BioDemultiplexer is a Julia package designed for demultiplexing reads based on barcodes. It efficiently classifies sequences from FASTQ files, making it a valuable tool in bioinformatics workflows.

## Usage
The primary function of this package is `execute_demultiplexing()`. It classifies sequences in an R2 FASTQ file by calculating similarity scores from R1 sequences and barcodes in a reference file. Usage is as follows:
```Julia 
execute_demultiplexing(file_R1, file_R2, bc_file, output_dir)
```

**Arguments**

* `file_R1::String`: Path to the input R1 FASTQ file.
* `file_R2::String`: Path to the input R2 FASTQ file.
* `bc_file::String`: Path to the reference barcode file in TSV format.
* `output_dir::String`: Path to the directory where demultiplexed files will be saved.

**Optional Arguments**

* `max_error_rate::Float64 = 0.22`: Maximum permissible error rate for sequence assignment.
* `min_delta::Float64 = 0.1`: Minimum difference between the highest and second-highest barcode similarity scores.
* `classify_both::Bool = false`: Set to true to classify both R1 and R2 sequences.
* `bc_rev::Bool = true`: Set to true to reverse the barcode for demultiplexing.

Demultiplexing can also be performed with a single FASTQ file:
```Julia
execute_demultiplexing(file_R1, bc_file, output_dir)
```

## Parallel computing
BioDemultiplexer supports parallel computing, allowing faster processing of large datasets:
```Julia
using Distributed
addprocs(n)# 'n' is the number of desired workers.
@everywhere using BioDemultiplexer
execute_demultiplexing(file_R1, file_R2, bc_file, output_dir)
``` 

## License
This package is licensed under the MIT License.