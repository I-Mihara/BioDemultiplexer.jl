# BioDemultiplexer
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://I-Mihara.github.io/BioDemultiplexer.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://I-Mihara.github.io/BioDemultiplexer.jl/dev/)
[![Build Status](https://github.com/I-Mihara/BioDemultiplexer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/I-Mihara/BioDemultiplexer.jl/actions/workflows/CI.yml?query=branch%3Amain)

BioDemultiplexer is a Julia package that provides functions for demultiplexing reads based on barcodes.
## Usage
Execute_demultiplexeing is the main function of this package. It is executed as follows.
```julia
function execute_demultiplexing(file_R1::String, file_R2::String, bc_file::String, output_dir::String)
```
### Arguments
* `file_R1::String`: the path to the input R1 fastq file
* `file_R2::String`: the path to the input R2 fastq file
* `bc_file::String`: the path to the reference barcode file in TSV format
* `output_dir::String`: the path to the output directory where the demultiplexed files will be written
### Optional Arguments
* `max_error_rate`: 
* `min_deltha`: 
* `bc_rev`: 
・単一ファイルに関して書く
## Parallel computing
This package is compatible with parallel computing.
## License
This package is licensed under the MIT License.