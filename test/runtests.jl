using BioDemultiplexer
using Test

function test_classification()
    query = "GG"
    ref = "AAAAAGGAA"

    score = semiglobal_alignment(query, ref)
    @test score == 1.0
end


function test_demultiplexing()
    # Test inputs
    tmpdir = mktempdir()

    input = joinpath(tmpdir, "test_input.fastq")
    input2 = joinpath(tmpdir, "test_input2.fastq")
    bc_tsv = joinpath(tmpdir, "test_bc.tsv")
    output_dir = joinpath(tmpdir, "test_output")

    open(input, "w") do io
        println(io, "@example_read1")
        println(io, "AATTATATGGGGGGATATT")
        println(io, "+")
        println(io, rand('A':'J', 19))
        println(io, "@example_read2")
        println(io, "AATTATATGCCGCCATATT")
        println(io, "+")
        println(io, rand('A':'J', 19))
    end

    open(input2, "w") do io
        println(io, "@example_read1")
        println(io, "ATGC")
        println(io, "+")
        println(io, "FFFF")
        println(io, "@example_read2")
        println(io, "TTTTT")
        println(io, "+")
        println(io, "FFF:F")
    end

    io = open(bc_tsv, "w")
    data = [
        ["ID", "Number", "Full_seq", "Full_annotation"],
        ["sample_1", 1, "AAAGGGGGGAAAT", "LLLBBBBBB3333"],
        ["sample_2", 2, "AACCCCCCAAA", "LLBBBBBB333"],
        ["sample_3", 3, "AAAGGCGGCAT", "LLLBBBBBB33"]
    ]
    for row in data
        row_str = join(string.(row), "\t")
        println(io, row_str)
    end
    close(io)

    # Run function
    execute_demultiplexing(input, input2, bc_tsv, output_dir)

    # Test outputsq
    @test isdir(output_dir)
    @test !isfile(joinpath(output_dir, "unknown.fastq"))
    @test !isfile(joinpath(output_dir, "sample_1.fastq"))
    @test isfile(joinpath(output_dir, "trimmed-sample_2.fastq"))
    @test isfile(joinpath(output_dir, "trimmed-sample_3.fastq"))
    open(joinpath(output_dir, "trimmed-sample_2.fastq"),"r") do out
        @test readline(out) == "@example_read1"
        @test readline(out) == "ATGC"
        @test readline(out) == "+"
        @test readline(out) == "FFFF"
        @test readline(out) == "@example_read2"
        @test readline(out) == "TTTTT"
        @test readline(out) == "+"
        @test readline(out) == "FFF:F"
    end
    open(joinpath(output_dir, "trimmed-sample_3.fastq"),"r") do out
        @test readline(out) == "@example_read1"
        @test readline(out) == "ATGC"
        @test readline(out) == "+"
        @test readline(out) == "FFFF"
        @test readline(out) == "@example_read2"
        @test readline(out) == "TTTTT"
        @test readline(out) == "+"
        @test readline(out) == "FFF:F"
    end
end
@testset "BioDemultiplexer.jl" begin
    @testset "classification tests" begin
        test_classification()
    end

    @testset "demultiplexing tests" begin
        test_demultiplexing()
    end
end
