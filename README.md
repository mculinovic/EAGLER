# ONTscaffolder
Scaffolding using Oxford Nanopore reads

## Requirements
- g++ (4.8.2. or higher)
- make
- [Burrows-Wheeler Aligner][1] (0.7.12 or higher)

## Dependencies

- [SeqAn Library][2]

## Installation

To install the ONTscaffolder run the following commands from the folder where you want to install the tool:

	git clone https://github.com/mculinovic/ONTscaffolder.git
	cd ONTscaffolder/
	git submodule update --init --recursive
	make

Running the `make` command will create 2 binaries, `debug/main` and `release/main`,the debug and the release version of the tool respectively. A specific version may be built by running `make debug` or `make release`.

## Usage

- `usage: $0 <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta>`

[1]: https://github.com/lh3/bwa "Burrows-Wheeler Aligner"
[2]: https://github.com/seqan/seqan "SeqAn Library"
