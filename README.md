# ONTscaffolder
Scaffolding using Oxford Nanopore reads.

## Requirements
- g++ (4.8.2. or higher)
- make
- [Burrows-Wheeler Aligner][1] (0.7.12 or higher)
- doxgen (optional)
- graphviz (optional)

## Dependencies

- [SeqAn Library][2]

## Installation

To install the ONTscaffolder run the following commands from the folder where you want to install the tool:

	git clone https://github.com/mculinovic/ONTscaffolder.git
	cd ONTscaffolder/
	git submodule update --init --recursive
	make

Running the `make` command without arguments will build the release version of the tool as the binary file `release/scaffolder`. 

To build the debug version of the tool use:

	make debug
	
To build both the debug and release versions use:

	make all
	
Once the release version has been build you may run the following command to install the tool to `/usr/local/bin`:

	make install
	
To delete all files generated during the build process, for both debug and release, use:
 
	make clean

To remove the installed executable use:

	make uninstall

## Usage

To run the tool please use the provided run script as show below:

	./run.sh <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta>
	
The implementation will automatically detect the number of hardware threads supported by the system 
	
###Arguments:

 1. **ont_reads.fasta**: FASTA file containing the Oxford Nanopore reads to be used in the scaffolding
 2. **draft_genome.fasta**: FASTA file containing the draft genome created by some NGS pipeline
 3. **output_file.fasta**: FASTA file with the extended and/or scaffolded contigs

## Scripts
 
Some utility scripts are available in the `scripts` folder. All scripts have been developed and tested with Python 3.4.3. For detailed usage instructions run the following command for the dseired script:

	python3 scripts/<script_name>.py --help
	
###Available scripts:

 1. **genome2contigs**: cuts a reference genome into multiple contigs 


[1]: https://github.com/lh3/bwa "Burrows-Wheeler Aligner"
[2]: https://github.com/seqan/seqan "SeqAn Library"
