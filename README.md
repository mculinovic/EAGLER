# ONTscaffolder

ONT Scaffolder is a scaffolding tool for Oxford Nanopore reads. The scaffolder takes as input a draft genome created by any NGS assembler and a set of Nanopore reads. The long reads are used to extend the contigs present in the NGS draft.

The tool should be compatible with most UNIX flavors and has been successfully test on the following operating systems:

- Mac OS X 10.10.3
- Ubuntu 14.04 LTS

## Requirements

To run the scaffolder correctly, all the executables for the programs listed below should be reachable from the `PATH` variable of your shell.

- g++ (4.8.2. or higher)
- [GNU Make][4]
- [Burrows-Wheeler Aligner][1] (0.7.12 or higher)
- [Doxygen][3] (optional)

## Dependencies

The scaffolder depends directly on 2 libraries: SeqAn and cpppoa. Both libraries will be automatically downloaded, configured and built by entering the [installation commands](#installation) in a terminal window.

- [SeqAn Library][2]
- [cpppoa][5]

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

## Documentation

The documentation for this tool was written to work with the doxygen documentation generator. To successfully generate the documentation, the doxygen executable must be in your `PATH` variable.

To create the documentation in HTML and LaTeX format run the following command from the root of the tool:

	make docs
	
HTML documentation is placed in `docs/html`, while the LaTeX documentation is placed in `docs/latex`. To view the HTML documentation open `docs/html/index.html` in any web browser. The PDF documentation is obtainable by compiling the generated LaTeX code with the provided makefile.

Use the following commands from the root of the project to create the PDF version of the documentation:

	cd docs/latex/
	make
	open refman.pdf
	
Please check the [links.md](links.md) file for links to datasets and additional information on the alignment toolbox used by the scaffolder.

## Usage

To run the tool please use the provided run script as show below:

	./run.sh <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta>
	
The implementation will automatically detect the number of hardware threads supported by the system 
	
###Arguments:

 1. **ont_reads.fasta**: FASTA file containing the Oxford Nanopore reads to be used in the scaffolding
 2. **draft_genome.fasta**: FASTA file containing the draft genome created by some NGS pipeline
 3. **output_file.fasta**: FASTA file with the extended and/or scaffolded contigs

## Scripts
 
Some utility scripts are available in the `scripts` folder. All scripts have been developed and tested with Python 3.4.3.

For detailed usage instructions run the following command for the dseired script:

	python3 scripts/<script_name>.py --help
	
###Available scripts:

 1. **genome2contigs**: cuts a reference genome into multiple contigs
 2. **extension_analysis**: runs various statistics on contig extensions produced by the scaffolder

## Contributors

- [Marko Culinovic](marko.culinovic@gmail.com)
- [Luka Sterbic](luka.sterbic@gmail.com)

[1]: https://github.com/lh3/bwa "Burrows-Wheeler Aligner"
[2]: https://github.com/seqan/seqan "SeqAn Library"
[3]: http://www.stack.nl/~dimitri/doxygen/ "Doxygen"
[4]: http://www.gnu.org/software/make/ "GNU Make"
[5]: https://github.com/mculinovic/cpppoa "cpppoa"
