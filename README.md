# CodonAlign
Generates a codon-based alignment from an amino acid alignment and CDS sequences

## Description

This simple tool takes an amino acid alignment and converts it to a codon alignment,
provided a list of in-frame CDS sequences that translate to the aligned amino acids.

## Dependencies

It depends on [Biopython](https://biopython.org/) libraries: AlignIO, SeqIO, Alphabet, and codonalign.
Installation of Biopython can be done through `pip`:

    # for a system-wide install
    pip install biopython
    # or for a local install
    pip install --user biopython

## Installation

    git clone https://github.com/santiagosnchez/CodonAlign.git
    # or 
    wget https://raw.githubusercontent.com/santiagosnchez/CodonAlign/master/CodonAlign.py
    
## Running the code:

For the help message type:

    python CodonAlign.py -h

An example would be:

    python CodonAlign.py -p prot.aln.fas -c cds.fas -o cds.aln.fas
    
This example will save the output to a file and print a message to the screen.
Please note that `codonalign` is currently printing a warning to STDERR. This message can be avoided by using the `-s` or `--stdout` option:

    python CodonAlign.py -p prot.aln.fas -c cds.fas -s 2> /dev/null 1> cds.aln.fas

