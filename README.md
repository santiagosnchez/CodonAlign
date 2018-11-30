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
    wget 
