#!/usr/bin/env python

import argparse
import sys

try:
    from Bio import AlignIO
    from Bio import SeqIO
    from Bio import Alphabet
    from Bio import codonalign
except ImportError:
    print("Biopython is required. Try with: pip install biopython")
    sys.exit()
else:
    # args
    parser = argparse.ArgumentParser(prog="CodonAlign.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    Converts an amino acid alignment to nucleotides respecting codon positions\n""",
    epilog="""
    Examples:
    python CodonAlign.py -p myAlignment.aa.fas -c myCDSseqs.nt.fas -o CDS.aln.fas
    python CodonAlign.py -p myAlignment.aa.fas -c myCDSseqs.nt.fas -s

    Note that amino acid and nucleotide sequences should correspondingly have the same label.\n""")
    parser.add_argument(
    '--prot', '-p', type=str, default="", required=True,
    help='an amino acid alignment in FASTA format.')
    parser.add_argument(
    '--cds', '-c', type=str, default="", required=True,
    help='''un aligned CDS sequences FASTA format. These must be in-frame and correspond
    to the untranslated version of the amino acid sequence.''')
    parser.add_argument(
    '--outfile', '-o', default="aln.cds.fasta", nargs="?", type=str,
    help='name for output file (default: %(default)s)')
    parser.add_argument(
    '--stdout', '-s',  action="store_true", default=False,
    help='if sequences should be printed to screen.')
    args = parser.parse_args()

    # code
    # read data
    aln = AlignIO.read(args.prot, format="fasta", alphabet=Alphabet.IUPAC.protein)
    seq = list(SeqIO.parse(args.cds, format="fasta", alphabet=Alphabet.IUPAC.unambiguous_dna))

    # get names
    names_aln = [ s.name for s in aln ]
    names_seq = [ s.name for s in seq ]

    # check number of seqs
    if len(names_aln) != len(names_seq):
        sys.exit("Number of sequences in both files does not match.")

    # check names
    if all([ name in names_aln for name in names_seq]) and all([ name in names_seq for name in names_aln]):
        # build alignment
        codon_aln = codonalign.build(aln, seq)
        # delete the <unknown description> label
        for i in range(len(codon_aln)): codon_aln[i].description = ""
        if args.stdout:
            # print to screen
            print(codon_aln.format("fasta"))
        else:
            # print to file
            AlignIO.write(codon_aln, args.outfile, "fasta")
            print(f"{len(names_aln)} aligned CDS sequences saved to {args.outfile}.")
    else:
        sys.exit("Amino acid and nucleotide sequences do not have the same labels.")
