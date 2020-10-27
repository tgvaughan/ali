#!/usr/bin/env python3

from Bio import AlignIO
from argparse import ArgumentParser, FileType
from sys import stdin, stdout
import re

def cat(args):
    """Concatenate sequences and print in fasta format."""
    for infile in args.files:
        print(format(AlignIO.read(infile, args.format),args.output_format), file=args.output)
        
def summary(args):
    """Print summary of alignment contents."""
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        print("{name}: {nseqs} sequences of length {length}".format(
            name=infile.name,
            nseqs=len(align),
            length=align.get_alignment_length()),
              file=args.output)

def names(args):
    """Print names of sequences in alignment."""
    align = AlignIO.read(args.file, args.format)
    for record in align:
        print(record.id, file=args.output)


def repname(args):
    """Replace patterns in sequence names."""
    p = re.compile(args.pattern)
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        for record in align:
            record.id = p.sub(args.replacement, record.id)

        print(format(align, args.output_format), file=args.output)
        

if __name__ == '__main__':
    parser = ArgumentParser(description='Sequence file manipulation tool.')

    parser.add_argument('-f', '--format', dest="format", default='fasta', type=str,
                        help="Input format for sequences.")

    parser.add_argument('-t', '--output_format', dest="output_format", default='fasta', type=str,
                        help="Format used for output alignment (where applicable).")

    parser.add_argument('-o', '--output', dest="output", default=stdout, type=FileType('w'),
                        help="Output file.")

    subparsers = parser.add_subparsers(dest='command', help='Command help', required=True)

    parser_cat = subparsers.add_parser('cat', help="Concatenate alignments vertically.")
    parser_cat.add_argument('files', nargs='*', default=[stdin], type=FileType('r'),
                            help="Alignments to concatenate.")
    parser_cat.set_defaults(func=cat)

    parser_summary = subparsers.add_parser('summary', help="Summarize alignments.")
    parser_summary.add_argument('files', nargs='*', default=[stdin], type=FileType('r'),
                            help="One or more alignments to summarize.")
    parser_summary.set_defaults(func=summary)

    parser_names = subparsers.add_parser('names', help="List names of sequences in alignment.")
    parser_names.add_argument('file', nargs='?', default=stdin, type=FileType('r'),
                            help="Alignment file to read.")
    parser_names.set_defaults(func=names)

    parser_repname = subparsers.add_parser('repname', help="Replace string in names.")
    parser_repname.add_argument('pattern', type=str, help="Pattern to replace")
    parser_repname.add_argument('replacement', type=str, help="Replacement string")
    parser_repname.add_argument('files', nargs='*', default=[stdin], type=FileType('r'))
    parser_repname.set_defaults(func=repname)

    args = parser.parse_args()
    args.func(args)