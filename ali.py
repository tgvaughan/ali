#!/usr/bin/env python3

from Bio import AlignIO
from argparse import ArgumentParser, FileType
from sys import stdin, stdout
import re

from random import sample

def cat(args):
    """Concatenate sequences and print in fasta format."""
    for infile in args.files:
        print(format(AlignIO.read(infile, args.format),args.output_format), file=args.output)
        
def summary(args):
    """Print summary of alignment contents."""
    total_seq_count = 0
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        total_seq_count += len(align)

        print("{name}: {nseqs} sequences of length {length}".format(
            name=infile.name,
            nseqs=len(align),
            length=align.get_alignment_length()),
              file=args.output)
    if len(args.files)>1:
        print("--- Total sequence count: {} ---".format(total_seq_count))

def names(args):
    """Print names of sequences in alignment."""
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        for record in align:
            print(record.description, file=args.output)

def rename(args):
    """Replace patterns in sequence names."""
    p = re.compile(args.pattern)
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        for record in align:
            record.description = p.sub(args.replacement, record.description)
            record.name = record.description
            record.id = ""

        print(format(align, args.output_format), file=args.output)

def grep(args):
    """Select sequences with pattern in name."""
    p = re.compile(args.pattern)
    if args.invert:
        pred = lambda s: p.search(s.id) == None
    else:
        pred = lambda s: p.search(s.id) != None
    
    for infile in args.files:
        align = AlignIO.read(infile, args.format)
        filtered = AlignIO.MultipleSeqAlignment(filter(pred, align))
        print(format(filtered, args.output_format), file=args.output)

def subsample(args):
    """Sample sequences with replacement."""
    sequences = sum([list(AlignIO.read(infile, args.format)) for infile in args.files], [])
    count = min(args.count, len(sequences))
    output_alignment = AlignIO.MultipleSeqAlignment(sample(sequences, count))
    print(format(output_alignment, args.output_format), file=args.output)

if __name__ == '__main__':
    parser = ArgumentParser(description='Sequence file manipulation tool.')

    parser.add_argument('-f', '--format', dest="format", default='fasta', type=str,
                        help="Input format for sequences. (Default fasta.)")

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

    parser_names = subparsers.add_parser('names', help="List names of sequences in alignments.")
    parser_names.add_argument('files', nargs='*', default=[stdin], type=FileType('r'),
                            help="Alignment files to read.")
    parser_names.set_defaults(func=names)

    parser_rename = subparsers.add_parser('rename', help="Replace string in names.")
    parser_rename.add_argument('pattern', type=str, help="Pattern to replace")
    parser_rename.add_argument('replacement', type=str, help="Replacement string")
    parser_rename.add_argument('files', nargs='*', default=[stdin], type=FileType('r'))
    parser_rename.set_defaults(func=rename)

    parser_grep = subparsers.add_parser('grep', help="Find sequences with matching names.")
    parser_grep.add_argument('pattern', type=str, help="Pattern to find")
    parser_grep.add_argument('files', nargs='*', default=[stdin], type=FileType('r'))
    parser_grep.add_argument('-v', dest='invert', action='store_true', help="Invert search.")
    parser_grep.set_defaults(func=grep)

    parser_subsample = subparsers.add_parser('subsample', help="Subsample alignment.")
    parser_subsample.add_argument('count', type=int, help="Number of sequences to sample")
    parser_subsample.add_argument('files', nargs='*', default=[stdin], type=FileType('r'))
    parser_subsample.add_argument('with-replacement', action='store_true', help='Sample with replacement.')
    parser_subsample.set_defaults(func=subsample)

    args = parser.parse_args()
    args.func(args)
