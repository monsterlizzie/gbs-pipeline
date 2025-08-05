#!/usr/bin/env python3
import argparse, sys, re


def get_targets(targets_file):
    """Read target text file into targets list"""
    targets = []
    with open(targets_file, 'r') as txt:
        for line in txt:
            targets.append(line.split('\n')[0])
    return targets


def write_line(line, target, flag, out):
    """Write line of matching targets"""
    if line[0] == '>':
        if '>{}\n'.format(target) == line:
            out.write(line)
            flag = 1
        else:
            flag = 0
    else:
        if flag:
            out.write(line)
    return flag


def write_fasta_file(fasta_file, target, output_prefix):
    """Write a FASTA file for the target"""
    with open(output_prefix + target + '_ref.fna', 'w') as out:
        with open(fasta_file, 'r') as fasta:
            flag = 0
            for line in fasta:
                flag = write_line(line, target, flag, out)


def write_target_fasta_files(targets, fasta_file, output_prefix):
    """Write a FASTA file for each target from targets list"""
    for target in targets:
        write_fasta_file(fasta_file, target, output_prefix)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get targets from res db.')
    parser.add_argument('--fasta_file', '-f', dest='fasta', required=True,
                        help='Input FASTA file.')
    parser.add_argument('--target_file', '-t', dest='target', required=True,
                        help='Input target text file.')
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output prefix.')
    return parser


def main():
    args = get_arguments().parse_args()

    # Get list of target names from target text file
    targets = get_targets(args.target)

    # Write FASTA file for each target specified 
    write_target_fasta_files(targets, args.fasta, args.output)


if __name__ == "__main__":
    sys.exit(main())
