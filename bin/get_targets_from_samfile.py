#!/usr/bin/env python3
import argparse, sys, re


def get_targets(targets_file):
    """Read target text file into targets list"""
    targets = []
    with open(targets_file, 'r') as txt:
        for line in txt:
            targets.append(line.split('\n')[0])
    return targets


def in_line(line, target):
    """Return True if headers and mappings of target sequence is in main SAM file, otherwise False"""
    if (bool(re.search('^@HD|^\@SQ.*' + target + '|^@PG', line))==True):
        return True
    elif (line.split('\t')[2] == target):
        return True
    else:
        return False


def write_sam_file(sam_file, target, id, output_prefix):
    """Write a SAM file for the ID and target"""
    with open(output_prefix + target + '_' + id + '_seq.sam', 'w') as out:
        with open(sam_file, 'r') as sam:
            for line in sam:
                if in_line(line, target):
                    out.write(line)


def write_target_sam_files(targets, sam_file, id, output_prefix):
    """Write a SAM file for each ID and target from targets list"""
    for target in targets:
        write_sam_file(sam_file, target, id, output_prefix)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get targets from sam file.')
    parser.add_argument('--sam_file', '-s', dest='sam', required=True,
                        help='Input sam file.')
    parser.add_argument('--target_file', '-t', dest='target', required=True,
                        help='Input target text file.')
    parser.add_argument('--id', '-i', dest='id', required=True,
                        help='Read ID.')
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output prefix.')
    return parser


def main():
    args = get_arguments().parse_args()

    # Get list of target names from target text file
    targets = get_targets(args.target)

    # Write SAM file for each ID and target specified
    write_target_sam_files(targets, args.sam, args.id, args.output)


if __name__ == "__main__":
    sys.exit(main())
