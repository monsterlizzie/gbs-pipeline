#!/usr/bin/env python3
import argparse, sys


def get_mismatch_and_depth(file):
    """Read mlst file"""
    with open(file, 'r') as txt:
        next(txt)
        for line in txt:
            values = line.split('\t')
            if values[1] != "failed":
                mismatch_depth = (str(values[9]), float(values[11]), "ST-{}".format(str(values[1])))
            else:
                mismatch_depth = 0
    return mismatch_depth


def write_alleles_file(alleles, output_file):
    with open(output_file, 'w') as out:
        for allele in alleles:
            out.write(allele + '\n')


def get_new_and_existing_alleles(mismatch_depth, min_read_depth, output_prefix):
    alleles = []
    if mismatch_depth[0] != '0':
        if mismatch_depth[1] >= min_read_depth:
            mismatches = mismatch_depth[0].split(';')
            alleles = ['Alleles found'] + [mismatch.split('/')[0] for mismatch in mismatches]
        else:
            alleles = [output_prefix + ': No new MLST alleles were found with sufficient read depth above ' + str(min_read_depth) + '.']
        write_alleles_file(alleles, "{}_new_mlst_alleles.txt".format(output_prefix))
    else:
        if mismatch_depth[1] >= min_read_depth:
            alleles = ['ID\tST', '{}\t{}'.format(output_prefix, mismatch_depth[2])]
        else:
            alleles = ['ID\tST', '{}\tNone found'.format(output_prefix)]
        write_alleles_file(alleles, "{}_existing_mlst_alleles.txt".format(output_prefix))


def get_arguments():
    parser = argparse.ArgumentParser(description='Get targets from sam file.')
    parser.add_argument('--mlst_results_file', '-m', dest='mlst', required=True,
                        help='Input SRST2 MLST results file.', type = str)
    parser.add_argument('--min_read_depth', '-d', dest='min_depth', required=True,
                        help='Minimum read depth threshold.', type = int)
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output file of alleles.', type = str)
    return parser


def main():
    args = get_arguments().parse_args()

    # Get contents of mlst results file
    mismatch_depth = get_mismatch_and_depth(args.mlst)

    # Get alleles from mismatch
    get_new_and_existing_alleles(mismatch_depth, args.min_depth, args.output)


if __name__ == "__main__":
    sys.exit(main())
