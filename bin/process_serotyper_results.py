#!/usr/bin/env python3
import argparse
import sys

replace_values = {
    'GBS-SBG:': '',
    'A': 'a',
    'B': 'b',
    'IIIa': 'III',
    'IIIb': 'III',
    'IIIc': 'III',
    'IIId': 'III',
    'III-1': 'III',
    'III-2': 'III',
    'III-3': 'III',
    'III-4': 'III'
}

def make_gene_list(input_file, depth_threshold):
    """Get features from SRST2 input file into dictionary depending on read depth threshold"""
    gene_list = []
    with open(input_file, 'r') as fg_file:
        next(fg_file) # Skip header row
        for line in fg_file:
            feature = line.split('\t')
            if float(feature[5]) <= depth_threshold:
                feature[3] = 'NT'
            gene_list.append(feature[2:-1])
    return gene_list


def write_outfile(gene_list, out_file):
    """Write serotype, match type status and average read depth to output file"""
    matched_alleles = []
    match_type = []
    serotype = []
    avgdepth = []
    for values in gene_list:
        status = 'imperfect' if values[4] != '' else 'identical'
        value = values[1]
        for key, item in replace_values.items():
            value = value.replace(key, item)
        matched_alleles = matched_alleles + [value]
        match_type = match_type + [value + '=' + status]
        serotype = serotype + [value]
        avgdepth = avgdepth + [values[3]]
    if len(serotype) > 1 and 'NT' in serotype:
        serotype = [x for x in serotype if x != 'NT']
    with open(out_file, 'w') as out:
        out.write('Matched_Allele'+'\t'+'Match_Type'+'\t'+'Serotype'+'\t'+'AvgDepth'+'\n'+'/'.join(matched_alleles)+'\t'+'/'.join(match_type)+'\t'+'/'.join(serotype)+'\t'+'/'.join(avgdepth)+'\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify fullgenes output of SRST2.')
    parser.add_argument('--srst2_output', '-s', dest='id', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--sero_db', '-b', dest='db', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output filename.')
    parser.add_argument('--min_read_depth', '-d', dest='depth', default = 0, type=float, required=False,
                        help='Minimum read depth where mappings with fewer reads are excluded. Default: 0.')

    return parser


def main():
    args = get_arguments().parse_args()
    db_name = ' '.join(args.db.split('.')[:-1])

    # Specift fullgenes file path from ID and database name
    fullgenes_file = args.id + '__fullgenes__' + db_name + '__results.txt'

    # Get feature dictionary
    gene_list = make_gene_list(fullgenes_file, args.depth)

    # Write tab-delimited output file with serotype features
    write_outfile(gene_list, args.output)


if __name__ == "__main__":
    sys.exit(main())
