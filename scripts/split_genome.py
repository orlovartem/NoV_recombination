import argparse
import copy
import os
import pandas as pd
from Bio import SeqIO


def split_genome(input_file, coord_file):
    '''
    Splits nucleotide sequences from input_file into ORFs using their
    coordinates from coord_file

    Input:
        input_file - file with nucleotide sequences in fasta-format
        coord_file - text file with coordinates of NTRs
    '''

    coord_df = pd.read_csv(coord_file, sep=",", index_col=0)
    orfs = list(coord_df.columns)

    print(coord_df)
    dict_orfs = {}
    for orf in orfs:
        dict_orfs[orf] = list()

    seqs = SeqIO.parse(open(input_file), format='fasta')

    for seq in seqs:
        # accession number of sequence
        acc = seq.id.split("_")[0]

        if acc == 'NC' or acc == 'AC':
            acc = '_'.join([acc, seq.id.split("_")[1]])
        print(acc)
        for orf in orfs:
            # checks whether all orfs have known coordinates
            if ((acc in coord_df.index) and (coord_df.loc[acc] != 'NA-NA').all()):

                # coordinates
                cd = coord_df[orf][acc]
                s, e = [int(x) for x in cd.split('-')]
                seq_orf = copy.deepcopy(seq)
                seq_orf.seq = seq.seq[s:e]
                dict_orfs[orf].append(seq_orf)
            else:
                print('{} not in coord file'.format(acc))
    out_temp = os.path.splitext(input_file)[0]
    for orf in dict_orfs.keys():
        outfile = out_temp + '_' + orf + '.fasta'
        SeqIO.write(dict_orfs[orf], outfile, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-coord", "--coord_file", type=str,
                        help="Csv-file with coordinates of 5\'UTR, coding region and 3\'UTR",
                        required=True)
    args = parser.parse_args()

    split_genome(args.input_file, args.coord_file)
