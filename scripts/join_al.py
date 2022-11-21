import argparse
import copy
from Bio import AlignIO


def join_alignments(input_file_n_list, out_name):
    '''
    input_file_n_list - list with filenames of alignments in fasta-format to join
    '''
    alignment = AlignIO.read(open(input_file_n_list[0]), "fasta")  # alignment object
    alignment_all = copy.deepcopy(alignment)

    for i in range(1, len(input_file_n_list)):
        print(input_file_n_list[i])
        alignment_temp = AlignIO.read(open(input_file_n_list[i]), "fasta")  # alignment object
        alignment_all = alignment_all + alignment_temp

    AlignIO.write(alignment_all, open(out_name, 'w'), "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_list", "--input_list", type=str,
                        help="List with names of input file", required=True)
    parser.add_argument("-out_name", "--out_name", type=str,
                        help="Name of output file", required=True)

    args = parser.parse_args()

    join_alignments(args.input_list.split(','), args.out_name)
