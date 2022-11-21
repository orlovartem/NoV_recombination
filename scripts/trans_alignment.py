import argparse
import os
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def translate_aln(input):
    out_file_n = os.path.splitext(input)[0] + '_tr.fasta'
    seqs = SeqIO.parse(open(input), 'fasta')
    # list with translated SeqRecords objects
    list_translated = []

    # dictionary with list of codons for sequences
    dict_codons = {}
    for seq in seqs:
        seq_aa = seq.translate()
        seq_aa.description = ''
        seq_aa.id = seq.id
        list_translated.append(seq_aa)
        if len(seq.seq) % 3 != 0:
            print("Len of seq in not a multiple of three for {}".format(seq.id))
        list_codons = re.findall(r'.{3}', str(seq.seq))
        if len(seq.seq) % 3 == 2:
            print(2)
            list_codons.append(str(seq.seq[-2:]))
        if len(seq.seq) % 3 == 1:
            print(1)
            list_codons.append(str(seq.seq[-1:]))
        dict_codons[seq.id] = list_codons
    SeqIO.write(list_translated, out_file_n, "fasta")
    return out_file_n, dict_codons


def reverse_translate_aln(input, dict_codons):
    out_file_n = os.path.splitext(input)[0] + '_rt.fasta'
    # list with SeqRecord objects
    reverse_tr_seqs = list()
    aa_seqs = list(SeqIO.parse(open(input), 'fasta'))
    max_length = len(aa_seqs[0].seq)*3 + 3
    print(max_length)
    for aa_seq in aa_seqs:
        new_seq = ''
        count = 0
        for aa in aa_seq.seq:
            if aa != '-':
                new_seq += dict_codons[aa_seq.id][count]
                count += 1
            else:
                new_seq += '---'
        dif_length = len(dict_codons[aa_seq.id]) - count
        if dif_length > 0:
            left_seq = ''.join(dict_codons[aa_seq.id][-dif_length:])
            new_seq += left_seq
        if len(new_seq) < max_length:
            print(len(new_seq), max_length)
            new_seq = new_seq + '-' * (max_length - len(new_seq))
        if len(new_seq) > max_length:
            print(len(new_seq), max_length)
        new_nt_seq = SeqRecord(Seq(new_seq), id=aa_seq.id)
        new_nt_seq.description = ''
        reverse_tr_seqs.append(new_nt_seq)
    SeqIO.write(reverse_tr_seqs, out_file_n, "fasta")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)

    args = parser.parse_args()

    trans_file, dict_codons = translate_aln(args.input_file)
    trans_file_aln = os.path.splitext(trans_file)[0] + '_aln.fasta'
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        os.system(('mafft --retree 1 ' + trans_file + ' > ' + trans_file_aln))
    else:
        os.system(('mafft --retree 1 ' + trans_file + ' > ' + trans_file_aln))
    reverse_translate_aln(trans_file_aln, dict_codons)
