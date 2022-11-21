import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from Bio import SeqIO

# list with ambiguous nucleotides
ambig_nt = ['n',  # A, T, G, C
            'r',  # Purine (A or G)
            'y',  # Pyrimidine (T or C)
            'k',  # Keto (G or T)
            'm',  # Amino (A or C)
            's',  # Strong interaction (3 H bonds) (G or C)
            'w',  # Weak interaction (2 H bonds) (A or T)
            'b',  # Not A (C or G or T)
            'd',  # Not C (A or G or T)
            'h',  # Not G (A or C or T)
            'v',  # Not T or U (A or C or G)
            'total'
            ]


def resolve_ambiguos(input_file, output_dir, window):
    '''
    Resolves ambiguous nucleotides in nucleotide sequences according to
    consensus in most related sequences to the region with ambiguous nt.
    Removes sequences which

    Input:
        input_file - path to file with alignment of nt sequences in fasta-format
        out_dir - path to directory of output alignment without amb nucleotides
    '''

    if sys.platform == 'win32' or sys.platform == 'cygwin':
        output_dir += "\\"
    else:
        output_dir += "/"

    # alignment with nucleotide sequences in fasta-format
    fasta_al = list(SeqIO.parse(open(input_file), "fasta"))

    # the length of window around the ambiguous nucleotide to cut from original sequence
    window = 100

    # list with SeqIO objects which sequences contain less ambiguous characters that specified threshold
    fasta_al_less_amb = []

    # list with SeqIO objects which sequences
    # correspond to slice surrounding ambiguous character
    list_slices = []

    for rec in fasta_al.copy():
        # total number of ambiguous nucleotides in sequence
        amb_total = len(re.findall(r"[nrykmswbdhv]", str(rec.seq)))
        if amb_total == 0:
            # adds records with no ambiguous characters to the new alignment
            fasta_al_less_amb.append(rec)
        else:
            # checking whether the number of ambiguous characters exceed specified threshold
            if (amb_total/len(re.sub("-", "", str(rec.seq)))) > 0.01:
                print(rec.name, 'exceeded threshold')
                continue
            elif (re.findall(r"nnnnn", str(rec.seq))):
                print(rec.name, 'has 5 or more n-nt in a row')
                continue
            else:
                print(rec.name)
                # if the number of ambiguous nucleotides doesn't exceed the threshold
                # add copy record to a new list
                fasta_al_less_amb.append(rec)
                # finds positions of ambiguous nucleotides in sequence
                starts = [m.start() for m in re.finditer(r"[nrykmswbdhv]", str(rec.seq))]

                # for each ambiguous nt creates a slice with length=window surrounding this nt
                i = 0
                # list with start and end positions of slices
                slices = []
                while i < len(starts):
                    # starts of ambiguous nt in current slice
                    current_starts = []
                    current_starts.append(str(starts[i]+1))
                    # takes the start of sequence if amb nt is closer than window/2
                    # to the beginning of seq
                    if starts[i] < window/2:
                        st = 0
                        e = window
                    else:
                        # takes the end of the sequence if amb nt is closer than window/2
                        # to the end of sequences
                        if len(rec.seq) - starts[i] < window/2:
                            st = len(rec.seq) - window
                            e = len(rec.seq)
                        # takes the slice [starts[i]-window/2, starts[i]+window/2]
                        else:
                            st = int(starts[i]-window/2)
                            e = int(starts[i]+window/2)

                    # creates slice object
                    cur_slice_rec = copy.deepcopy(rec[st:e])
                    cur_slice_rec.description = ''

                    # appends sliced sequence to the list
                    list_slices.append(cur_slice_rec)

                    slices.append([st, e])

                    if i+1 < len(starts):
                        for j in range(i+1, len(starts), 1):
                            if st+int(window/5) < starts[j] and starts[j] < e-int(window/5):
                                current_starts.append(str(starts[j]+1))
                                i += 1
                                if j == len(starts) - 1:
                                    i += 1
                                continue
                            else:
                                i = i + 1
                                break
                    else:
                        i += 1
                    cur_slice_rec.id = rec.id + "_" + ":".join([str(st+1)]+current_starts+[str(e+1)])

    # filename for fasta-file with slices
    file_name_slices = os.path.splitext(input_file)[0] + "_slices.fasta"
    # writes slices to fasta_file
    SeqIO.write(list_slices, file_name_slices, "fasta")

    file_name_less_amb = os.path.splitext(input_file)[0] + "_less_amb.fasta"
    with open(file_name_less_amb, 'w') as file_less_amb:
        SeqIO.write(fasta_al_less_amb, file_less_amb, "fasta")
    file_less_amb.close()

    # commands for creating local database and blast slices against it
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        makeblast_command = 'makeblastdb -in {} -dbtype nucl -out {}local_db'.format(file_name_less_amb, output_dir)
        print(makeblast_command)
        blastn_command = 'blastn -db {out_path}local_db -query {input} -outfmt 6 -out \
                            {out_path}blast.out -strand plus -evalue 1e-20 -word_size 7'.format(
                            input=file_name_slices, out_path=output_dir)
        print(blastn_command)
    else:
        makeblast_command = 'makeblastdb -in {} -dbtype nucl -out {}local_db'.format(file_name_less_amb, output_dir)
        blastn_command = 'blastn -db {out_path}local_db -query {input} -outfmt 6 -out \
                            {out_path}blast.out -strand plus -evalue 1e-20 -word_size 7'.format(
                            input=file_name_slices, out_path=output_dir)
    subprocess.call(makeblast_command, shell=True)

    # blast against reference sequences
    subprocess.call(blastn_command, shell=True)

    # dataframe with blast results
    blast_output = pd.read_csv(output_dir+'blast.out', sep='\t', header=None,
                               names=['qseqid',
                                      'sseqid',
                                      'pident',
                                      'length',
                                      'mismatch',
                                      'gapopen',
                                      'qstart',
                                      'qend',
                                      'sstart',
                                      'send',
                                      'evalue',
                                      'bitscore'])
    blast_output.head()

    # create dictionary with sequences instead of list
    fasta_al_less_amb = SeqIO.to_dict(fasta_al_less_amb)

    # flag indicates whether the sequence has been resolved
    flag = 0

    current_seq_id = ''

    for row in blast_output.iterrows():
        # changes flag when meets new sequence with amb nt
        if row[1]['qseqid'] != current_seq_id:
            current_seq_id = row[1]['qseqid']
            # id of sequence with amb nt according to original fasta
            current_seq_id_orig = "_".join(row[1]['qseqid'].split('_')[:-1])
            # start position of slice (enumeration from 0)
            start = int(row[1]['qseqid'].split('_')[-1].split(':')[0]) - 1
            # end positions of slice
            end = int(row[1]['qseqid'].split('_')[-1].split(':')[-1]) - 1
            # ambiguous positions within slice
            left_amb_pos = row[1]['qseqid'].split('_')[-1].split(':')[1:-1]
            left_amb_pos = [int(x)-1 for x in left_amb_pos]

            flag = 0
            # skips the row if amb nts have been resolved
        if flag != 1:
            if row[1]['sseqid'] == '_'.join(row[1]['qseqid'].split('_')[:-1]):
                continue
            else:
                # relative start of ambiguous character in window
                rel_amb_pos = [x - start for x in left_amb_pos]
                # position corresponding to ambiguous character in reference sequence
                ref_pos = [int(row[1]['sstart']) - 1 + x for x in rel_amb_pos]
                # nucleotides in reference sequence in positions that are ambiguous
                ref_res_nuc = [fasta_al_less_amb[row[1]['sseqid']].seq[x] for x in ref_pos]

                left_amb_pos_copy = left_amb_pos.copy()
                # changes amb nt to the ones in the reference sequence
                for i in range(len(left_amb_pos)):
                    if ref_res_nuc[i] not in ambig_nt:
                        fasta_al_less_amb[current_seq_id_orig].seq = fasta_al_less_amb[current_seq_id_orig].seq[:left_amb_pos[i]]+ref_res_nuc[i]+fasta_al_less_amb[current_seq_id_orig].seq[left_amb_pos[i]+1:]
                        left_amb_pos_copy.remove(left_amb_pos[i])
                    else:
                        continue
                left_amb_pos = left_amb_pos_copy[:]
                if len(left_amb_pos) == 0:
                    flag = 1
                else:
                    print(current_seq_id, 'reference has amb')
    print(file_name_less_amb)
    SeqIO.write(fasta_al_less_amb.values(), file_name_less_amb, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input xml-file generated in BEAUTi", required=True)
    parser.add_argument("-pout", "--path_out", type=str,
                        help="Output directory. If not defined the output files will be saved \
                        in 'years' folder in the directory of input file")
    parser.add_argument("-w", "--window", type=str,
                        help="window size")
    args = parser.parse_args()

    if not args.path_out:
        args.path_out = os.path.split(args.input_file)[0]

    resolve_ambiguos(args.input_file, args.path_out, args.window)
