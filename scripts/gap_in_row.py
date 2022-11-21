from Bio import SeqIO
import argparse
import re


def gap_in_row(input_file, gap_count):
    OUTPUT_FILE_NAME = ('.'.join(input_file.split('.')[:-1]) + '_{gaps}gp.fasta').format(gaps=gap_count)
    gaps_regexp = (r'.*\-{gap_count}.*'.format(gap_count='{'+str(gap_count)+'}'))
    gaps = re.compile(gaps_regexp)
    with open(input_file, 'r') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
        out_records = []
        count_all = 0
        count_rem = 0
        for rec in records:
            count_all += 1
            if re.match(gaps, str(rec.seq)):
                count_rem += 1
                print(rec.id.split('_')[0], 'removed')
            else:
                out_records.append(rec)
    with open(OUTPUT_FILE_NAME, 'w') as out_file:
        SeqIO.write(out_records, out_file, 'fasta')
    print('\n{count} records ({gaps} gaps in a row) removed'.format(count=count_rem, gaps=gap_count))
    print('{count1} ---> {count2}'.format(count1=count_all, count2=(count_all-count_rem)))
    print(OUTPUT_FILE_NAME)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-gap_count", "--gap_count", type=int,
                        help="Amount of gaps in a row", required=True)
    args = parser.parse_args()
    gap_in_row(args.input_file, args.gap_count)
