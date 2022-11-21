import argparse
import csv
import os
import re
import sys
from textwrap import wrap


def parse_gb(input_file, min_length, max_length, field_names=[]):
    '''
    For each entry in GenBank file retrieves nucleotide sequence, collection
    date and country generates file with sequences in fasta format.
    Sequences' names are in the following format:
        GenbankAccessionNumber_country_collectionYear
    Saves the file in the directory of input_file

    input_file - file with sequences in GenBank format
    min length, max length - minimal and maximal lenghts ofsequence in originfield
    field_names - list with qualifiers to be retrieved
    '''
    # File with countries/cities/other features and their abbreviations
    # Should be located at the same directory as script
    COUNTRY_MAP_FILE = os.path.join(sys.path[0], "country_map.csv")
    CITY_MAP_FILE = os.path.join(sys.path[0],"city_map.csv")

    # Output file
    OUTPUT_FILE = '.'.join(input_file.split('.')[:-1]) + '.fasta'

    # Length tresholds for sequences extracted from origin field

    MIN_ORIGIN_SIZE = min_length
    MAX_ORIGIN_SIZE = max_length

    # fields in GB entry which will be extracted
    if field_names == []:
        FIELD_NAMES = ["strain", "isolate",  "country", "collection_date"]
    else:
        FIELD_NAMES = field_names

    # File with abbreviations of countries
    country_map = read_csv(COUNTRY_MAP_FILE)

    if "host" in FIELD_NAMES:
        host_map = read_csv(os.path.join(sys.path[0], "host_map.csv"),
                            strip_it=True)
        host_map = compile_feature_map(host_map)

    # RegExp for parsing INPUT_FILE
    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)")  # accession number
    features = re.compile(r"^\s+\/([a-z_]+)=\"([^\"]+)\"$")  # qualifiers from FEATURE source field
    origin = re.compile(r"\s+\d+\s+([atgcnrykmswbdhvu\s]+)$")  # origin field (contains nt sequence)
    end_line = re.compile(r"^//\s+$")  # the end line of entry
    clean_from = re.compile(r"[:;.,/\s_]")
    clean_to = "-"

    def check_year(stri):
        year0 = re.compile(r"[0-9]{4}")
        m = year0.search(stri)
        if m:
            return year0.search(m.group()).group()

    # Parsing INPUT_FILE

    tests = []  # stores entries
    tests_nodate = []  # list for entries with no collection date

    test_accession = ""  # accession number of entry
    test_features = {}  # dictionary for qualifiers from FEATURE source field
    test_origin = ""  # ORIGIN field

    if not os.path.exists(input_file):
        print('{filename}: No such file or directory'.format(filename=input_file))
        return

    with open(input_file, "r") as in_f:
        for line in in_f:

            m = accession.match(line)  # finds ACCESSION field using RegExp
            if (m):
                test_accession = m.group(1)
            # check whether /source has started or has ended
            if re.match(r"^\s{5}[a-zA-z_]+[\s0-9<>\.]+", line) or re.match(r"^[a-zA-z]+", line):
                s = 0  # we are outside source field
            if re.match(r"^\s+source[0-9\.\s]+", line):
                s = 1  # inside source field
            if s == 1:
                m = features.match(line)  # finds all features in field 'FEATURES source' using RegExp
                if (m):
                    test_features[m.group(1)] = m.group(2)

            m = origin.match(line)  # ORIGIN field  (contains na sequence)
            if(m):
                test_origin += re.sub(r"\s+", "", m.group(1))

            if(end_line.search(line)):  # the end of entry
                if ("collection_date" or "collected_by") in test_features:
                    test_features["collection_date"] = check_year(test_features["collection_date"])  # reformates collection date

                if "country" in test_features:
                    test_features["country"] = map_feature(test_features["country"], country_map)  # replaces country by its abbreviation

                if "host" in FIELD_NAMES and "host" in test_features:
                    test_features["host"] = map_feature_reg(test_features["host"], host_map)

                if ("collection_date" or "collected_by") not in test_features:  # adds entries without collection date to test_nodate list
                    tests_nodate.append({
                        "accession": test_accession,
                        "features": test_features,
                        "origin": test_origin
                        })

                else:
                    if test_accession.split('_')[0] == 'NC':
                        test_accession = 'NC-' + test_accession[3:]
                    tests.append({
                        "accession": test_accession,
                        "features": test_features,
                        "origin": test_origin
                        })

                # resets variables
                test_accession = ""
                test_features = {}
                test_origin = ""
    in_f.close()

    # adds entries without collection date to the end of tests list
    for ent in tests_nodate:
        tests.append(ent)

    # Writing sequences in fasta-format to OUTPUT_FILE
    out = open(OUTPUT_FILE, "w+")

    for test in tests:
        aaccession, f, origin = test["accession"], test["features"], test["origin"]
        if not (MIN_ORIGIN_SIZE < len(origin) < MAX_ORIGIN_SIZE):
            continue

        # new name of entry = GB accession number + country + year
        output_name = ">%s" % aaccession

        for field_name in FIELD_NAMES:
            if field_name not in f:
                f[field_name] = 'NA'
            output_name += "_%s" % re.sub(clean_from, clean_to, f[field_name])
        output_name += "\n"
        out.write(output_name)

        # writes sequence from ORIGIN field
        out.write('\n'.join(wrap(origin, 60)) + '\n')

    out.close()
    return OUTPUT_FILE


def read_csv(file_name, strip_it=True):
    if not os.path.exists(file_name):
        return {}

    def strip(value):
        return strip_it and value.strip() or value

    with open(file_name) as csvfile:
        reader = csv.DictReader(csvfile,
                                delimiter=",",
                                fieldnames=["base", "new"])
        result = {}
        for row in reader:

            result[strip(row["base"])] = strip(row["new"])
        csvfile.close()
        return result


def compile_feature_map(feature_map):
    # feature_map is list wih tuples (compiled regular exp, value)
    feature_map_comp = []
    for key in feature_map:
        key_comp = re.compile(key)
        feature_map_comp.append((key_comp, feature_map[key]))
    return feature_map_comp


def map_feature_reg(feature, feature_map):
    '''
    feature - str - value in qualifier
    feature map - list with tuples, e.g. (reg_exp, value)
    If reg exp  matches feature, return value, else return feature
    '''
    for k, v in feature_map:
        if k.search(feature):
            return v
    print('Couldn\\t find key for {}'.format(feature))
    return feature


def map_feature(feature, feature_map):
    '''
    feature_map - dictionary, e.g. feature_map['Italia']='ITA'
    feature - possible key from feature_map
    return value if key is in feature_map
    '''
    for k, v in feature_map.items():
        if feature.lower().startswith(k.lower()):
            return v
    return feature


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-min", "--min_length", type=int,
                        help="Minimal length of sequence.\
                        Sequences shorter than min length will not be included in the final dataset",
                        required=True)
    parser.add_argument("-max", "--max_length", type=int,
                        help="Maximal length of sequence. \
                        Sequences longer than max length will not be included in the final dataset",
                        required=True)

    parser.add_argument("-f", "--features", type=str,
                        help="string with qualifiers to retrieve from GenBank annotation,\
                        e.g. 'country,host,collection_date'",
                        required=True)
    args = parser.parse_args()
    features = args.features.split(',')

    parse_gb(args.input_file, args.min_length, args.max_length, features)
