#!/usr/bin/env python3
import argparse
import re
import sys
import time
import urllib.error
import urllib.request

# TODO: translate nucleotides to protien (all 6 frames) for querying similar to tblastn
# TODO: output option for saving results file
# TODO: output option for displaying results, output delimiter

# maybe add fasta validation or at the very least more friendly errors


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='smqcli',
        description='Find sequence motifs.')

    parser.add_argument(
        'motifs',
        help='Supports python\'s regular expression formatting. ' + \
             '"X" is also substituted as a wildcard.')

    sequence_type_group = parser.add_mutually_exclusive_group(required=False)
    sequence_type_group.add_argument(
        '-p', '--protein',
        action='store_true',
        help='Sequences and motifs are read as protiens. This option ' + \
             'is mutually exclusive to "-n"')
    sequence_type_group.add_argument(
        '-n', '--nucleotide',
        action='store_true',
        help='NOT IMPLEMENTED. Sequences and motifs are read as nucleotides. ' + \
             'This option is mutually exclusive to "-p"')

    option_group = parser.add_mutually_exclusive_group(required=False)
    option_group.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress stdout.')
    option_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Output additional information to stdout or stderr.')

    input_group = parser.add_argument_group(
        title='input group')
    input_group.add_argument(
        '-r', '--raw-sequence',
        nargs='*',
        help='Raw sequence(s) separated by commas or spaces.')
    input_group.add_argument(
        '-i', '--ifile', '--inputfile', '--input-file', '--fasta',
        nargs='*',
        help='File containing fasta formatted sequence. ' + \
             'Multiple files should be seperated by commas or spaces. ' + \
             '*Commas in file names will cause unintended results.*',
        metavar='FILENAME')
    input_group.add_argument(
        '-a', '--accession',
        nargs='*',
        help='GenBank acessions number(s) will be pulled from NCBI. ' + \
             'GI numbers (genInfo Identifier) also work. Multiple' + \
             'accession numbers should be seperated by commas or spaces.')

    output_group = parser.add_argument_group(
        title='output group')
    output_group.add_argument(
        '-s', '--standard',
        action='store_true',
        help='DEFAULT. Output is sent to stdout.')
    output_group.add_argument(
        '-o', '--ofile', '--outputfile', '--output-file',
        nargs='*',
        help='NOT IMPLEMENTED. Will output to file in .csv format.',
        metavar='FILENAME')
    # TODO: move from output group to regular option
    output_group.add_argument(
        '-d', '--delimiter',
        default='\t',
        help='This will be the delimiter for output. Special ' + \
        ' characters must be escaped. The default is \\t (tab). ' + \
        'For escaped charcters such as tab, use the syntax $\'\\t\'')
    args = parser.parse_args()

    if (args.protien is False and
            args.nucleotide is False):
        args.protien = True

    if (args.raw_sequence is None and
            args.ifile is None and
            args.accession is None):
        # Get the argument names for input_group and format for error.
        input_arg_names = []
        for store_action in input_group._group_actions:
            input_arg_names.append(store_action.option_strings[0])
        parser.error(
            'At least one argument in the input group {0} is required.'.format(
            '(' + ' | '.join(input_arg_names) + ')'))

    if (args.ofile is None and
            args.standard is False):
        args.standard = True

    sequences = {} # dictionary where {key: value} is {accession_number: sequence}
    matches = [] # list of lists with information on each match. will refactor into Match class
    delimiter = args.delimiter

    # prep motifs
    compiled_motifs = []
    for raw_motif in args.motifs.split(','):
        compiled_motifs.append(motif_to_regex(raw_motif))

    # prep sequences
    if args.raw_sequence:
        count = 0
        for sequence in split_args(args.raw_sequence):
            sequences['raw_sequence' + str(count)] = sequence
            count += 1
    if args.ifile:
        for file in split_args(args.ifile):
            try:
                f = open(file)
            except FileNotFoundError:
                print('FileNotFoundError on "{0}"'.format(file))
                sys.exit()
            else:
                with f:
                    add_fasta_to_sequences(f.read(), sequences)
    if args.accession:
        if args.verbose:
            num_accessions = len(split_args(args.accession))
            count_accessions = 1
            print(
                'Searching entrez for accession number. Will wait 0.5s ' + \
                'between each search so entrez server doesn\'t get angry')
        for an in split_args(args.accession):
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' + \
                  '?db=nuccore&rettype=fasta&id=' + an
            req = urllib.request.Request(url)
            try: fasta = urllib.request.urlopen(req).read().decode('utf-8')
            except urllib.error.URLError as e:
                if args.verbose:
                    print("error on {0}: {1}".format(an, e.reason))
                    count_accessions += 1
            except urllib.error.HTTPError as e:
                if args.verbose:
                    print("{0} error on {1}: {2}".format(e.code, an, e.reason))
                    count_accessions += 1
            else:
                add_fasta_to_sequences(fasta, sequences)
                # TODO: add argument for sleep time
                time.sleep(0.5)
                if args.verbose:
                    print("[%s/%s] %s "%(count_accessions, num_accessions, an))
                    count_accessions += 1

    for motif in compiled_motifs:
        for key in sequences.keys():
            for hit in motif.finditer(sequences[key]):
                # TODO: replace 'matches' with object
                matches.append([key, motif, hit])

    # TODO: add arguments for output format. e.g. --format=id,motif,hit,location
    if not args.quiet:
        #print("ID{0}MOTIF{0}HIT{0}START{0}END".format(delimiter))
        for match in matches:
            print(
                "{1}{0}{2}{0}{3}{0}{4}{0}{5}".format(
                delimiter,
                match[0],
                match[1].pattern,
                match[2].group(0),
                match[2].span()[0],
                match[2].span()[1]))

def motif_to_regex(raw_motif):
    # X to wildcard
    motif = re.sub('[Xx]', '.', raw_motif) # X to wildcard
    return re.compile(motif, re.IGNORECASE)

def split_args(arg_list):
    # split args on comma
    if len(arg_list) == 1:
        return arg_list[0].split(',')
    return arg_list

def split_fasta(fasta):
    # split on the first newline
    return fasta.split('\n', 1)

def add_fasta_to_sequences(fasta, sequences):
    for fasta_form in fasta.split('>')[1:]:
        # split a fasta formatted sequence and store it in data structure
        id = split_fasta(fasta_form)[0]
        sequence = split_fasta(fasta_form)[1].replace('\n', '')
        sequences[id] = sequence

if __name__ == "__main__":
    main()
