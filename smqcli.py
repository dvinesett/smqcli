#!/usr/bin/env python3
import argparse
import getopt
import re
import requests
import sys
import time

# TODO: cli input: choose between raw sequence, sequence in fasta file, sequence from ncib accession number (gi number works too)
# TODO: translate nucleotides to protien (all 6 frames) for querying similar to tblastn
# TODO: output option for saving results file
# TODO: output option for displaying results, output delimiter


def motif_to_regex(raw_motif):
    motif = re.sub('[Xx]', '.', raw_motif) # X to wildcard
    return re.compile(motif, re.IGNORECASE)

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
        '-p', '--protien',
        action='store_true',
        help='Sequences and motifs are read as protiens. This option ' + \
             'is mutually exclusive to "-n"')
    sequence_type_group.add_argument(
        '-n', '--nucleotide',
        action='store_true',
        help='NOT IMPLEMENTED. Sequences and motifs are read as nucleotides. ' + \
             'This option is mutally exclusive to "-p"')

    option_group = parser.add_mutually_exclusive_group(required=False)
    option_group.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='NOT IMPLEMENTED. Suppress stdout.')
    option_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='NOT IMPLEMENTED. Output additional information to stdout or stderr.')

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
             '*Commas in filenames will cause unintended results.*',
        metavar='FILENAME')
    input_group.add_argument(
        '-a', '--accession',
        nargs='*',
        help='GenBank acessions number(s) will be pulled from NCBI. ' + \
             'GI nubmers (genInfo Identifier)also work. Multiple' + \
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

    print(args)
    print()

    motifs = [] # list of unprocessed, unvalidated motifs, input by user
    accessions = [] # list of unvalidated accession numbers, input by user
    sequences = {} # dictionary where {key: value} is {accession_number: sequence}
    matches = [] # list of lists with information on each match. will refactor into Match class

    motifs = ['VXV','VXXV'] # arg
    accessions = ['ALJ99974', 'AIU94630'] # arg

    # TODO: validate motif

    compiled_motifs = [re.compile(motif) for motif in motifs]

    print('Will wait 0.5s between each search so entrez server doesn\'t get angry')
    num_accessions = len(accessions)
    count_accessions = 1
    for an in accessions:
        url = 'http://www.ncbi.nlm.nih.gov' + \
              '/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=' + an
        request = requests.get(url) # TODO: refactor using stdlib
        # TODO: implement error handling for bad url
        split = request.text.split('\n', 1) # split on first newline
        sequences[an] = split[1].replace('\n', '')

        for motif in motifs:
            for hit in motif_to_regex(motif).finditer(sequences[an]):
                matches.append([an, motif, hit, url])
                # TODO: replace 'matches' with object

        print("[%s/%s] %s "%(count_accessions, num_accessions, an))
        time.sleep(0.5)
        count_accessions += 1

    print('done looking for hits')

    print()
    for match in matches:
        print("%s,%s,%s,%s,%s" % (
            match[0], match[1], match[2].span(),
            match[2].group(0), match[3]))
        # TODO: use .format rather than %


# http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=940373999


if __name__ == "__main__":
    main()

# questions
# are nucleotide motifs a thing
# is there a standard format for the 6 frames of a protien sequence
#   e.g. ALJ99974.1 - ALF99974.6
