#!/usr/bin/env python3
import argparse
import re
import sys
import time
import urllib.error
import urllib.request


# TODO: translate nucleotides to protein (all 6 frames) for querying similar to tblastn
# TODO: output option for saving results file
# TODO: output option for displaying results, output delimiter

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='smqcli',
        description='Find sequence motifs.')

    parser.add_argument(
        'motifs',
        help='Supports python\'s regular expression formatting. ' +
             '"X" is also substituted as a wildcard.')

    sequence_type_group = parser.add_mutually_exclusive_group(required=False)
    sequence_type_group.add_argument(
        '-p', '--protein',
        action='store_true',
        help='Sequences and motifs are read as proteins. This option ' +
             'is mutually exclusive to "-n"')
    sequence_type_group.add_argument(
        '-n', '--nucleotide',
        action='store_true',
        help='NOT IMPLEMENTED. Sequences and motifs are read as nucleotides. ' +
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
        help='File containing fasta formatted sequence. ' +
             'Multiple files should be separated by commas or spaces. ' +
             '*Commas in file names will cause unintended results.*',
        metavar='FILENAME')
    input_group.add_argument(
        '-a', '--accession',
        nargs='*',
        help='GenBank accession number(s) will be pulled from NCBI. ' +
             'GI numbers (genInfo Identifier) also work. Multiple' +
             'accession numbers should be separated by commas or spaces.')

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
        help='This will be the delimiter for output. Special ' +
             ' characters must be escaped. The default is \\t (tab). ' +
             'For escaped characters such as tab, use the syntax $\'\\t\'')
    args = parser.parse_args()

    if (args.protein is False and
                args.nucleotide is False):
        args.protein = True

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

    sequences = {}  # dictionary where {key: value} is {accession_number: sequence}
    matches = []  # list of lists with information on each match. will refactor into Match class
    delimiter = args.delimiter

    # prep motifs
    compiled_motifs = []
    # split on space
    for raw_motif in args.motifs.split():
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
                'Searching entrez for accession number. Will wait 0.5s ' +
                'between each search so entrez server doesn\'t get angry')
        for an in split_args(args.accession):
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' + \
                  '?db=nuccore&rettype=fasta&id=' + an
            req = urllib.request.Request(url)
            try:
                fasta = urllib.request.urlopen(req).read().decode('utf-8')
            except urllib.error.HTTPError as e:
                if args.verbose:
                    print("{0} error on {1}: {2}".format(e.code, an, e.reason))
                    count_accessions += 1
            except urllib.error.URLError as e:
                if args.verbose:
                    print("error on {0}: {1}".format(an, e.reason))
                    count_accessions += 1
            else:
                add_fasta_to_sequences(fasta, sequences)
                # TODO: add argument for sleep time
                time.sleep(0.5)
                if args.verbose:
                    print("[%s/%s] %s " % (count_accessions, num_accessions, an))
                    count_accessions += 1

    for motif in compiled_motifs:
        for key in sequences.keys():
            for hit in motif.finditer(sequences[key]):
                # TODO: replace 'matches' with object
                matches.append([key, motif, hit])

    # TODO: add arguments for output format. e.g. --format=id,motif,hit,location
    if not args.quiet:
        for match in matches:

            # hide wildcard ranges
            match_str = match[1].pattern[4:-2]
            out_str_list = []
            original_str = sequences[match[0]]
            tup = tuple(match[2].span(i) for i in range(2, len(match[2].groups()) + 1))
            i = 0
            if len(tup) > 0:
                out_str_list.append(original_str[match[2].span(1)[0]:tup[0][0]])
                while i < len(tup):
                    out_str_list.append("-({0})-".format(tup[i][1] - tup[i][0]))
                    if i == 0 and len(tup) != 1:
                        out_str_list.append(original_str[tup[0][1]:tup[1][0]])
                    elif i == len(tup) - 1:
                        pass
                    else:
                        out_str_list.append(original_str[tup[i][1]:tup[i+1][0]])
                    i += 1
                out_str_list.append(original_str[tup[i-1][1]:match[2].span(1)[1]])
                match_str = ''.join(out_str_list)
            ####################

            print(
                "{1}{0}{2}{0}{3}{0}{4}{0}{5}".format(
                    delimiter,
                    match[0],
                    # slicing is to get rid of (?=()) around motif
                    match_str,
                    match[2].span(1)[0],
                    match[2].span(1)[1],
                    match[2].string))



def motif_to_regex(raw_motif: str):
    # TODO tags for deletions '?'
    motif = raw_motif
    # remove parenthesis
    motif = re.sub(r'[\(\)]', '', motif)
    # X to wildcard
    motif = re.sub(r'[Xx]', '.', motif)
    # puts wildcards in groups
    # TODO format this verbosely
    # NOTE use these as test cases
    # V.V
    # V.+V
    # V.*V
    # V.{3}
    # V.{2,}
    # V.{,4}
    # V.{2,4}V
    # V.{2,4}.{5,6}V
    # V.{2,4}V.{5,6}V
    # V.{2,4}VV.{5,6}V
    # V.{2,4}VV.{5,}V.{3,7}.*....*.+.
    motif = re.sub(r'(\.[\{\+\*\?]((\d+)?(,\s*)?(\d+)?\})?[\+\*\?]?)+',
                   r'(\1)',
                   motif, re.VERBOSE)
    # Looks for overlapping results. I benchmarked it with a couple of use
    # cases and it only took a 10% performance hit
    return re.compile(r'(?=({0}))'.format(motif), re.IGNORECASE)


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
        desc = split_fasta(fasta_form)[0]
        sequence = split_fasta(fasta_form)[1].replace('\n', '')
        sequences[desc] = sequence


if __name__ == "__main__":
    main()
