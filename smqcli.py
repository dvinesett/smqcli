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
# TODO: add option for regular expression queries rather than custom PROSITE-like one
#       as well as regular output
# TODO: reduce memory overhead on regex search by matching against [i:something] rather than [i:] ...
#       ... in various places

# PACKAGE_NAME = sys.modules[__name__]
PACKAGE_NAME = 'smqcli'

PROTEIN_CHARS = [
    'A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q',
    'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
    'P', 'S', 'T', 'W', 'Y', 'V'
]
PROTEIN_CHARS_STR = ''.join(PROTEIN_CHARS)


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
    compiled_patterns = []
    for raw_motif in args.motifs.split():
        compiled_patterns.append(motif_to_pattern(raw_motif))

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

    for pattern in compiled_patterns:
        # useful for debugging
        # print(",".join([str(subpattern) for subpattern in pattern.subpatterns]) +
        #       "; len: " + str(len(pattern.subpatterns)))
        for key in sequences.keys():
            for hit in pattern.findall(sequences[key]):
                # TODO: don't use this tuple any more. maybe add as attributes to _Match
                matches.append((key, pattern, hit))

    # TODO: add arguments for output format. e.g. --format=id,motif,hit,location
    if not args.quiet:
        for match in matches:
            print("{1}{0}{2}{0}{3}{0}{4}{0}{5}".format(
                delimiter,
                match[0],
                "".join([str(subpattern) for subpattern in match[1].subpatterns]),
                match[2].string,
                match[2].start,
                match[2].end
            ))


def motif_to_regex(raw_motif):
    # X to wildcard
    motif = re.sub('[Xx]', '.', raw_motif)  # X to wildcard
    return re.compile(motif, re.IGNORECASE)


def motif_to_pattern(raw_motif):
    return _Pattern(raw_motif)


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


class _Pattern:
    def __init__(self, raw_motif: str):
        self.subpatterns = []
        raw_motif = raw_motif.upper()
        i = 0
        while i < len(raw_motif):
            char = raw_motif[i]
            # chars
            if char in PROTEIN_CHARS:
                self.subpatterns.append(_CharSubpattern(char))
                i += 1
            # wildcard .X
            elif char == '.' or char == 'X':
                self.subpatterns.append(_WildcardSubpattern())
                i += 1
            # choice []
            elif char == '[':
                beg = i + 1
                try:
                    end = i + raw_motif[i:].index(']')
                except ValueError:
                    raise ValueError("Invalid motif. Could not find a closing bracket ']'")
                choices = set(raw_motif[beg:end])
                for choice in choices:
                    if choice not in PROTEIN_CHARS:
                        raise ValueError("Invalid motif. Characters in bracket starting at "
                                         "index {0} must be included in {1}".format(
                                            i, PROTEIN_CHARS_STR))
                self.subpatterns.append(_ChoiceSubpattern(choices, True))
                i = end + 1
            # not choice {}
            elif char == '{':
                beg = i + 1
                try:
                    end = i + raw_motif[i:].index('}')
                except ValueError:
                    raise ValueError("Invalid motif. Could not find a closing brace '}'")
                choices = set(raw_motif[beg:end])
                for choice in choices:
                    if choice not in PROTEIN_CHARS:
                        raise ValueError("Invalid motif. Characters in brace starting at "
                                         "index {0} must be included in {1}".format(
                                            i, PROTEIN_CHARS_STR))
                self.subpatterns.append(_ChoiceSubpattern(choices, False))
                i = end + 1
            # range ()
            elif char == '(':
                beg = i + 1

                try:
                    end = i + raw_motif[i:].index(')')
                except ValueError:
                    raise ValueError("Invalid motif. Could not find a closing parentheses ')'")

                try:
                    subpattern = self.subpatterns.pop()
                except IndexError:
                    raise ValueError("Invalid motif. Must have a subpattern before a range '()'")

                try:
                    range_groups = re.match(r"^(\d+)?(,)?(\d+)?$", raw_motif[beg:end]).groups()
                except AttributeError:
                    raise ValueError("Invalid motif. Range syntax must follow: 'x', 'x,y', 'x,', or ',y'")

                # init rgX variables so I don't have to worry about NameError
                rg0, rg1, rg2 = None, None, None

                # range_groups values from str to int
                try:
                    rg0 = int(range_groups[0])
                except TypeError:
                    pass
                rg1 = range_groups[1]
                try:
                    rg2 = int(range_groups[2])
                except TypeError:
                    pass

                # {x,y}
                if rg0 is not None and rg1 is not None and rg2 is not None:
                    self.subpatterns.append(_RangeSubpattern(subpattern, slice(rg0, rg2)))
                # {,y}
                elif rg0 is None and rg1 is not None and rg2 is not None:
                    self.subpatterns.append(_RangeSubpattern(subpattern, slice(0, rg2)))
                # {x,}
                elif rg0 is not None and rg1 is not None and rg2 is None:
                    raise TypeError("Don't search for range of patterns {0} in sequence with unbound stop.".format(
                        (rg0, rg2)))
                    # self.subpatterns.append(_RangeSubpattern(subpattern, slice(range_groups[0], None)))
                # {x}
                elif rg0 is not None and rg1 is None and rg2 is None:
                    self.subpatterns.append(_RangeSubpattern(subpattern, slice(rg0)))
                else:
                    raise NotImplementedError("The developer doesn't think you should be here. "
                                              "Contact him to rub it in.")
                i = end + 1
            else:
                raise ValueError("Invalid sequence expression. The input received is '{0}'. This error commonly "
                                 "occurs when the sequence expression is not contained by quotes".format(raw_motif))

    def __repr__(self):
        return "<{0}._Pattern object; {1}>".format(PACKAGE_NAME, ([str(subpattern) for subpattern in self.subpatterns]))

    def __len__(self):
        return len(self.subpatterns)

    def findall(self, matching_string):
        """ Return list of match objects for a Pattern in a sequence
        :param matching_string:
        :return: list of match objects
        """
        matches = []
        # loop through matching_string
        char_index = 0
        while char_index < len(matching_string) - len(self.subpatterns):
            submatches = []
            submatches_length = 0
            pattern_index = 0
            while pattern_index < len(self.subpatterns):
                # this is useful for debugging
                # print("chr {0}\tchrstr {1}\tmch {2}\tmchstr {3}\tptn {4}\tptnstr {5}".format(
                #         char_index, matching_string[char_index], char_index + submatches_length,
                #         matching_string[char_index + submatches_length], pattern_index,
                #         self.subpatterns[pattern_index]))
                match = self.subpatterns[pattern_index].match(
                    matching_string[(char_index + submatches_length):])
                if not match:
                    submatches.clear()
                    break
                submatches.append(match)
                submatches_length += len(match)
                pattern_index += 1

            if pattern_index == len(self.subpatterns):
                start = char_index
                stop = char_index + submatches_length
                matches.append(_Match(matching_string[start:stop], start, stop))
            char_index += 1
        return matches


class _Subpattern:
    def __init__(self, min_, max_):
        self.min = min_
        self.max = max_

    def match(self, s: str):
        raise NotImplementedError("Implemented in subclasses.")


class _CharSubpattern(_Subpattern):
    def __init__(self, char_: str):
        super().__init__(1, 1)
        if len(char_) == 1:
            self.char = char_
        else:
            raise ValueError("'char_' is a str but must have a length of 1.")

    def __str__(self):
        return self.char

    def match(self, sequence: str):
        if sequence[0] != self.char:
            return None
        else:
            return _Match(self.char, 0, 1)


class _WildcardSubpattern(_Subpattern):
    def __init__(self):
        super().__init__(1, 1)

    def __str__(self):
        return 'x'

    def match(self, sequence: str):
        if sequence == '':
            return None
        return _Match(sequence[0], 0, 1)


# includes [] and {}
class _ChoiceSubpattern(_Subpattern):
    def __init__(self, choices_: set, includes_=True):
        super().__init__(1, 1)
        self.choices = choices_
        self.includes = includes_

    def __str__(self):
        choices_str = "".join([str(choice) for choice in self.choices])
        if self.includes:
            # [x]
            return "[{0}]".format(choices_str)
        else:
            # {x}
            return "\{{0}\}".format(choices_str)

    def match(self, sequence: str):
        char = sequence[0]
        if self.includes:
            if char in self.choices:
                return _Match(char, 0, 1)
        elif not self.includes:
            if char not in self.choices:
                return _Match(char, 0, 1)
        else:
            return None


class _RangeSubpattern(_Subpattern):
    def __init__(self, subpattern, range_: slice):
        # TODO: refactor by implementing constants for each range type:
        if range_.step is not None and range_.step != 1:
            raise ValueError("'range_.step' is '{0}'. It should be 1 or None.".format(range_.step))
        self.subpattern = subpattern
        self.range_slice = range_
        min_range = None
        max_range = None
        if self.range_slice.start:
            min_range = subpattern.min * self.range_slice.start
        if self.range_slice.stop:
            max_range = subpattern.max * self.range_slice.stop
        if not min_range and not max_range:
            raise TypeError("Empty slice. This shouldn't happen. Please contact the developer.")
        elif not min_range and max_range:
            super().__init__(None, max_range)
            self.range_max_len = max_range
        elif min_range and not max_range:
            # self.range_max_len = None
            # super().__init__(min_range, None)
            raise NotImplementedError("Don't search for a range of patterns in sequence with unbound stop.")
        else:
            self.range_max_len = max_range
            super().__init__(min_range, max_range)

    def __str__(self):
        if self.range_slice.start:
            # x(2,5)
            return "{0}({1},{2})".format(self.subpattern, self.range_slice.start, self.range_slice.stop)
        else:
            # x(5)
            return "{0}({1})".format(self.subpattern, self.range_slice.stop)

    def match(self, sequence: str):
        i = 0
        match_str = ''
        total_match_len = 0
        while i < self.range_max_len:
            match = self.subpattern.match(sequence[i:])
            if not match:
                return None
            match_str += match.string
            match_range = match.end - match.start
            total_match_len += match_range
            i += 1
        return _Match(sequence[0:total_match_len], 0, i)


class _Match:
    def __init__(self, string_: str, start_: int, end_: int):
        self.string = string_
        self.start = start_
        self.end = end_

    def __str__(self):
        return "<{0}._Match object; span={1}, match={2}>".format(
            PACKAGE_NAME, (self.start, self.end), self.string)

    def __len__(self):
        return self.end - self.start


if __name__ == "__main__":
    main()
