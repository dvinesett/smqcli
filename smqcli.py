#!/usr/bin/env python3
import re, requests, sys, time

# TODO: cli input: choose between raw sequence, sequence in fasta file, sequence from ncib accession number (gi number works too)
# TODO: decide the best way to structure cli option formatting
# TODO: translate nucleotides to protien (all 6 frames) for querying similar to tblastn
# TODO: output option for saving results file
# TODO: output option for displaying results, output delimiter

def motif_to_regex(raw_motif):
    motif = re.sub('[Xx]', '.', raw_motif)
    return re.compile(motif, re.IGNORECASE)

motifs = ['VXV','VXXV'] # arg
accessions = ['ALJ99974', 'AIU94630'] # arg
sequences = {}
matches = []

# TODO: validate motif

compiled_motifs = [re.compile(motif) for motif in motifs]

print('Will wait 0.5s between each search so entrez server doesn\'t get angry')
num_accessions = len(accessions)
count_accessions = 1
for an in accessions:
    url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=' + an
    request = requests.get(url)
    # TODO: implement error handling for bad url
    split = request.text.split('\n', 1) # split on fist newline
    sequences[an] = split[1].replace('\n', '')

    for motif in motifs:
        for hit in motif_to_regex(motif).finditer(sequences[an]):
            matches.append([an, motif, hit, url]) # TODO: replace 'matches' with object

    print("[%s/%s] %s "%(count_accessions, num_accessions, an))
    time.sleep(0.5)
    count_accessions += 1

print('done looking for hits')

print()
for match in matches:
    print("%s,%s,%s,%s,%s" % (match[0], match[1], match[2].span(), match[2].group(0), match[3]))


# http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=940373999
