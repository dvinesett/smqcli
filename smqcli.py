#!/usr/bin/env python3
import re
import requests
import time

def motif_to_regex(raw_motif):
    motif = re.sub('[Xx]', '.', raw_motif)
    return re.compile(motif, re.IGNORECASE)

motifs = ['VXV','VXXV'] # arg
accessions = ['ALJ99974'] # arg
sequences = {}
matches = []
#sequence = 'aaattatagggatatata'

# TODO: validate motif

#matches = {re.compile(motif): None for motif in motifs} # initialize keys for dict with empty value
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

    print("[%s/%s] ... "%(count_accessions, num_accessions))
    time.sleep(0.5)
print('done querying web')

for an in sequences:
    #for match in list(matches):#.finditer(sequence): # list(matches) are the keys of dict
    for motif in motifs:
        for hit in motif_to_regex(motif).finditer(sequences[an]):
            #print("%s:\t%s\t%s" % (an, hit.span(), hit.group(0)))
            matches.append([an, motif, hit])
            #matches[match.span()] = match.group(0)
print('done looking for hits')

print()
for match in matches:
    print("%s\t%s\t%s\t%s" % (match[0], match[1], match[2].span(), match[2].group(0)))

#for match in match.finditer(sequence):
#    matches[match.span()] = match.group(0)

# http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=940373999
