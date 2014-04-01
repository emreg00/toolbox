#!/usr/bin/env python
#
# Sequence related utilitiesa - written during development 
# of the test project of GFang
#
# Emre Guney, 20/12/2012, distributed under GNU-GPL.

from bisect import bisect_left, bisect_right

def read_fasta_file(fasta_file):
    f = open(fasta_file)
    for i, line in enumerate(f):
        line = line.strip()
        if line == "":
            continue
        # Only consider fasta header lines
        if line[0] == ">":
            gene, positions, strand = parse_fasta_header(line)
            # For this example inter-genic regions are skipped and the gene is assumed to span whole region
            start = positions[0][0]
            end = positions[-1][1]
    f.close()
    return

def get_closest_distances(positions, start, end, k):
    """
    Returns the k closest distances to start and end within positions array.

    >>> get_closest_distances([ 3, 6, 12, 14 ], 7, 9, 2)
    [-1, 3]
    >>> get_closest_distances([ 3, 6, 8, 9, 10, 12, 14 ], 7, 9, 3)
    [0, 0, -1]
    """
    distances = []
    zeros = []
    # Trick to only use sublist (positions are sorted)
    lower = bisect_left(positions, start)
    upper = bisect_right(positions, end)
    lower -= 10
    if lower < 0: lower = 0
    upper += 10
    positions = positions[lower:upper]
    for pos in positions:
        if pos <= start:
            distances.append(pos-start)
        elif pos > start and pos <= end:
            zeros.append(0)
        elif pos > end:
            distances.append(pos-end)
        else:
            raise ValueError("Logical error")
    distances.sort(lambda x,y: cmp(abs(x), abs(y)))
    zeros.extend(distances)
    return zeros[:k]

def parse_fasta_header(line):
    """
    Returns gene_name, [(start, end), ..], strand for a given fasta header line.

    >>> parse_fasta_header(">lcl|NC_000913.2_cdsid_NP_417358.2 [gene=xanQ] [protein=xanthine permease] [protein_id=NP_417358.2] [location=3022373..3023773]")
    ('xanQ', [(3022373, 3023773)], '+')
    
    >>> parse_fasta_header(">lcl|NC_000913.2_cdsid_NP_414616.1 [gene=leuA] [protein=2-isopropylmalate synthase] [protein_id=NP_414616.1] [location=complement(81958..83529)]")
    ('leuA', [(81958, 83529)], '-')

    >>> parse_fasta_header(">lcl|NC_000913.2_cdsid_NP_417367.1 [gene=prfB] [protein=peptide chain release factor RF-2] [protein_id=NP_417367.1] [location=complement(join(3033206..3034228,3034230..3034304))]")
    ('prfB', [(3033206, 3034228), (3034230, 3034304)], '-')
    """
    import re

    # Regular expressions to match id and location
    #exp_id = re.compile("\[protein_id=([a-zA-Z0-9_\.]+)\]")
    exp_id = re.compile("\[gene=([a-zA-Z0-9]+)\]")
    exp_loc = re.compile("\[location=([a-zA-Z0-9_\.(),]+)\]")

    positions = []
    strand = '+'
    protein_id = None
    m = exp_id.search(line)
    if m:
        protein_id = m.group(1)
    start, end = None, None
    m = exp_loc.search(line)
    if m:
        loc_str = m.group(1)
        if loc_str.startswith("complement"):
            strand = '-'
            loc_str = loc_str[11:-1]
        if loc_str.startswith("join"):
            loc_str = loc_str[5:-1]
            for pair in loc_str.split(","):
                start, end = map(int, pair.split(".."))
                positions.append((start, end))
        else:
            start, end = map(int, loc_str.split(".."))
            positions.append((start, end))
    return protein_id, positions, strand

#def flip(c):
#    d = { 'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A' }
#    return d[c]

def find_all_occurrences_of_motif(motif, fasta_file):
    """
    Returns all occurrences of the motif in a fasta file.

    >>> find_all_occurrences_of_motif("GATC", "mock_genome.fasta")
    [61, 68, 138]
    """
    f = open(fasta_file)
    # Skip fasta header
    f.readline()
    # For remainder motifs coming from previous line
    previous_line = "XXX"
    all_occurrences = []
    for i, line in enumerate(f):
        line = line.strip()
        if line == "":
            continue
        # All indices will be shifted by 3 due to consideration of remainders
        sequence = previous_line[-3:] + line
        previous_line = line
        # Assumes properly formatted FASTA
        if len(line) != 70:
            print "Line %d in FASTA file does not contain 70 characters" % (i+2)
        base_position = -3 + i*70
        occurrences = search_motif(motif, sequence)
        for k in occurrences:
            all_occurrences.append(base_position+k)
    f.close()
    return all_occurrences

def search_motif(motif, sequence):
    """
    Returns occurrences (start indices) of the motif in the sequence.

    >>> search_motif("GATC", "ACAACGATCGCCAGCAGATCGTGAGGC")
    [5, 16]
    >>> search_motif("GATC", "GATCGCCAGCAGATC")
    [0, 11]
    >>> search_motif("GATC", "GCTCGCCAGCATATC")
    []
    >>> search_motif("GATC", "GAT")
    []
    """
    occurrences = []
    i_limit = len(sequence) - 4
    for i, c in enumerate(sequence):
        if c == motif[0]:
            if i > i_limit:
                continue
            if sequence[i:i+4] == motif:
                occurrences.append(i)
    return occurrences


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main()


