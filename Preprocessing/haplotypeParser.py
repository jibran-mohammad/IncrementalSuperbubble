"""
haplotypeParser.py

This module parses:

1) ordered_haplotypes.txt
2) unique_elements.txt

Outputs:

hap_order:
    list of haplotypes in insertion order

hap_edges:
    dictionary:
        haplotype → list of edges

This module is independent of any superbubble tool.
"""


def read_haplotype_order(file):

    order=[]

    with open(file) as f:

        for line in f:

            line=line.strip()

            if line.startswith("# HAPLOTYPE"):

                hap=line.split(":")[1].strip()

                order.append(hap)

    return order


def parse_unique_elements(file):

    hap_edges={}

    current=None

    with open(file) as f:

        for line in f:

            line=line.strip()

            if line.startswith("# HAPLOTYPE"):

                current=line.split(":")[1].strip()

                hap_edges[current]=[]

            elif line.startswith("EDGE"):

                edge=line.split()[1]

                hap_edges[current].append(edge)

    return hap_edges