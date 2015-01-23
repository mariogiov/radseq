#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys

from Bio import SeqIO

# This object represents the location of the cut site in the sequence.
# It's basically just a number (int) with the added enzyme attribute to track
# which enzyme made the cut.
class Position(int):
    def __new__(cls, position, enzyme):
        obj = super(Position, cls).__new__(cls, position)
        obj.enzyme = enzyme
        return obj


def double_digest(input_fasta, lower_bound, upper_bound, enzyme_names):
    print("Checking desired enzymes... ", file=sys.stderr, end="")
    # This looks a little tricky but basically I'm just loading the
    # enzyme objects using the names that are supplied by the user.
    # This is the same as doing something like:
    #   from Bio.Restriction import EcoRI
    # but done dynamically.
    enzymes = __import__("Bio.Restriction", fromlist=enzyme_names[:2])
    try:
        enz_one = getattr(enzymes, enzyme_names[0])
    except AttributeError:
        raise ValueError("No such enzyme available: '{}'".format(enzyme_names[0]))
    try:
        enz_two = getattr(enzymes, enzyme_names[1])
    except AttributeError:
        raise ValueError("No such enzyme available: '{}'".format(enzyme_names[1]))
    print("successfully loaded enzmyes '{}' and '{}'.".format(enz_one, enz_two), file=sys.stderr)
    # Here I combine the enzymes to create a RestrictionBatch object, which
    # I can use to perform the in silico double-digest using the 'search' method
    r_batch = enz_one + enz_two

    print("Reading input from fasta file '{}'...".format(input_fasta), file=sys.stderr)
    with open(input_fasta, 'rU') as f:
        seq_records = list(SeqIO.parse(f, "fasta"))
    print("Calculating length of input fasta sequences...", file=sys.stderr)
    total_length_genome = sum(map(len, seq_records))
    print("Total length fasta sequences: {:,} bp".format(total_length_genome), file=sys.stderr)

    # Here we go through all the sequences loaded from the fasta file and
    # see where the enzymes would produce cuts. I then store these positions
    # as Position objects (from class above) so I can sort them by position
    # but keep the information about which enzyme made the cut
    frag_lens_dig_by_both = []
    for seq_record in seq_records:
        digest_dict = r_batch.search(seq_record.seq)
        pos_by_enzyme = []
        for enzyme_name, pos_list in digest_dict.iteritems():
            pos_by_enzyme.extend( [ Position(x, enzyme_name) for x in pos_list ] )
        pos_by_enzyme_iter = iter(sorted(pos_by_enzyme))

        try:
            next_elt = pos_by_enzyme_iter.next()
        except StopIteration:
            # No sequences digested
            continue
        while True:
            try:
                cur_elt = next_elt
                next_elt = pos_by_enzyme_iter.next()
                if cur_elt.enzyme != next_elt.enzyme:
                    # If two adjacent positions are cut by different enzymes,
                    # it produces a fragment having both cut sites (which is
                    # what we want); the difference in the positions is equal
                    # to the size of the fragment, so we record this.
                    frag_lens_dig_by_both.append(next_elt - cur_elt)
            except StopIteration:
                # End of the list of cut sites
                break

    # Here we restrict this list to fragments in our size range of interest
    fragments_in_range = [ length for length in frag_lens_dig_by_both \
                           if length >= lower_bound and length <= upper_bound ]
    length_fragments_in_range = sum(fragments_in_range)
    fraction_genome_repr = length_fragments_in_range / float(total_length_genome)

    print("Number of fragments in range {:,}-{:,}: {:,}".format(lower_bound,
                                                          upper_bound,
                                                          len(fragments_in_range)),
                                                   file=sys.stderr)
    print("Fraction of genome represented: {:.5f} ({:.2f}%) "
          "({:,}/{:,} bp)".format(fraction_genome_repr, fraction_genome_repr * 100,
                              length_fragments_in_range, total_length_genome),
                       file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Perform a double digest of a fasta sequence "
                                     "and determine both what percentage of fragments "
                                     "fall into the target range as well as what "
                                     "percentage of fragments in this range will "
                                     "represent usable reads. Requires Biopieces "
                                     "and thus all its dependencies.")
    parser.add_argument("-i", "--input-fasta", required=True,
            help=("The path to the input fasta file containing the sequence "
                  "to be digested."))
    parser.add_argument("-l", "--lower-bound", type=int, required=True,
            help=("The lower bound (in bp) of the size range of interest."))
    parser.add_argument("-u", "--upper-bound", type=int, required=True,
            help=("The upper bound (in bp) of the size range of interest."))
    parser.add_argument("enzyme_names", nargs=2,
            help=("The names of the two enzymes to use. Use standard naming for "
                  "example EcoRI (not ecori or EcoR1). Exactly two enzymes are "
                  "required."))

    args = vars(parser.parse_args())
    double_digest(**args)
