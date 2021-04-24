#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:14:26 2020

@author: grant
"""
import sys
from bisect import insort as insort
from bisect import bisect_left as bisect
from random import randrange, choice, sample
from itertools import product
from math import log10, sqrt
from datetime import datetime

# files needed:
#   BLAST output
#   marker correlations from PLINK
#   F2 map data

# let's start with a single chromosome

#%% class definitions


class Marker(object):
    """
    A container for markers, ones they are mapped.
    """

    def __init__(self, name, chr_, pos_):
        """
        Takes three arguments:
            * name : the name of the marker
            * chr_ : the chromosome number, or "*", or None
            * pos_ : the positon on the chromosome, or "*", or None
        """
        self.name = name

        if chr_ is None or chr_ == "*":
            self.chr_ = "*"

        if pos_ is None or pos_ == "*":
            self.pos_ = "*"

        self.chr_ = chr_
        self.pos_ = pos_

    def set_F2(self, map_chr, map_pos):
        if map_chr is None or map_pos is none:
            self.map_chr_ = "*"
            self.map_pos_ = "*"
        else:
            self.map_chr_ = map_chr
            self.map_pos_ = map_pos

        """if self.map_chr !- self.chr_:
            self.map_pos_ = "*"
        else:
            self.map_pos_ = map_pos"""

    def plinkmap(self):
        """
        Convert to a plink map string.
        """
        return "{0}\t{1}\t{2}\t{3}".format(self.chr_,
                                           self.marker,
                                           self.map_pos_,
                                           self.pos_)


    def __str__(self):
        """Convert to str."""
        return "{0}\t{1}\t{2}".format(self.name,
                                      self.chr_,
                                      self.pos_)
    __repr__ = __str__


#%% function defs


def get_map():
    """
    Reads in the F2 map file, to get an idea of where to start with the
    SNP loci.
    """
    print("Reading in intraspecific F2 map.")
    print()

    with open("newplink.map") as fp:
        F2_map = {line[1]: tuple([int(line[0]), float(line[2])])
                  if line[0].isnumeric() else tuple([None, None])
                  for line in (entry.strip().split("\t")
                               for entry
                               in fp.readlines())}

    print("\tTotal number of markers read in: {0}".format(len(F2_map)))
    print("\tNumber of markers for each chr: ")

    for i in range(1, 27):
        print("\t{0}\t{1:>5}".format(str(i), sum((1 for v in F2_map.values()
                                                  if v[0] == i))))
    else:
        print("\t{0}\t{1:>5}".format("n/a", sum((1 for v in F2_map.values()
                                                 if v[0] is None))))
    print()

    return F2_map


def get_corr(loci_names, cutoff=0.2):
    """
    Reads in the correlations files. Takes a single argument:
        * loci names: an iterable containing the loci names to add to dict
    Returns a marker : {marker: corr} dict
    """
    print("Reading in marker correlations.")

    corr_dict = {i: dict() for i in loci_names}
    with open("sorted.ld") as fp:
        for line in fp.readlines():
            line = line.strip().split()
            marker1 = line[0]
            marker2 = line[1]
            corr = float(line[2])
            if corr < cutoff:
                continue

            if marker1 in loci_names and marker2 in loci_names:
                corr_dict[marker1][marker2] = corr

    print("\tNumber of markers with no correlated markers (R2 >= {0}): {1}"
          .format(cutoff, sum((1 for v in corr_dict.values() if not v))))

    return corr_dict


def get_blast(loci_names, filename):
    """
    Reads in the BLASTn result file. Takes a single argument:
        * loci names: an iterable containing the loci names to add to dict
    """
    print("Reading in BLASTn hits.")

    # maps new D chromosome names to old ones
    D_dict = {1: 15,
              2: 14,
              3: 17,
              4: 22,
              5: 19,
              6: 25,
              7: 16,
              8: 24,
              9: 23,
              10: 20,
              11: 21,
              12: 26,
              13: 18}

    # initialize the empty dictionary using the marker names from F2_map
    blast_dict = {i: [] for i in loci_names}
    with open(filename) as fp:
        for line in fp.readlines():
            line = line.strip().split()
            marker = line[0].split("_")[0]

            if marker not in blast_dict:
                continue

            # skip if its anchored to a scaffold
            if line[1].startswith("scaffold"):
                continue
            # convert to 26 chromosome notation
            chromosome = int(line[1][1:]) if line[1][0] == "A" \
                else D_dict[int(line[1][1:])]
            if line[8] < line[9]:
                # deal with the case where its on the positive strand
                start_pos = int(line[8]) - int(line[6]) + 1
            else:
                # deal with the case where the match is on the negative strand
                start_pos = int(line[8]) + int(line[6]) - 1
            # start_pos = int(line[8])
            len_match = int(line[3])
            escore = float(line[10])
            add_tup = tuple([chromosome, start_pos, len_match, escore])
            # only add if not in lst
            if not len([i for i in blast_dict[marker]
                        if i[0] == chromosome and i[1] == start_pos]):
                blast_dict[marker].append(add_tup)
            elif [i for i in blast_dict[marker] if i[0] == chromosome and
                  i[1] == start_pos][0][2] < len_match:
                # do add if they start at the same place, but this match
                # is longer
                blast_dict[marker].remove([i for i in blast_dict[marker]
                                           if i[0] == chromosome and
                                           i[1] == start_pos][0])
                blast_dict[marker].append(add_tup)
            else:
                # do not add
                pass

        # print some summary statistics from BLAST
        print("\tresults from BLASTn:")
        print("\tmultiple hits: ", sum((1 for i in blast_dict.values()
                                        if len(i) > 1)))
        print("\tone hit:       ", sum((1 for i in blast_dict.values()
                                        if len(i) == 1)))
        print("\tzero hits:     ", sum((1 for i in blast_dict.values()
                                        if len(i) == 0)))
        print()

        return(blast_dict)


def find_RBM(loci_names, corr, cutoff=0.2):
    """
    Function to find reciprocal best matches in the correlation file. Args:
        loci_names: iterable with the loci names
        corr: the dictionary with the locus : (marker, corr) dictionary
    """
    # create a dictionary to find the highest for each marker
    highest_corr = {i: None for i in loci_names}

    # loop through all the loci that have an R2 > cutoff for at least one pair
    for locus in loci_names:
        # skip loci with no correlations
        if not corr[locus]:
            continue
        # find the highest correlation
        max_corr = max(corr[locus].values())
        # find the loci with the max correlation
        max_loci = [i for i, j in corr[locus].items()
                    if j == max_corr and i in loci_names]
        # add the locus if theres only one that meets this criteria
        if len(max_loci) == 1 and max_corr > cutoff:
            highest_corr[locus] = max_loci[0]
        elif len(max_loci) > 1:
            pass
        else:
            # set the locus to none if there isn't a unique solution
            highest_corr[locus] = None

    # build the dictionary to output
    out_dict = {i: None for i in loci_names}

    # loop thorugh all the loci
    for locus in highest_corr:
        # get the second locus for the match-up
        locus2 = highest_corr[locus]
        # skip it if either one doesn't have a best match
        if (highest_corr[locus] is None) or (highest_corr[locus2] is None):
            continue
        # assuming they both have a match, see if they are equal
        elif locus == highest_corr[locus2]:
            out_dict[locus] = locus2
        else:
            # do nothing
            pass

    # returns the dictionary
    # contains every locus, with either a None or a marker name
    return out_dict


def choose_RBM_maps(RBM_map, F2_map, blast_dict, noties=True):
    """
    Function to find the best place to anchor the RBM pairs. Args:
        * RBM_map : the output from find_RBM, filtered to only include RBMs
        * F2_map : output dict from get_map
        * blast_dict : output from get_blast
        * noties: disables the tiebreaker
    The way the function is written now, tie breaks won't work iff
    the distances between markers are the same in multiple places.
    print("\tFinding the best BLAST hit for these matches, {0} tiebreaks."
          .format("without" if noties else "with"))    """

    # build a set of RBMs, then convert to lst
    RBM_set = set()
    for tup in RBM_map.items():
        RBM_set.add(tuple(sorted(tup)))

    RBM_lst = [sorted(i) for i in RBM_set]
    RBM_matches = [(F2_map[i][0], F2_map[j][0]) for i, j in RBM_lst]
    # dont add in a match if the RBMs are on different places in the F2 map
    rm_lst = []
    same_lst = []
    nochr_lst = []
    for n in range(len(RBM_matches)):
        # leave it be if either one is unmapped
        if not (RBM_matches[n][0] and RBM_matches[n][1]):
            nochr_lst.append(RBM_lst[n])
        elif RBM_matches[n][0] != RBM_matches[n][1]:
            rm_lst.append(RBM_lst[n])
        # the ones mapped to the same chr are the ones we want to work with
        elif RBM_matches[n][0] == RBM_matches[n][1]:
            same_lst.append(RBM_lst[n])

    # remove the offending matches
    for i in rm_lst:
        RBM_lst.remove(i)

    out_dict = dict()

    # find the closest two blast hits for each pair.
    # first, do the matches known to be on the same chromosome
    difference = {}

    for marker1, marker2 in same_lst:
        hits1 = filter_blast(marker1, blast_dict, F2_map[marker1][0], None,
                             n=20)
        hits2 = filter_blast(marker2, blast_dict, F2_map[marker2][0], None,
                             n=20)
        # do nothing if there are no blast hits
        if not hits1 or not hits2:
            continue

        difference = {abs(i[1][1] - j[1][1]): (i[0], j[0])
                      for i, j in product(enumerate(hits1), enumerate(hits2))}

        if noties and len(difference) > 1:
            continue

        idx1, idx2 = difference[min(difference)]

        out_dict[marker1] = hits1[idx1]
        out_dict[marker2] = hits2[idx2]

    for marker1, marker2 in nochr_lst:
        chr_ = F2_map[marker1][0] or F2_map[marker2][0]
        # if one of them is on the F2 map, depend on that chr
        if chr_ is not None:
            hits1 = filter_blast(marker1, blast_dict,
                                 chr_, None,
                                 n=20)
            hits2 = filter_blast(marker2, blast_dict,
                                 chr_, None,
                                 n=20)

            if not hits1 or not hits2:
                continue

            difference = {abs(i[1][1] - j[1][1]): (i[0], j[0])
                          for i, j in product(enumerate(hits1),
                                              enumerate(hits2))}
            if noties and len(difference) > 1:
                continue

            idx1, idx2 = difference[min(difference)]
            out_dict[marker1] = hits1[idx1]
            out_dict[marker2] = hits2[idx2]
        # figure out which chromosome it should be on
        else:
            if noties:
                continue
            # try each chromosome separately
            trial_dict = {i: None for i in range(1, 27)}
            for num in trial_dict.keys():
                hits1 = filter_blast(marker1, blast_dict,
                                     num, None,
                                     n=20)
                hits2 = filter_blast(marker2, blast_dict,
                                     num, None,
                                     n=20)

                if not hits1 or not hits2:
                    trial_dict[num] = (None, None, None)
                    continue

                difference = {abs(i[1][1] - j[1][1]): (i[0], j[0])
                              for i, j in product(enumerate(hits1),
                                                  enumerate(hits2))}

                idx1, idx2 = difference[min(difference)]
                trial_dict[num] = (min(difference), hits1[idx1], hits2[idx2])

            # pick the best chromosome
            good_lst = [i for i in trial_dict.values() if i[0] is not None]

            if not good_lst:
                continue
            else:
                best_choice = min(good_lst, key=lambda x: x[0])
                out_dict[marker1] = best_choice[1]
                out_dict[marker2] = best_choice[2]

    return out_dict


def find_neighbors(start_locus, correlations, excluded={}, cutoff=0.8):
    """
    Get the neighbors for a given locus. Takes two arguments:
        start_locus : the locus to lookup in the correlations dict
        correlations : the dictionary with the correlations
        excluded : skips these loci if theyve already been accounted for
    Returns a list of tuples of the form (locus, R2), sorted on R2.
    """
    return [i for i in sorted(correlations[start_locus].items(),
                              key=lambda x: x[1],
                              reverse=True)
            if i[0] not in excluded and i[1] >= cutoff]


def filter_blast(locus, dict_, chr_, pos_=None, n=20):
    """
    Gets the n highest entries for a blast search of that probe seq. Args:
        * locus : the locus to look up in the blast dict
        * dict : the dictionary mapping locus name to blast hits on chr
        * chr_ : the chr to restrict entries to, can be None
        * n : the number of blast hits to return
    Returns a list of blast hits, or an empty list, but usually a single hit.
    """
    out_lst = []

    # control flow to handle if chr_ is None
    if chr_ is None:
        for entry in dict_[locus]:
            out_lst.append(entry)

        return out_lst

    for entry in dict_[locus]:
        if entry[0] == chr_:
            out_lst.append(entry)

    # return the empty lst
    if not out_lst:
        return out_lst

    # find the lowest value
    min_val = min(out_lst, key=lambda x: x[2])[2]

    # determine if that value is unique, for the e-value
    min_vals = [i for i in out_lst if i[2] == min_val]

    if n == 1:
        # if there is exactly one match on that chromosome, use it!
        if len(min_vals) == 1:
            return [min_vals[0]]
        elif len(min_vals) > 1:
            if pos_:
                # select the anchor pos closest to this marker
                [sorted(out_lst, key=lambda x: abs(x[1] - pos_))[0]]
            else:
                # otherwise, just pick a random marker
                return [choice(out_lst)]
        # if there is more than one or none, return an empty lst
        else:
            return []
    else:
        if len(min_vals) == 1:
            return out_lst
        elif len(min_vals) > 1:
            if pos_:
                # sort the markers by distance to pos, but dont select one
                sorted(out_lst, key=lambda x: abs(x[1] - pos_))
            else:
                # return the list, sorted by chr position
                return sorted(out_lst, key=lambda x: x[1])[:20]
        # if there is more than one or none, return an empty lst
        else:
            return []


def insert_marker(chrs_dict, chr_, marker_name, pos, used_set):
    """
    Uses bisect to insert a tuple of the following form, (pos, name).
    Has the following args:
        * chrs_dict: a dictionary with each chr_ as a key
        * chr_: the chromosome being inserted into
        * marker_name: the name of the probe locus
        * pos: the bp coordinate of the marker
    Returns the chrs_dict, and used_set.
    """
    if marker_name in used_set:
        print(marker_name, "already inserted")
        return chrs_dict, used_set

    insort(chrs_dict[chr_], tuple([pos, marker_name]))
    used_set.add(marker_name)

    return chrs_dict, used_set


def build_chr(corr_dict, blast_dict, marker_dict,
              anchor, anchor_chr, anchor_pos,
              used_markers,
              cutoff=0.8):
    # better option: build a list of newly added markers
    # get index for potential marker starting points
    # options = [i for i in range(len(chrs_dict[anchor_chr]))]

    # loop until random choices no longer add new markers
    # while options:
    # remove that options from the current pool
    # current_marker = (chrs_dict[anchor_chr]
    #                  [options.pop(randrange(len(options)))][1])

    # find neighbors based on correlation
    neighbors = find_neighbors(anchor, corr_dict, used_markers, cutoff)

    # find eligible blast hits
    neighbors_blast = {neighbor[0]: filter_blast(neighbor[0],
                                                 blast_dict,
                                                 anchor_chr,
                                                 anchor_pos,
                                                 n=20)
                       for neighbor in neighbors}

    # filter to only include succesful searches
    neighbors_blast = {i: j for i, j in neighbors_blast.items() if j}

    # loop through each of the neighboring loci
    for i in neighbors_blast:
        # if there is a single hit  leftover (!!), attempt to add it
        if len(neighbors_blast[i]) == 1:
            marker_dict, \
                used_markers = insert_marker(marker_dict,
                                             neighbors_blast[i][0][0],
                                             i,
                                             neighbors_blast[i][0][1],
                                             used_markers)
        else:
            pass

    return marker_dict, used_markers


#%% main section


def build_map(F2_map_, blast_dict_, cutoff_, blast_dict_alt_):
    # first, read in the F2 map, the whole map (filter later)
    F2_map = F2_map_

    # read in the BLAST results
    blast_dict = blast_dict_

    # set the R2 cutoff to allow to be imported
    cutoff = cutoff_
    corr_dict = get_corr(F2_map, cutoff)

    print()

    # get the best matches
    best_matches = find_RBM(F2_map, corr_dict, cutoff_)
    best_matches = {i: j for i, j in best_matches.items() if j is not None}
    print("\t{0} sets of best matches".format(int(len(best_matches) / 2)))

    print()

    best_hits = choose_RBM_maps(best_matches, F2_map, blast_dict)

    # build the marker dictionary, to be filled in sorted order
    marker_dict = {i: [] for i in range(1, 27)}
    rev_marker_dict = dict()
    used_markers = set()

    print()

    # 1 : INSERT RBMs first
    # go through and insert the eligible best matches
    for marker, blast_results in best_hits.items():
        marker_dict, used_markers = insert_marker(marker_dict,
                                                  blast_results[0],
                                                  marker,
                                                  blast_results[1],
                                                  used_markers)

    print("Mapped {0} markers from best matches.".format(len(used_markers)))
    print()

    # maybe, we should try extending from other RBMs first?
    # like, maybe do
    tobeparsed = sample([i for i in used_markers],
                        int(len(used_markers) / 5))
    print("Extending from {0} markers....".format(len(tobeparsed)))

    for marker in tobeparsed:
        anchor = marker
        anchor_chr = blast_results[0]
        anchor_pos = blast_results[1]
        marker_dict, used_markers = build_chr(corr_dict,
                                              blast_dict,
                                              marker_dict,
                                              marker,
                                              anchor_chr,
                                              anchor_pos,
                                              used_markers,
                                              0.9)

    print("After extending from some markers, total {0}."
          .format(len(used_markers)))
    print()

    print("Now mapping other markers from the F2_map.")
    print()

    # find the ones on the F2 map that havent been anchored yet
    F2_unmapped = [i for i, j in F2_map.items()
                   if j[0] is not None and i not in used_markers]

    # find the matching blast hits on the proper chromosome
    blast_hits = {marker: filter_blast(marker,
                                       blast_dict,
                                       F2_map[marker][0])
                  for marker in F2_unmapped}

    # if theres exactly one hit, go for it
    for marker, hits in blast_hits.items():
        if len(hits) == 1:
            marker_dict, used_markers = insert_marker(marker_dict,
                                                      hits[0][0],
                                                      marker,
                                                      hits[0][1],
                                                      used_markers)

    print("Used F2 map to finish up more markers, total {0}."
          .format(len(used_markers)))

    print()

    # start with just the F2 mapped markers
    tobedone_set = {i for i, j in F2_map.items() if j[0] is not None} \
        - used_markers

    print("Number unanchored, but in F2 map with entry: {0}"
          .format(len(tobedone_set)))

    # 3 : INSERT
    # loop through the other markers with known chromosome...
    # this time, by using the distance to anchored markers

    while tobedone_set:
        # get a random element from this set
        anchor = choice([i for i in tobedone_set])
        tobedone_set.remove(anchor)
        anchor_chr = F2_map[anchor][0]

        # get its BLAST matches
        anchor_blast = filter_blast(anchor, blast_dict, anchor_chr)

        # if only one blast hit, go ahead and insert it
        if len(anchor_blast) == 1:
            chosen_score = anchor_blast[0]
            marker_dict, used_markers = insert_marker(marker_dict,
                                                      chosen_score[0],
                                                      anchor,
                                                      chosen_score[1],
                                                      used_markers)
            continue
        elif len(anchor_blast) == 0:
            # continue if there are no BLAST hits
            continue

        # deal with the case that there are multiple hits
        chr_bp = [i for i, j in marker_dict[anchor_chr]]
        chr_corr = [corr_dict[anchor][j] if j in corr_dict[anchor] else 0
                    for i, j in marker_dict[anchor_chr]]

        # calculate a score for each insertion point, using correlations
        scores = list()
        for result in anchor_blast:
            # get the bp position for this blast hit
            bp = result[1]

            # find the point +/- 5 Mbp
            before_bp = bp - 5000000
            before_bp = before_bp if before_bp > 0 else 0

            after_bp = bp + 5000000
            after_bp = after_bp if after_bp < chr_bp[-1] else chr_bp[-1]

            before_idx = bisect(chr_bp, before_bp)
            after_idx = bisect(chr_bp, after_bp)

            # calculate a score for this insertion point
            bp_vec = [1 / sqrt(abs(i - bp) + .00001)
                      for i in chr_bp[before_idx:after_idx]]

            bp_vec = sum([i * j for i, j in
                         zip(bp_vec, chr_corr[before_idx:after_idx])])

            # append it to the scoring matrix
            scores.append(tuple([bp_vec, result]))

        scores = sorted(scores, key=lambda x: x[0], reverse=True)
        # tiebreakers:
        # evalue tiebreaker
        # retain if evalue within 2
        # re-sort by e-value

        scores_e = sorted(scores, key=lambda x: x[1][3], reverse=False)

        eval_filt = [i for i in scores_e
                     if (-log10(scores_e[0][1][3]) - -log10(i[1][3])) < 2]

        # do not add if the metric gives two best places to insert

        if eval_filt[0] != scores[0]:
            continue

        # if e-val filter only contains one element, insert that one
        if len(eval_filt) == 1:
            chosen_score = eval_filt.pop(0)[1]
        elif len(eval_filt) > 1:
            chosen_score = scores[0][1]

        marker_dict, used_markers = insert_marker(marker_dict,
                                                  chosen_score[0],
                                                  anchor,
                                                  chosen_score[1],
                                                  used_markers)

    print()

    last_markers = {i for i in F2_map} - used_markers

    print("Working through remaining {0} markers, one by one..."
          .format(len(last_markers)))

    # 4 : INSERT
    # loop through the remaining markers, picking a chromosome based on
    # longest stretch of LD

    prev_size = len(last_markers)
    new_blast_dict = False
    tie_breaking = False

    # update the reverse marker dictionary
    for chr_, lst_ in marker_dict.items():
        for loc, marker in lst_:
            rev_marker_dict[marker] = Marker(marker, chr_, loc)

    while True:
        # set the flags and loop variables
        if not last_markers:
            last_markers = {i for i in F2_map} - used_markers
            if len(last_markers) == prev_size or not last_markers:
                # once all the good choices are taken up, start tie breaking!
                if not tie_breaking:
                    print("\t***enabling tie breaking!")
                    tie_breaking = True
                elif not new_blast_dict:
                    new_blast_dict = True
                    print("\t***reading in the bad blast hits....")
                    blast_dict = blast_dict_alt_
                else:
                    break
            else:
                prev_size = len(last_markers)

        # pick a marker to use
        anchor = choice([i for i in last_markers])
        last_markers.remove(anchor)

        # get all the blast hits
        anchor_blast = filter_blast(anchor, blast_dict, None)

        # pick a chromosome if its on the F2_map
        anchor_F2_chr = F2_map[anchor][0]

        # subset the blast hits if the chromosome is known from genetic map
        if anchor_F2_chr is not None:
            anchor_blast = [i for i in anchor_blast if i[0] == anchor_F2_chr]

        # if only one blast hit, go ahead and insert it
        if len(anchor_blast) == 1:
            chosen_score = anchor_blast[0]
            marker_dict, used_markers = insert_marker(marker_dict,
                                                      chosen_score[0],
                                                      anchor,
                                                      chosen_score[1],
                                                      used_markers)
            # XXX
            rev_marker_dict[anchor] = Marker(anchor, chosen_score[0],
                                             chosen_score[1])
            continue
        elif len(anchor_blast) == 0:
            # continue if there are no BLAST hits
            continue

        # deal with the case that there are multiple hits
        # calculate a score for each insertion point, using correlations
        scores = list()

        # number of correlated markers on each chromosome:
        ct_corrdict = {i: 0 for i in range(1, 27)}
        for marker, corr in corr_dict[anchor].items():
            if marker in used_markers:
                ct_corrdict[rev_marker_dict[marker].chr_] += 1

        for result in anchor_blast:
            # just use the term anchor_chr so dont have to replace all inst
            anchor_chr = result[0]
            bp = result[1]

            # get the data needed for scoring
            chr_bp = [i for i, j in marker_dict[anchor_chr]]
            chr_corr = [corr_dict[anchor][j] if j in corr_dict[anchor] else 0
                        for i, j in marker_dict[anchor_chr]]

            # find the point +/- 5 Mbp
            before_bp = bp - 5000000
            before_bp = before_bp if before_bp > 0 else 0

            after_bp = bp + 5000000
            after_bp = after_bp if after_bp < chr_bp[-1] else chr_bp[-1]

            before_idx = bisect(chr_bp, before_bp)
            after_idx = bisect(chr_bp, after_bp)

            # calculate a score for this insertion point
            # try with log10
            bp_vec = [1 / log10(abs(i - bp) + 10)
                      for i in chr_bp[before_idx:after_idx]]

            bp_vec = sum([i * j for i, j in
                         zip(bp_vec, chr_corr[before_idx:after_idx])])

            # append it to the scoring matrix
            scores.append(tuple([bp_vec, result]))

        # sort from highest to lowest score
        scores = sorted(scores, key=lambda x: x[0], reverse=True)

        # tiebreakers:
        # evalue tiebreaker
        # retain if evalue within 2
        # re-sort by e-value

        scores_e = sorted(scores, key=lambda x: x[1][3], reverse=False)

        eval_filt = [i for i in scores_e
                     if (-log10(scores_e[0][1][3]) - -log10(i[1][3])) < 2]

        # do not add if the metric gives two best places to insert
        # and tie breaking hasnt been triggered yet
        if (eval_filt[0] != scores[0]) and not tie_breaking:
            continue
        elif eval_filt[0] == scores[0]:
            chosen_score = eval_filt[0][1]
        # now actually deal with tie breakers
        # XXX might wanna mess with threshholds here
        elif tie_breaking:
            # tie break purely on marker correlation
            # this won't help map our markers with no friends....
            scores = [i for i in scores if i[0] > 0.01]
            if len(scores) >= 1:
                chosen_score = scores[0][1]
                # when this is done, automatically reset randomization,
                # try to map RBMs
                last_markers = {}
            else:
                # allow it to insert if the marker is not correlated
                if len(corr_dict[anchor]) == 0:
                    chosen_score = eval_filt[0][1]
                else:
                    continue

        marker_dict, used_markers = insert_marker(marker_dict,
                                                  chosen_score[0],
                                                  anchor,
                                                  chosen_score[1],
                                                  used_markers)

        rev_marker_dict[anchor] = Marker(anchor, chosen_score[0],
                                         chosen_score[1])

    print()
    print("total done ", len(used_markers))

    last_markers = {i for i in F2_map} - used_markers

    return marker_dict, last_markers


#%% bootstrapping loop

if __name__ == "__main__":
    pass
    start_time = datetime.now()

    # read in the constants
    F2_map_ = get_map()
    blast_dict_ = get_blast(F2_map_, "2alleles_all_more_filt.out")
    blast_dict_alt_ = get_blast(F2_map_, "blast_unmapped.txt")

    bootstrap_dict = {i: {("*", "*"): 0} for i in F2_map_}
    all_unmapped = set()

    reps = int(sys.argv[1])

    for x in range(0, reps):
        print()
        print("----------------------")
        print("starting rep {0}".format(x))
        print()

        cutoff = randrange(20, 80, 1) / 100
        marker_dict_, unmapped = build_map(F2_map_, blast_dict_, cutoff,
                                           blast_dict_alt_)
        all_unmapped |= unmapped
        marker_set = set([i for i in F2_map_])

        for chr_, full_map in marker_dict_.items():
            for place in full_map:
                marker = place[1]
                marker_set.remove(marker)
                pos = place[0]
                pos_tup = tuple([chr_, pos])

                if pos_tup not in bootstrap_dict[marker]:
                    bootstrap_dict[marker][pos_tup] = 1
                else:
                    bootstrap_dict[marker][pos_tup] += 1
        else:
            for marker in marker_set:
                bootstrap_dict[marker][("*", "*")] += 1
    else:
        for key, value in bootstrap_dict.items():
            for pos, count in value.items():
                bootstrap_dict[key][pos] = count / reps

    with open("bootstrap.txt", "w") as outfile:
        for marker in bootstrap_dict:
            outfile.write("{0}\t{1}\n"
                          .format(marker,
                                  "\t".join(["{0}\t{1}\t{2:.3f}"
                                             .format(str(k[0]),
                                                     str(k[1]),
                                                     v)
                                             for k, v
                                             in bootstrap_dict[marker]
                                             .items()])))
    end_time = datetime.now()
    print("Duration {0}".format(end_time - start_time))

    with open("unmapped.txt", "w") as outfile:
        for marker in all_unmapped:
            outfile.write(marker + "\n")


#%% helper code

"""
# helper for testing LD assumptions
print("Printing LD file...")
with open("ldcorr.txt", "w") as outfile:
    for marker1, marker2 in best_matches.items():
        if marker1 not in best_hits or marker2 not in best_hits:
            continue
        chr1 = best_hits[marker1][0]
        pos1 = best_hits[marker1][1]
        chr2 = best_hits[marker2][0]
        pos2 = best_hits[marker2][1]
        corr = corr_dict[marker1][marker2]
        mappos1 = F2_map[marker1][1] or "*"
        mappos2 = F2_map[marker2][1] or "*"
        outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n"
                  .format(marker1,
                          marker2,
                          chr1,
                          pos1,
                          chr2,
                          pos2,
                          corr,
                          mappos1,
                          mappos2))


# helper for writing out chr for new mappings
with open("F2_test/"+anchor+".txt", "w") as outfile:
    for option in anchor_blast:
        outfile.write("{0}\n"
                      .format("\t".join(str(i) for i in option)))
    for marker_, bp_, pos_, corr_ in zip(chr_markers, chr_bp,
                                         chr_pos, chr_corr):
        outfile.write("{0:.2f}\t{1:d}\t{2}\t{3:.4f}\n"
try:
    idx = int(input("please choose an index.... "))
    chosen_score = scores[idx][1]
except (ValueError, IndexError):
    print("not inserting, index out of range or char input")
    counter = 1
else:
    # do the insertion here i think
    marker_dict, used_markers = insert_marker(marker_dict,
                                              chosen_score[0],
                                              anchor,
                                              chosen_score[1],
                                              used_markers)
finally:
    pass

 with open("testfile.txt", "w") as outfile:
        for chr_ in marker_dict:
            for pos, marker in marker_dict[chr_]:
                if marker not in F2_map or F2_map[marker][0] is None:
                    map_chr = "*"
                    map_pos = "*"
                else:
                    map_chr = F2_map[marker][0]
                    map_pos = F2_map[marker][1]
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chr_,
                                                                 pos,
                                                                 map_chr,
                                                                 map_pos,
                                                                 marker))


if print_flag is True:
    print_flag = False
    print(anchor, file=fp)
    print(scores, file=fp)
    print(eval_filt, file=fp)
    print(F2_map[anchor], file=fp)
    print({i: j for i, j in ct_corrdict.items() if j > 0}, file=fp)
    print(corr_dict[anchor], file=fp)
    print([j for i, j in rev_marker_dict.items() if
           i in corr_dict[anchor]], file=fp)
    print("\n", file=fp)
    continue

best_matches = find_RBM(last_markers, corr_dict)
best_matches = {i: j for i, j in best_matches.items()
                if j is not None}

best_hits = choose_RBM_maps(best_matches, F2_map, blast_dict,
                            noties=False)

for marker, blast_results in best_hits.items():
    # insert the marker
    marker_dict, used_markers = insert_marker(marker_dict,
                                              blast_results[0],
                                              marker,
                                              blast_results[1],
                                              used_markers)
    # and remembers to add it into the reverse dict!
    rev_marker_dict[marker] = Marker(marker, blast_results[0],
                                     blast_results[1])

"""
