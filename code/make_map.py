#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:40:49 2020

@author: grant
"""
from combiner import Marker, get_map, insert_marker
from more_itertools import grouper


class possiblePositions(object):
    """Container for the bootstrap.txt output."""

    def __init__(self, marker, line):
        self.marker = marker
        self.hits = dict()
        self.fail_rate = 0

        mappings = [*grouper(line, 3)]

        sum_total = sum((float(i[2]) for i in mappings if [0] != "*"))

        for tup in mappings:
            if tup[0] == "*":
                self.fail_rate = float(tup[2])
            else:
                chr_ = int(tup[0])
                bp = int(tup[1])
                rate = float(tup[2]) / sum_total
                self.hits[(chr_, bp)] = rate

    def choose_pos(self):
        if self.fail_rate > .333:
            return ("*", "*")
        else:
            options = sorted([*self.hits.items()], key=lambda x: x[1],
                             reverse=True)
            if len(options) == 1:
                return options[0][0]
            else:
                if (options[0][1] - options[1][1]) / options[1][1] > 0.18:
                    return options[0][0]
                else:
                    return ("*", "*")

F2_map = get_map()
positions_dict = dict()

chrs_dict = {i : [] for i in range(1,27)}
chrs_dict["*"] = []
used_set = set()

with open("bootstrap.txt", "r") as readfile:
    for line in readfile.readlines():
        line = line.strip().split("\t")
        marker = line.pop(0)
        positions_dict[marker] = possiblePositions(marker, line)

for k, v in positions_dict.items():
    chr_, pos_ = v.choose_pos()
    chrs_dict, used_set = insert_marker(chrs_dict,
                                        chr_,
                                        k,
                                        pos_,
                                        used_set)

with open("bsplink.map", "w") as outfile:
    for chromosome, inserted in chrs_dict.items():
        for pos_, marker in inserted:
            field1 = chromosome
            field2 = marker
            field3 = F2_map[marker]
            field3 = "*" if field3[0] is None or field3[0] \
                != chromosome else field3[1]
            field4 = pos_
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(str(field1),
                                                        field2,
                                                        str(field3),
                                                        str(field4)))
