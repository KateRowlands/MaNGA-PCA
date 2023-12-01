#!/usr/bin/env python
# encoding: utf-8
#
# MaNGA_Utils.py
#
# Created by Kate Rowlands on 14 March 2017.

import numpy as np

def spx_map(binid):
    #Make a mask so we don't double count spaxels which are in a bin (and therefore have the same value)

    binids = np.sort(np.unique(binid))
    map = np.zeros(shape=(binid.shape[1],binid.shape[0]))

    #Store positions of (first spaxel of) binids which have one or more spaxels
    #First element of -1 is excluded
    for x in binids[1:]:
        #print(x)
        storey, storex = np.where(binid == x)

        #Pick a value near the middle of the bin
        id_x = int(len(storex)/2)
        id_y = int(len(storey)/2)

        #print(id_y, id_x)

        #For central spaxel
        middle_spx_x = storex[id_x]
        middle_spx_y = storey[id_y]
        #print(middle_spx_y, middle_spx_x)

        map[middle_spx_y, middle_spx_x] = 1

    return map
