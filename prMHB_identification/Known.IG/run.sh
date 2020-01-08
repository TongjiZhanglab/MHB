#!/bin/bash

sort -k1,1 -k2,2n known.imprinting.gene.tss | uniq > known.imprinting.gene.sort.tss
