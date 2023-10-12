'''
This script was written in Python 3.9 and adapted from https://gist.github.com/philippbayer/0052f5ad56121cd2252a1c5b90154ed1.

It extracts the motif output from NLR Annotator (https://github.com/steuernb/NLR-Annotator) and converts it into CNL or TNL subfamily classification.

'''

import os
files = os.listdir(".")

# This motif table is from https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-75
motifs_d = {1:'nb_arc_cnl_or_tnl',2:'nb_arc_cnl',3:'nb_arc_cnl_or_tnl',4:'nb_arc_cnl_or_tnl',5:'nb_arc_cnl_or_tnl',6:'nb_arc_cnl',7:'linker_cnl_or_tnl',8:'linker_cnl_or_tnl',9:'lrr_cnl_or_tnl',10:'nb_arc_cnl_or_tnl',11:'lrr_cnl_or_tnl',12:'nb_arc_cnl_or_tnl',13:'tir_tnl',14:'monocot',15:'tir_tnl',16:'prenb_cnl',17:'prenb_cnl',18:'tir_tnl',19:'lrr_cnl_or_tnl',20:'monocot'}

class_dict = {frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'CNL',frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CNL',frozenset(['nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'prenb_cnl', 'nb_arc_cnl']):'CNL', frozenset(['monocot', 'nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCNL', frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'NL', frozenset(['nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN',frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']):'N', frozenset(['tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'NL', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CNL', frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'TCNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl']):'CN', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl']):'N', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl']):'NL', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']):'NL', frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl']):'TCNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']):'CN', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CNL'}

from collections import defaultdict

o = open("meta.tsv", "w")

o.write("taxa\tnlr\tstatus\tstart\tend\tstrand\tmotif\tspecies\tclass")

for file in files:
    if ".nlr.txt" in file:
        sp = file.removesuffix(".nlr.txt")
        
        with open(file, "r") as i:
            for line in i:
                ll = line.rstrip().split('\t')
                motifs = ll[-1].split(',')   # extracts the motif output from NLR Annotator
                this_domains = set()
                for m in motifs:   # iterates over the motifs
                    m  = int(m)
                    this_domains.add(motifs_d[m])
                this_domains = frozenset(this_domains)
                if this_domains in class_dict:    # classifying based on the motif configurations
                    this_class = class_dict[this_domains]
                else:    # otherwise return unclassified
                    this_class = 'Unclassified'
                
                o.write("\n"+line.rstrip()+"\t"+sp+"\t"+this_class)
