#! /usr/local/bin/env python
# -*- coding: utf-8 -*-
import gzip
import sys
import os
import numpy as np
from scipy import stats

# patient infomation
p_list = "/Volumes/areca42TB/tcga/all_patient/gwas_focal_patient_list.tsv"
if(not os.path.exists(p_list)):
    print("ERROR::not exist "+p_list+"\n")
    sys.exit()
p_col = {}
p_inf = {}
with open(p_list) as pl:
    col = pl.readline().rstrip().split("\t")
    for col_num, c in enumerate(col):
        p_col[c] = col_num
    # main table row
    for l in pl:
        line = l.rstrip().split("\t")
        all_v, race = int(line[p_col["all_variant"]]), line[p_col["race"]]
        p_inf[line[p_col["patient_id"]]] = {}
        p_inf[line[p_col["patient_id"]]]["vnum"] = all_v
        p_inf[line[p_col["patient_id"]]]["race"] = race


# return vnum class array
def vnum_class_add(vnum, race="white"):
    if(vnum > 14 and race == "white"):
        return(np.array([0, 0, 0, 1]))
    elif(vnum > 12 and race == "white"):
        return(np.array([0, 0, 1, 1]))
    elif(vnum > 10 and race == "white"):
        return(np.array([0, 1, 1, 1]))
    elif(vnum > 8 and race == "white"):
        return(np.array([1, 1, 1, 1]))
    else:
        return(np.array([0]*4))


# return AN and AC ,,,
def allele_count(genotype, vnum, race="white"):
    if(genotype == -1 or race != "white"):
        return(np.array([0]*8))
    elif(genotype == 0):
        return(np.array([2, 0, 0, 0, 0, 0, 0, 0]))
    elif(genotype == 1):
        arr = np.array([2, 1, 1, 0])
        return(np.hstack((arr, vnum_class_add(vnum))))
    elif(genotype == 2):
        arr = np.array([2, 2, 0, 1])
        return(np.hstack((arr, vnum_class_add(vnum)*2)))
    else:
        print("ERROR::what genotype ??"+genotype+"\n")
        sys.exit()


def fisher_test(ac_p, ac, an_p, an):
    # print(ac_p, ac, an_p, an)
    table = np.array([[ac_p, an_p-ac_p], [ac-ac_p, an-an_p-ac+ac_p]])
    adds, p = stats.fisher_exact(table)
    return(p)


def fisher_by_row(ac_arr, vnum_an):
    # print(ac_arr.tolist(), vnum_an.tolist())
    p_10 = fisher_test(ac_arr[4], ac_arr[1], vnum_an[0], ac_arr[0])
    p_5 = fisher_test(ac_arr[5], ac_arr[1], vnum_an[1], ac_arr[0])
    p_3 = fisher_test(ac_arr[6], ac_arr[1], vnum_an[2], ac_arr[0])
    p_1 = fisher_test(ac_arr[7], ac_arr[1], vnum_an[3], ac_arr[0])
    return(np.array([p_10, p_5, p_3, p_1]))


infile = "/Volumes/areca42TB/tcga/array_genotype/"
# infile += "test.tsv.gz"
infile += "gentype_for_control_rare_gwas.tsv.gz"
if(not os.path.exists(infile)):
    print("ERROR::not exist "+infile+"\n")
    sys.exit()
outf = "/Volumes/areca42TB/tcga/array_genotype/gwas_all_patient.tsv.gz"
outwf = "/Volumes/areca42TB/tcga/array_genotype/gwas_white_patient.tsv.gz"
out = gzip.open(outf, mode='wt')
out.write("posi\tan\tac\thet\thom\tac_10\tac_5\tac_3\tac_1\t")
out.write("p_10\tp_5\tp_3\tp_1\n")
outw = gzip.open(outwf, mode='wt')
out.write("posi\tan\tac\thet\thom\tac_10\tac_5\tac_3\tac_1\t")
out.write("p_10\tp_5\tp_3\tp_1\n")
with gzip.open(infile, mode="rt") as genoty:
    coln = {}
    col = genoty.readline().rstrip().split("\t")
    col.pop(0)
    for cn, c in enumerate(col):
        coln[cn] = c

    # main infile line
    ln = 0
    for l in genoty:
        ln += 1
        line = l.rstrip().split("\t")
        vnum_class_an = np.array([0]*4)  # 10%(>8), 5%(>10), 3%(>12), 1%(>14)
        vnum_class_anw = np.array([0]*4)
        ac_all = np.array([0]*8)  # an, ac, ac_het, ac_hom, 10%, 5%, 3%, 1%
        ac_w = np.array([0]*8)
        position = line.pop(0)
        for cn, ge in enumerate(line):
            vnum, race = p_inf[coln[cn]]["vnum"], p_inf[coln[cn]]["race"]
            if(int(ge) != -1):
                vnum_class_an += vnum_class_add(vnum)*2
                vnum_class_anw += vnum_class_add(vnum, race)*2
            ac_all += allele_count(int(ge), vnum)
            ac_w += allele_count(int(ge), vnum, race)
        p_value_all = fisher_by_row(ac_all, vnum_class_an)
        print(*ac_all.tolist(), *p_value_all.tolist(), sep='\t', file=out)
        p_value_w = fisher_by_row(ac_w, vnum_class_anw)
        print(*ac_w.tolist(), *p_value_w.tolist(), sep='\t', file=outw)

out.close()
outw.close()
