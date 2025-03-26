#!/bin/bash
PROJECT=${HOME}/study/bioinfo3
DATA=${PROJECT}/data_linkedomics
WORKDIR=${PROJECT}/fig1

MAF=${DATA}/Mutation_results.maf # input


# 총 124 column, 136개 KRAS mutation 
awk -F"\t" -v OFS=$'\t' '
    BEGIN {print "case_id","mutation","KRAS_VAF"}

    !($16 in id) && FNR>1 {id[$16]=$16} # Save all case_id

    $1=="KRAS" { # Save KRAS mutations
    denom=$41+$42
    denom=($41+$42==0) ? 1 : ($41+$42)
    vaf=$42/denom
    split($37,var,".")
    # print $16,var[2],vaf
    kras[$16]=var[2]"\t"vaf
    next
    } 

    END {
    for (i in id) {
        if (i in kras) {
            print i,kras[i]
        } else {
            print i,"NA","NA"
        }
    }}' ${MAF} > ${WORKDIR}/data/1d-caseid_variant_vaf.tsv