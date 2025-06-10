#!/bin/bash
PROJECT=${HOME}/study/bioinfo3
DATA=${PROJECT}/data_linkedomics
WORKDIR=${PROJECT}/fig1

MAF=${DATA}/Mutation_results.maf # input

if true; then
    # 1. WES KRAS VAF
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



#2. WGS CNV
# For the mutations of each case, ($52 in the MAF file), or extract from the gene_level mutation file
awk -F"\t" -v OFS=$'\t' '
BEGIN {print "case_id","KRAS","KRAS_VAF","mutation_count","CDKN2A","SMAD4","TP53"}

NR==FNR {
    mut[$1]=$2"\t"$3"\t"$4;next
}

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

$1=="TP53" { # Save TP53 mutations
    tp53[$16]=$9
}

$1=="CDKN2A" { # Save CDKN2A mutations
    cdk[$16]=$9
}
$1=="SMAD4" { # Save SMAD4 mutations
    smad[$16]=$9
}
FNR > 1 && $9!~/Silent/ {
    cnt[$16]+=1
}

END {
for (i in id) {
    kras_id="NA\tNA"; SMAD4_id="NA"; TP53_id="NA"; CDKN2A_id="NA";cnt_id=0
    if (i in kras) {
        kras_id=kras[i]
    }
    if (i in smad) {
        smad_id=smad[i]
    }
    if (i in tp53) {
        tp53_id=tp53[i]
    }
    if (i in cdk) {
        cdk_id=cdk[i]
    }
    if (i in cnt) {
        cnt_id=cnt[i]
    }
    # print i,kras_id,cnt_id,cdk_id,smad_id,tp53_id
    print i,kras_id,cnt_id,mut[i]
}}' <(awk -F"\t" '
$1~/SMAD4|CDKN2A/ || $1=="TP53" || FNR==1 {
    for (i=1;i<=NF;i++) {
        if (i in record==0) {
            record[i]=$i
        } else {
            record[i]=record[i]"\t"$i
        }
    }
}
END {
    for (i in record) {
        r = record[i]
        gsub(/missense_variant/, "Missense", r)
        gsub(/stop_gained/, "Nonsense/frameshift", r)
        gsub(/frameshift_variant/, "Nonsense/frameshift", r)
        gsub(/inframe_deletion/, "In-frame ins/del", r)
        gsub(/inframe_insertion/, "Inframe ins/del", r)

        gsub(/splice_acceptor_variant,/, "", r)
        gsub(/splice_donor_variant,/, "", r)
        gsub(/splice_region_variant,/, "", r)
        gsub(/coding_sequence_variant,/, "", r)
        gsub(/intron_variant,/, "", r)
        gsub(/synonymous_variant,/, "", r)

        gsub(/,splice_acceptor_variant/, "", r)
        gsub(/,splice_donor_variant/, "", r)
        gsub(/,splice_region_variant/, "", r)
        gsub(/,coding_sequence_variant/, "", r)
        gsub(/,intron_variant/, "", r)
        gsub(/,synonymous_variant/, "", r)

        gsub(/splice_acceptor_variant;/, "", r)
        gsub(/splice_donor_variant;/, "", r)
        gsub(/splice_region_variant;/, "", r)
        gsub(/coding_sequence_variant;/, "", r)
        gsub(/intron_variant;/, "", r)
        gsub(/synonymous_variant;/, "", r)

        gsub(/splice_acceptor_variant/, "", r)
        gsub(/splice_donor_variant/, "", r)
        gsub(/splice_region_variant/, "", r)
        gsub(/coding_sequence_variant/, "", r)
        gsub(/intron_variant/, "", r)
        gsub(/synonymous_variant/, "", r)

        split(r, arr, "\t")
        for (j=2;j<=length(arr);j++) {
            split(arr[j], var, ";")
            if (length(var)>1) {
                arr[j]=var[1]

            }
            split(arr[j], var, ",")
            if (length(var)>1) {
                arr[j]=var[1]
            }
            if (arr[j]=="") {
                arr[j]="WT"
            }
        }
        r = arr[1]
        for (j=2;j<=length(arr);j++) {
            r = r"\t"arr[j]
        }

        print r
    }
}' ${DATA}/Mutation_gene_level.cgt)  ${MAF} > ${WORKDIR}/data/1c-caseid_kras_vaf_smad4_tp53_cdk2a_mutation.tsv

awk -F"\t" ' BEGIN {OFS="\t"}
NR==FNR { cnv[$1]=$2"\t"$3;next}
FNR==1 {print $0,"amp","del"}
FNR>1 {print $0,cnv[$1]}
' <(awk -F"\t" -v OFS=$'\t' -v thr=$"0.2" '
FNR == 1{for (i=2;i<=NF;i++) {id[i]=$i;amp[i]=0;del[i]=0}}
FNR > 1 {
    for (i=2;i<=NF;i++) {
        if ($i>thr) {
            # amp[i]++
            ceil_pow = ($i == int($i)) ? 2^$i : int(2^$i) + 1
            amp[i]+=ceil_pow
        } else if ($i<-thr) {
            # del[i]-=1
            ceil_pow = ($i == int($i)) ? 2^$i : int(2^$i) + 1
            del[i]-=ceil_pow
        }
    }
}
END {
    for (i in id) {
        print id[i],amp[i],del[i]
    }
}' ${DATA}/SCNA_log2_gene_level.cct) ${WORKDIR}/data/1c-caseid_kras_vaf_smad4_tp53_cdk2a_mutation.tsv > ${WORKDIR}/data/1c-caseid_kras_vaf_smad4_tp53_cdk2a_mutation_amp_del.tsv
fi

if false;then
awk -F"\t" -v OFS=$'\t' -v thr=$"0.4" '
FNR == 1{
        for (i=2;i<=NF;i++) {id[i]=$i;amp[i]=0;del[i]=0}
    }
FNR > 1 {
        for (i=2;i<=NF;i++) {
            if ($i>thr) {
                # amp[i]++
                amp[i]+=2**$i
            } else if ($i<-thr) {
                # del[i]--
                del[i]-=2**$i
            }
        }
    }
END {
    for (i in id) {
        print id[i],amp[i],del[i]
    }
}' ${DATA}/SCNA_log2_gene_level.cct >  ${WORKDIR}/data/1c_cnv-caseid_amp_del.tsv
fi

# awk -F"\t" '{if ($NF<0) {val=-$NF} else {val=$NF}
#     ; chridx[$1]+=val;cnt[$1]++
# }
# END {for (i in chridx) {print i,chridx[i]}}
# ' ${DATA}/SCNA_log2_segment_level.cct



# Mutation comparison between $52 vs $9
# awk -F"\t" ' {split($52,mut,",");for (i in mut){print mut[i],$9}}' data_linkedomics/Mutation_results.maf| sort | uniq -c| sort -k3