#!/bin/bash

Rscript fig7a.r
python fig7_clsuter.py H.txt $1
awk '{print $2}' $1 | sort | uniq -c
awk 'NR==FNR {a[$1];next;} ($1 in a) {print}' gt_member1.list <(awk '($2 == 1) {print $1}' $1) | wc -l
awk 'NR==FNR {a[$1];next;} ($1 in a) {print}' gt_member1.list <(awk '($2 == 2) {print $1}' $1) | wc -l
awk 'NR==FNR {a[$1];next;} ($1 in a) {print}' gt_member2.list <(awk '($2 == 1) {print $1}' $1) | wc -l
awk 'NR==FNR {a[$1];next;} ($1 in a) {print}' gt_member2.list <(awk '($2 == 2) {print $1}' $1) | wc -l
