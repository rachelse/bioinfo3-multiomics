import sys
import numpy as np
import matplotlib.pyplot as plt

# [Sample, start, end, log2_copyratio]
def process_seg_level(file):
    samples = []
    start_end = []
    fold_change = []
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
        prev_sample = ""
        while line:
            line = line.strip("\n").split("\t")
            sample = line[0]
            if prev_sample != sample:
                samples.append(sample)
                start_end.append([])
                fold_change.append([])
                prev_sample = sample
            start_end[-1].append([int(line[2]), int(line[3])])
            fold_change[-1].append(float(line[-1]))
            line = f.readline()

    for i in range(len(fold_change)):
        fold_change[i] = np.array(fold_change[i])
    
    return samples, start_end, fold_change

def process_gene_level(file):
    lines = []
    with open(file) as f:
        line = f.readline()
        line = line.strip("\n").split("\t")
        samples = line[1:]
        # print(f"Sample count: {len(samples)}")
        for i in range(len(samples)):
            lines.append([])
        line = f.readline()
        while line:
            line = line.strip("\n").split("\t")
            assert len(line) == 141, f"Error: {len(line)}, {line[0]}"
            for i in range(1, len(line)):
                if line[i] != '':
                    lines[i-1].append(float(line[i]))
                else:
                    lines[i-1].append(0)
            line = f.readline()

    for i in range(len(lines)):
        lines[i] = np.array(lines[i])
    
    return samples, lines

if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage: python cnv_fig1.py <segment/gene> <file>")
        sys.exit(1)
    if sys.argv[1] == "segment":
        samples, start_end, fold_change = process_seg_level(sys.argv[2])
    else:
        samples, fold_change = process_gene_level(sys.argv[2])
    
    sample2val = {}
    neg_sample2val = {}
    for i in range(len(samples)):
        pos_idx = np.where(fold_change[i] > 0.2)[0]
        neg_idx = np.where(fold_change[i] < -0.2)[0]
        # Round up
        fold_change[i] = np.power(2, fold_change[i])
        fold_change_all = np.ceil(fold_change[i]).astype(int)
        # fold_change[i][fold_change[i] > 4] = 4
        # fold_change[i][fold_change[i] <= 1] += 1
        sample2val[samples[i]] = np.sum(fold_change_all[pos_idx])
        # neg_sample2val[samples[i]] = np.sum(fold_change[i][neg_idx])
        # sample2val[samples[i]] = len(pos_idx)
        neg_sample2val[samples[i]] = len(neg_idx)

    # sort the sample2val by value
    sample2val = dict(sorted(sample2val.items(), key=lambda item: item[1], reverse=True))
    # for sample in sample2val:
    #     print(sample, sample2val[sample], neg_sample2val[sample])

    # Read the histology.real
    name_list = []
    with open("histology.real") as f:
        line = f.readline()
        while line:
            line = line.strip().split("\t")
            name_list.append(line[0])
            line = f.readline()

    sample2val_sorted = {}
    neg_sample2val_sorted = {}
    for name in name_list:
        if name in sample2val:
            sample2val_sorted[name] = sample2val[name]
        if name in neg_sample2val:
            neg_sample2val_sorted[name] = -neg_sample2val[name]
        else:
            print(f"Warning: {name} not found in the CNV file")
            exit(-1)

    # Plot the bar plot for the sample2val_sorted and neg_sample2val_sorted in the same plot
    fig, ax = plt.subplots()
    width = 0.35
    x = np.arange(len(sample2val_sorted))
    ax.bar(x, sample2val_sorted.values(), width, label="CNV amplification", color='r')
    ax.bar(x, neg_sample2val_sorted.values(), width, label="CNV deletion", color='b')
    # ax.set_xticks(x)
    # ax.set_xticklabels(sample2val_sorted.keys(), rotation=90, fontsize=2)
    # No xticks
    ax.set_xticks([])
    ax.legend()
    plt.savefig("cnv_fig1.png", dpi=300)