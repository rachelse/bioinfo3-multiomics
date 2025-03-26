import sys
import matplotlib.pyplot as plt

x = []
y = []
# types = ["G12C","G12D","G12R","G12S","G12V","G13D","NA","Q61H","Q61R"]
# colors = ['r','g','b','c','m','y','k','BLUE','orange']
# type2color = {}
# for i in range(len(types)):
#     type2color[types[i]] = colors[i]
# color2type = {}
# for i in range(len(types)):
#     color2type[colors[i]] = types[i]
# this = []
# colormap = []
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip("\n").split()
        x.append(line[0])
        # this.append(line[1])
        # colormap.append(type2color[line[1]])
        y.append(float(line[1]))

# Sort by y
# x, y, this, colormap = zip(*sorted(zip(x, y, this, colormap), key=lambda x: x[1], reverse=True))

# # Read the histology.real
#     name_list = []
#     with open("histology.real") as f:
#         line = f.readline()
#         while line:
#             line = line.strip().split("\t")
#             name_list.append(line[0])
#             line = f.readline()

# Draw the figure
fig, ax = plt.subplots()
# Reverse the color
ax.imshow([y], cmap='Blues', interpolation='nearest')
ax.set_yticks([0])
ax.set_yticklabels([''])
# ax.set_xticks(range(len(sample2idx_sorted)))
# ax.set_xticklabels(list(sample2idx_sorted.keys()), rotation=90, fontsize=2)
# No  x ticks
ax.set_xticks([])
# colorbar
# cbar = plt.colorbar(ax.imshow([list(sample2idx_sorted.values())], cmap='Greys', interpolation='nearest'), fraction=0.046, pad=0.04)
plt.savefig("deconvolution.png", bbox_inches='tight', dpi=500)