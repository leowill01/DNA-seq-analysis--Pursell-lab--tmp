#!/anaconda3/bin/python

# Import libraries
from matplotlib import pyplot as plt
import numpy as np
import sys
from matplotlib_venn import venn3, venn3_circles

# Take arguments
a_only = int(sys.argv[1])
b_only = int(sys.argv[2])
a_and_b = int(sys.argv[3])
c_only = int(sys.argv[4])
a_and_c = int(sys.argv[5])
b_and_c = int(sys.argv[6])
a_b_and_c = int(sys.argv[7])
plt.title(str(sys.argv[8]))
set1 = sys.argv[9]
set2 = sys.argv[10]
set3 = sys.argv[11]

# Print out arguments entered
print(sys.argv)

# Set up Venn
v = venn3(
subsets=(a_only, b_only, a_and_b, c_only, a_and_c, b_and_c, a_b_and_c),
set_labels=(set1, set2, set3))

# Adjust 3rd set label position
lbl = v.get_label_by_id("C")
x, y = lbl.get_position()
lbl.set_position((x+0.2, y))

# v.get_patch_by_id('100').set_alpha(1.0) # what does this do
# v.get_patch_by_id('100').set_color('#f06449')

# Customize circle lines
c = venn3_circles(
subsets=(a_only, b_only, a_and_b, c_only, a_and_c, b_and_c, a_b_and_c),
linewidth=2
)

c[0].set_lw(2)
c[0].set_ls('dotted')

c[1].set_lw(2)
c[1].set_ls('dashed')

c[2].set_lw(2)

# Show plot
plt.show()
