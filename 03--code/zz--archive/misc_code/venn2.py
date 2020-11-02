#!/anaconda3/bin/python

# USAGE: python venn2.py [a_only] [b_only] [a_and_b] [plot title] [set1set1Name] [set2set2Name]

# Import libraries
from matplotlib import pyplot as plt
import numpy as np
import sys
from matplotlib_venn import venn2, venn2_circles

# Enter arguments
a_only = int(sys.argv[1])
b_only = int(sys.argv[2])
a_and_b = int(sys.argv[3])
plt.title(str(sys.argv[4]))
set1Name = str(sys.argv[5])
set2Name = str(sys.argv[6])

# Print out arguments entered
print(sys.argv)

# Set up Venn diagram
v = venn2(subsets=(a_only, b_only, a_and_b),
set_labels=(set1Name, set2Name))

# Adjust labels
lblA = v.get_label_by_id("A")
x, y = lblA.get_position()
lblA.set_position((x-0.5, y+0.3))

lblB = v.get_label_by_id("B")
x, y = lblB.get_position()
lblB.set_position((x+0.45, y+0.3))

# Customize circle lines
c = venn2_circles(
subsets=(a_only, b_only, a_and_b),
linewidth=1
)


# Show plot
plt.show()
