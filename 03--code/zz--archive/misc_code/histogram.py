import pandas as pd
import matplotlib.pyplot as plt

# use a newline-separated text file with AF values
in_text_file = '/Users/leo/Desktop/temp/m33_N_MEF_m17_T_het.AF.txt'

allele_freq = pd.read_csv(in_text_file, delimiter='\n')

allele_freq.hist(bins=10)
plt.xlabel('Allele Frequency')
plt.ylabel('Number of Calls')
plt.title('AF: m33_N_MEF v. m17_T_het')
plt.grid(False)
plt.xlim(left=0, right=100)
plt.savefig('/Users/leo/Desktop/image.png')
