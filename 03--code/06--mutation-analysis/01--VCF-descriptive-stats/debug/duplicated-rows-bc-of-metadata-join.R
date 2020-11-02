# DEBUG: duplicated sample filepaths
test_absfp_tesv = read_tsv("/Volumes/Pursell-Lab-HD/Pursell-lab/02--projects/project_02--DNA-seq-analysis/04--analysis/06--mutation-analysis/BRCA1-related--Jackson/2020-10-01--mutect2-only-1st-run/03--results-desc-stats/results--2020-10-02-152206/tmp/absoluteFilepathsVCFs.tsv")

# the info table includes duplicates - likely from when I joined the metadata table that included the same samples from multiple sequencing sets. Ill have to change the make info table script to remove duplicates in the metadata.
