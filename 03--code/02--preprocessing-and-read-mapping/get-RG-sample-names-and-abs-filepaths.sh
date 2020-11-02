idev -c 1
module load samtools

cd "/lustre/project/zpursell/leo/Pursell-lab-HD--MIRROR/Pursell-lab/02--projects/project_02--DNA-seq-analysis/04--analysis/02--preprocessing-and-read-mapping/RRM1-related/05--ribont-mut-mice-targetseq--Chabes/02--map-reads-and-clean"

printf "bam_rg_sample_name\tabs_filepath\n" > bam-RG-sample-names-and-abs-filepaths.tsv
for i in $(find . -name "*.bam") ; do
	# get BAM read group sample name
	echo "Processing: " "$(basename "$i")"
	bam_RG_sample_name="$(
		samtools view -H "$i" | \
		grep '^@RG' | \
		cut -f 2 | \
		sed 's/ID://'
	)"
	echo "$bam_RG_sample_name"

	# get BAM absolute filepath
	bam_abs_filepath=$(
		echo "$(cd "$(dirname "$i")"; pwd -P)/$(basename "$i")"
	)
	echo "$bam_abs_filepath"

	# print both to a TSV line
	echo -e "${bam_RG_sample_name}\t${bam_abs_filepath}" >> bam-RG-sample-names-and-abs-filepaths.tsv
done