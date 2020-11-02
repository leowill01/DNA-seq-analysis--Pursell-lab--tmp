cd "/lustre/project/zpursell/leo/Pursell-lab-HD--MIRROR/Pursell-lab/02--projects/project_02--DNA-seq-analysis/04--analysis/03--variant-calling--postprocessing--annotation/BRCA1-related/all-samples-as-of-2020-09-20/mutect/03--tumor-normal-somatic-variants--Mutect2/samples-with-matched-normals"

in_table="tumor-normal-table.tsv"

headerless_table="headerless-table.tsv"

sed 1d "$in_table" > "$headerless_table"


while read line ; do
	arg=$(echo "$line" | cut -f 2)
	echo "command -option "$arg""
	arg=$(echo "$line" | cut -f 2)
	echo "command -option "$arg""
done < "$headerless_table"

# sed pipe doesnt work?
sed 1d "$in_table" | \
while read line; do
	arg=$(echo "$line" | cut -f 2)
	echo "command -option "$arg""
done