ncol(nmf_res_contrib_norm_mouse)
nrow(nmf_res_contrib_norm_mouse)
length(params$info_table_cluster_categories_colors)

# why this error:
	# Error in check.length("fill") :
	# 	'gpar' element 'fill' must not be length 0
# when plotting:
	plotSigContribCustom(nmf_res_contrib_norm_human, speciesName = "human")
# when this works:
	plotSigContribCustom(nmf_res_contrib_norm_mouse, speciesName = "mouse")

View(nmf_res_contrib_norm_human)
View(nmf_res_contrib_norm_mouse)

# problem is when a function inside pheatmap() calls the gpar() function and uses the local annotation_colors var

# view pheatmap() to find where annotation_colors is assigned ------------------------
pheatmap

# # ...
# 		  annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
# # ...
# 	if (length(annotation) != 0) {
# 		annotation_colors = generate_annotation_colours(annotation,
# 														annotation_colors, drop = drop_levels)
# 	}
# 	else {
# 		annotation_colors = NA
# 	}
# # ...
# 					   annotation_colors = annotation_colors, annotation_legend = annotation_legend,

# define local functions to pheatmap() that are used to define annotation_colors
# def is.na2 function
is.na2 = function(x){
	if(is.list(x) | length(x) > 1){
		return(FALSE)
	}
	if(length(x) == 0){
		return(TRUE)
	}

	return(is.na(x))
}

# def generate anno colors func
generate_annotation_colours_TEST = function(annotation, annotation_colors, drop = F){
	if(is.na2(annotation_colors)){
		annotation_colors = list()
	}
	count = 0
	for(i in 1:length(annotation)){
		annotation[[i]] = annotation[[i]][!is.na(annotation[[i]])]
		if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
			if (is.factor(annotation[[i]]) & !drop){
				count = count + length(levels(annotation[[i]]))
			}
			else{
				count = count + length(unique(annotation[[i]]))
			}
		}
	}

	factor_colors = dscale(factor(1:count), hue_pal(l = 75))

	oldseed = NULL
	if (exists(".Random.seed"))
		oldseed = get(".Random.seed", pos=.GlobalEnv)

	set.seed(3453)

	cont_counter = 2
	for(i in 1:length(annotation)){
		if(!(names(annotation)[i] %in% names(annotation_colors))){
			if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
				n = length(unique(annotation[[i]]))
				if (is.factor(annotation[[i]]) & !drop){
					n = length(levels(annotation[[i]]))
				}
				ind = sample(1:length(factor_colors), n)
				annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
				l = levels(as.factor(annotation[[i]]))
				l = l[l %in% unique(annotation[[i]])]
				if (is.factor(annotation[[i]]) & !drop){
					l = levels(annotation[[i]])
				}

				names(annotation_colors[[names(annotation)[i]]]) = l
				factor_colors = factor_colors[-ind]
			}
			else{
				annotation_colors[[names(annotation)[i]]] = brewer_pal("seq", cont_counter)(5)[1:4]
				cont_counter = cont_counter + 1
			}
		}
	}

	if(!is.null(oldseed)){
		assign(".Random.seed", oldseed, pos=.GlobalEnv)
	}
	else{
		remove(.Random.seed, pos=.GlobalEnv)
	}

	return(annotation_colors)
}

# try to generate colors
generate_annotation_colours_TEST(annotations, ann_colors)
# ERROR:
# Error in `[[<-.data.frame`(`*tmp*`, i, value = c("ovary", "lung", "lung",  :
# 												 	replacement has 311 rows, data has 312


# RESTARTED ----------------------------

# still getting same error
plotSigContribCustom(nmf_res_contrib_norm_mouse,
					 speciesName = "Mouse") # works

plotSigContribCustom(nmf_res_contrib_norm_human,
					 speciesName = "Human") # fails
# Error in check.length("fill") :
# 	'gpar' element 'fill' must not be length 0
# Called from: check.length("fill")

# view the norm contrib matrices --------------------------------

View(nmf_res_contrib_norm_mouse)
View(nmf_res_contrib_norm_human) # seems like same structure

# check structure --------------------------------
str(nmf_res_contrib_norm_mouse)
# num [1:37, 1:2] 0.811 1 0.915 0.75 0.232 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:37] "7887-MG--TN-set2-VAF005" "7895-Int--TN-set1-VAF005" "7895-MG--TN-set1-VAF005" "7895-MG--TN-set2-VAF005" ...
# ..$ : chr [1:2] "1" "2"
head(nmf_res_contrib_norm_mouse)
#								   1            2
# 7887-MG--TN-set2-VAF005  0.8111356 1.888644e-01
# 7895-Int--TN-set1-VAF005 0.9999994 6.382972e-07
# 7895-MG--TN-set1-VAF005  0.9149815 8.501853e-02
# 7895-MG--TN-set2-VAF005  0.7499607 2.500393e-01
# 7895-tumor               0.2324242 7.675758e-01
# 8671-tumor               0.5555739 4.444261e-01


str(nmf_res_contrib_norm_human)
# num [1:275, 1:2] 0.052 0.0672 0.0629 0.0728 0.7849 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:275] "TCGA-04-1357" "TCGA-05-4382" "TCGA-05-4424" "TCGA-06-0649" ...
# ..$ : chr [1:2] "1" "2"
head(nmf_res_contrib_norm_human)
# 					    1         2
# TCGA-04-1357 0.05195940 0.9480406
# TCGA-05-4382 0.06716020 0.9328398
# TCGA-05-4424 0.06293137 0.9370686
# TCGA-06-0649 0.07278898 0.9272110
# TCGA-06-5416 0.78487376 0.2151262
# TCGA-06-5858 0.07574448 0.9242555

# seems like same structure...

# plot with regular pheatmap() function ----------------------------
pheatmap(deNovoNMFResMouse$contribution, cluster_cols = F)
pheatmap(deNovoNMFResHuman$contribution, cluster_cols = F)
# both work

# try with normalized contrib matrices
pheatmap(nmf_res_contrib_norm_mouse, cluster_cols = F)
pheatmap(nmf_res_contrib_norm_human, cluster_cols = F)
# both work

# stepwise `pheatmap()` check ----------------------------
# stepwise complicate the plot to see where failure occurs:
pheatmap(mat = nmf_res_contrib_norm_human,
		 cluster_cols = F,
		 annotation_row = annotations) # FAILS: check `annotations`

# check `annotations`
annotations # lots of NA info for human samples

# stepwise check mouse
pheatmap(mat = nmf_res_contrib_norm_mouse,
		 cluster_cols = F,
		 annotation_row = annotations)
# WORKS - THERE IS SOMETHING WRONG WITH `annotations` wrt. THE HUMAN SAMPLES

# fix `annotations` ----------------------------
# see which annotations are being selected when only human samples are plotted:
glimpse(annotations)
# Rows: 312
# Columns: 8
# $ species      <chr> "human", "human", "human", "human", "human", "human", "human", "…
# $ is_tumor     <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",…
# $ tissue       <chr> "ovary", "lung", "lung", "brain", "brain", "brain", "ovary", "br…
# $ BRCA1_status <chr> "mut", "mut", "mut", "mut", "mut", "mut", "mut", "mut", "mut", "…
# $ p53_status   <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
# $ stress_type  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
# $ IR_status    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
# $ MMC.dose     <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …

str(annotations)
# 'data.frame':	312 obs. of  8 variables:
# $ species     : chr  "human" "human" "human" "human" ...
# $ is_tumor    : chr  "1" "1" "1" "1" ...
# $ tissue      : chr  "ovary" "lung" "lung" "brain" ...
# $ BRCA1_status: chr  "mut" "mut" "mut" "mut" ...
# $ p53_status  : chr  NA NA NA NA ...
# $ stress_type : chr  NA NA NA NA ...
# $ IR_status   : chr  NA NA NA NA ...
# $ MMC.dose    : chr  NA NA NA NA ...

# filter for only human samples
annotations %>%
	filter(species == "human")
# `p53_status`, `stress_type`, `IR_status`, and `MMC.dose` are completely `NA`

# make ad hoc annotations subset just to use with humans excluding the completely `NA` variables
annotations_human = annotations %>%
	filter(species == "human") %>%
	select(species, is_tumor, tissue, BRCA1_status)

annotations_mouse = annotations %>%
	filter(species == "human") %>%
	select(any_of(params$info_table_cluster_categories_mouse))

# try plot with just ad hoc human annotations
pheatmap(mat = nmf_res_contrib_norm_human,
		 cluster_cols = F,
		 annotation_row = annotations_human)
# FUCKING WORKS omfg

# does multispp work?
pheatmap(mat = nmf_res_contrib_norm_multispp,
		 cluster_cols = F,
		 annotation_row = annotations) # WORKS
# it must be unable to plot when an ENTIRE column is NA, but if there is at least one row with a value in that column, it will plot the annotation

# SOLUTION: make species specific subsets from `annotations` ---------------------

# also try to fix annotation colors ------------------------------------------
# plot with master annotation colors
pheatmap(mat = nmf_res_contrib_norm_human,
		 cluster_cols = F,
		 annotation_row = annotations_human,
		 annotation_colors = params$info_table_cluster_categories_colors)
# works - no need to change

# try to dynamically select annotaions based on species-specific cluster categories in params

annotations_human = annotations %>%
	filter(species == "human") %>%
	select(any_of(params$info_table_cluster_categories_human))
# plot again
pheatmap(mat = nmf_res_contrib_norm_human,
		 cluster_cols = F,
		 annotation_row = annotations_human,
		 annotation_colors = params$info_table_cluster_categories_colors)
