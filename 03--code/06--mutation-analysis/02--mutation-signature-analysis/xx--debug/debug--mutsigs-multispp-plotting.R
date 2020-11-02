# ----
view(fitResHuman$contribution[selectCtbSigsHuman,])
view(cancer_signatures[, selectCtbSigsHuman])

view(fitResMouse$contribution[selectCtbSigsMouse,])
view(cancer_signatures[, selectCtbSigsMouse])

# the human fit res dont have a sample name.. :
fitResHuman = fit_to_signatures(mutMatHuman, cancer_signatures)
mutMatHuman
# ... but they do? so where does it go missing when subsetting 'selectCtbSigsHuman'?

fitResHuman$contribution %>% head()
fitResHuman$contribution[selectCtbSigsHuman,] %>% head()

# do the mouse samples lose sample names?
fitResMouse$contribution %>% head()
fitResMouse$contribution[selectCtbSigsMouse,] %>% head() # NO!

# so somewhere when subsetting for contributing signatures and when there's only 1 sample
# let's look at the structures:
# human
	class(fitResHuman$contribution %>% head())
	class(fitResHuman$contribution[selectCtbSigsHuman,] %>% head()) # this is a numeric vector when it should be an array
#mouse
	class(fitResMouse$contribution %>% head())
	class(fitResMouse$contribution[selectCtbSigsMouse,] %>% head()) # correctly a numeric matrix


as.matrix(fitResHuman$contribution) %>% head() # same result
as.matrix(fitResHuman$contribution[selectCtbSigsHuman,, drop = F], dimnames = list(rownames()))

names(fitResHuman$contribution) # NULL
colnames(fitResHuman$contribution) #

# ----
# trying to solve error:
# 	Error in grid.newpage() : non-finite location and/or size for viewport

# error thrown for human data
plotCOSMICContributionCustom(fit_res = fitResHuman,
							 selectCtbSigs = selectCtbSigsHuman,
							 speciesName = "human")
	# Error during wrapup: non-finite location and/or size for viewport
	# Error: no more error handlers available (recursive errors?); invoking 'abort' restart

# does it happen for mouse? YES
plotCOSMICContributionCustom(fit_res = fitResMouse,
							 selectCtbSigs = selectCtbSigsMouse,
							 speciesName = "mouse")
