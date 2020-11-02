plotCOSMICContributionCustom = function(fit_res, selectCtbSigs, speciesName = NULL) {
	# plot relative contribution barplot
	pc3 = plot_contribution(contribution = fit_res$contribution[selectCtbSigs,, drop = F],
							signatures = cancer_signatures[,selectCtbSigs],
							mode = "relative")


fitResHuman[["contribution"]] %>% head()
fitResHuman[["reconstructed"]] %>% head()

plotCOSMICContributionCustom(fit_res = fitResHuman, selectCtbSigs = selectCtbSigsHuman, speciesName = "human")

# human
plot_contribution(contribution = fitResHuman$contribution[selectCtbSigsHuman,, drop = F],
				  signatures = cancer_signatures[,selectCtbSigsHuman],
				  mode = "absolute") +
	labs(title = "Mutation Contribution of COSMIC Signatures (Relative)") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
# still throwing error

# ok, compare human args with mouse
View(fitResMouse)
glimpse(fitResMouse)
View(fitResHuman)
glimpse(fitResHuman)
# no apparent functional differences

View(fitResMouse$contribution)
glimpse(fitResMouse$contribution)
View(fitResHuman$contribution)
glimpse(fitResHuman$contribution)
# no apparent functional differences

View(fitResMouse$reconstructed)
glimpse(fitResMouse$reconstructed)
View(fitResHuman$reconstructed)
glimpse(fitResHuman$reconstructed)
# no apparent functional differences

# try mouse
plot_contribution(contribution = fitResMouse$contribution[selectCtbSigsMouse,, drop = F],
				  signatures = cancer_signatures[,selectCtbSigsMouse],
				  mode = "absolute") +
	labs(title = "Mutation Contribution of COSMIC Signatures (Relative)") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
# still throwing error

# FIXED: error was being thrown from the `width = ` arg for ggsave
# width = 4 + 0.2 * ncol(fit_res$contribution[selectCtbSigs,]),
4 + 0.2*ncol(fitResMouse$contribution[selectCtbSigsMouse,])
