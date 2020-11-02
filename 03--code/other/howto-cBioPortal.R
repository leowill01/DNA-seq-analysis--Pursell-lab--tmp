# cBioPortal How-To
library('cgdsr')

# Create CGDS object
mycgds = CGDS("https://www.cbioportal.org/")

test(mycgds)

# Get list of cancer studies at server
stud = getCancerStudies(mycgds)
str(stud)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[2,1]
my_study_case_lists = getCaseLists(mycgds,mycancerstudy)
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
genetic_profiles = getGeneticProfiles(mycgds,mycancerstudy)
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[3,1]

# Get data slices for a specified list of genes, genetic profile and case list
gen_profile_data = getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

# documentation
help('cgdsr')
help('CGDS')
