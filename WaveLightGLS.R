# WaveLightGLS
#
# Method for analysing GLS data using wavelet
# A. Roy
#----------------------------------------------------

### LOADING DATA ----------------

# Read Metadata and select relevant data
metadata = "./data/Metadata_GLS.csv"
metadata <- read.csv(metadata, header=TRUE, sep=",")

# load graphic data
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
world <- fortify(wrld_simpl)

eez<-readOGR("./data/World_EEZ.shp", "World_EEZ")
EEZ <-  fortify(eez)
EEZ_br <- EEZ[which(EEZ$id==163 & !EEZ$hole),]

### PERFORM ANALYSIS ON EVERY BIRD ----------------

# get list of gls files
list.gls = list.files("./data/gls/", pattern=".lig", all.files=FALSE,full.names=FALSE)
list.gls = gsub('_000.lig', '', list.gls)

# gls = list.gls[[1]]
for (gls in list.gls){
  
  lig <- paste0("./data/gls/", gls, "_000.lig")
  act <- paste0("./data/gls/", gls, "_000.act")
  tem <- paste0("./data/gls/", gls, "_000.tem")
  
  try(
    rmarkdown::render('./WaveLightGLS.Rmd',
                      output_format = "html_document",
                      output_file =  paste0(gls, ".html"),
                      output_dir = "./results")
  )
}

