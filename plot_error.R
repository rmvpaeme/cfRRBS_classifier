library(Rmisc)
library(tidyverse)
require(reshape2)
require(data.table)
library(optparse)

#opt$directory <- "/Users/rmvpaeme/Repos/cfRRBS_classifier_v0.5/classifySamples/output/classification/"
#opt$outdir <- "./classifySamples/output/plots/"
#opt$normal <- "wbc,normal"

option_list = list(
  make_option(c("-d", "--directory"), type="character", default="./classifySamples/output/classification/", 
              help="directory for *deconv_output.csv files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./classifySamples/output/plots/", 
              help="output directory for figures [default= %default]", metavar="character"),
  make_option(c("-n", "--normal"), type="character", default=NULL, 
              help="comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment (after the deconvolution step) e.g. cfdna,wbc", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

normal <- unlist(strsplit(opt$normal, ","))

myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.6,
                 dodge.width = 0.2,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

data_path <- opt$directory
files<-list.files(data_path,recursive=TRUE)
files<-files[grep("deconv_output.csv", files)]

data_list <- list()
for(i in 1:length(files)){
  print(files[i])
  tmp<-data.table::fread(paste0(data_path,files[i]), header=T, sep=",",  data.table=FALSE)
  data_list[[i]] <- tmp
}

# Find what the predicted tumor is within each iteration and plot that
classFreq <- data_list
classFreq <- lapply(classFreq, melt)
classFreq <- classFreq %>% lapply(. %>% filter((V1 != "normal")) %>% filter((V1 != "wbc") ))
grouped <- classFreq %>% lapply(. %>% dplyr::group_by(variable) %>%
                                  filter(value == max(value))  %>%
                                  arrange(variable, V1)) #%>% dplyr::summarise(value = max(value)))
grouped <- do.call(rbind, grouped)

grouped$variable <- gsub("snake", "", grouped$variable )

p1 <- ggplot(grouped, aes(x = V1)) + geom_bar() + facet_wrap(~variable, scales="free") + coord_flip()
ggsave(paste0(opt$outdir, "/", format(Sys.time(),'%Y%m%d_%H%M%S_'), "countplot.png"), p1, width = 16, height = 16, limitsize = FALSE)

grouped <- as.data.table(grouped)
grouped <- grouped[grouped[, .I[value == max(value)], by=variable]$V1]
grouped <- distinct(grouped)

classificationResults_plasma <- as.data.frame(grouped)
colnames(classificationResults_plasma) <- c("ClassificationResult", "SampleID", "eTFx")

# Calculate tumor percentage with SD. 
full_results = do.call(rbind, data_list)
full_results <- full_results %>% melt()

tgc <- summarySE(full_results, measurevar="value", groupvars=c("V1", "variable")) %>% filter((!V1 %in% c(normal)))

grouped <- as.data.table(tgc)

grouped <- grouped[grouped[, .I[value == max(value)], by=variable]$V1]
classificationResults_plasma <- as.data.frame(grouped)
colnames(classificationResults_plasma) <- c("ClassificationResult", "SampleID", "len", "tumorFx", "sd", "se", "ci")

pd <- position_jitterdodge(jitter.width = NULL, jitter.height = 0,
                           dodge.width = 0.75, seed = 1)

classificationResults_plasma$SampleID <- gsub("snake", "", classificationResults_plasma$SampleID)

p2 <- ggplot(classificationResults_plasma, aes(x = SampleID, y = tumorFx*100, col = ClassificationResult)) +
  geom_point(position = myjit) + theme_bw() +    geom_errorbar(aes(ymin=tumorFx*100-sd*100, ymax=tumorFx*100+sd*100), width=.1, position=myjit) +
  theme(panel.grid.major.y=element_line(linetype = "dashed",color="gray88"), axis.text.y=element_text(size=8)) + 
  labs(y = "Estimated TFx (%)", x = "SampleID", col = "Classification call") + coord_flip()
ggsave(paste0(opt$outdir, "/", format(Sys.time(),'%Y%m%d_%H%M%S_'), "errorplot.png"), p2)
write_tsv(classificationResults_plasma, paste0(opt$outdir, "/", format(Sys.time(),'%Y%m%d_%H%M%S_'), "aggregated_results.tsv"))
