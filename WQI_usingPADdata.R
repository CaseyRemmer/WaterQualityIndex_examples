


#  Load packages & data
# --------------------------------------------------------- #
# Load packages
library('factoextra')
library('FactoMineR')
library('gridExtra')
library('ggplot2')

lat.longs<-read.csv("C:/Users/crrem/Dropbox/documents/UWaterloo/PAD/water chem/5 years pca/GIS/PAD_site_trimmed_for water chem.csv")

data <- read.csv("C:/Users/Casey/Dropbox/documents/UWaterloo/PAD/water chem/5 years pca/waterchem_15to19_logtransformed_M2M5and18 deleted.csv")
## this data has already been log-transformed
head(data)# PCA ordination
lake.pca <- PCA(data[2:18],
                ind.sup = 842:975,		# Rows of passive data NOT TO BE USED in PCA
                scale.unit = TRUE,		# scale incoming data
                ncp = 5,			# number of PCs to retain
                graph = FALSE)

spring15_data<-data %>% filter (grepl('M15', Site)) %>% filter (!grepl('S|J|j', Site))

spring2015_pca <- PCA(spring15_data[2:18],
                      scale.unit = TRUE,		# scale incoming data
                      ncp = 5,			# number of PCs to retain
                      graph = FALSE)


spring2015_pca$var$coord

# print(lake.pca)		# This gives a list of all the information in the PCA

# Apply 'varimax' rotation to PCA variables and samples
loadings <- lake.pca$var$coord %*% diag(sqrt(lake.pca$eig[,1:1]), 5, 5)	# make sure to change number of axes
var.vmax <- varimax(loadings)$loadings
ind.vmax <- scale(lake.pca$ind$coord) %*% varimax(loadings)$rotmat
lake.pca$var$coord <- var.vmax
lake.pca$ind$coord <- ind.vmax


##for spring PCA
loadings <- spring2015_pca$var$coord %*% diag(sqrt(lake.pca$eig[,1:1]), 5, 5)	# make sure to change number of axes
var.vmax <- varimax(loadings)$loadings
ind.vmax <- scale(spring2015_pca$ind$coord) %*% varimax(loadings)$rotmat
spring2015_pca$var$coord <- var.vmax
spring2015_pca$ind$coord <- ind.vmax

# --------------------------------------------------------- #
# Plot Outputs
# --------------------------------------------------------- #
# Set plot aesthetics
theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))		

# Eigenvalues &  Scree plot
round(lake.pca$eig, digits=2)	# Eigenvalues

scree.lab <- sprintf("%0.2f", round(lake.pca$eig[,1:1], digits = 2))
scree.p <- fviz_eig(lake.pca, choice = "eigenvalue", geom = "bar", addlabels = F, label = scree.lab) +
  geom_text(label = scree.lab[1:10], nudge_y=1) + theme +
  geom_hline(yintercept=1, linetype=2, color="grey")
windows(5,5); scree.p

# Colours: extra data in blue, baseline data in grey, vectors in red
# Check out link below for more info on making plots
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# Make plots: PC dimenions 1 & 2
# PCA baseline plot
p1 <- fviz_pca_biplot(lake.pca, title="PCA Dimensions 1 & 2",
                      invisible="ind.sup",	# Hide river data (added below)
                      repel=TRUE,	# Avoid label overlap
                      
                      # Vector Parameters
                      label="",			# No vector labels (too many vars)
                      col.var="red",		# Vector/label colour
                      
                      # Ordination Data Parameters
                      col.ind="darkgrey",	# Point colour
                      pointshape=1,		# Point shape (pch values)
                      
) + theme				# Plot aesthetic parameters (above)

# Add passive data to plots
p1_r <- fviz_add(p1, lake.pca$ind.sup$coord,
               geom="point",		# Geometry ("point", "text", "arrows", combinations)
               color="#0072B2",		# Colour
               shape=15,			# Point shape (pch values)
               addlabel=FALSE,		# Label points (T/F)
               labelsize=3,		# Label text size
               repel=TRUE,			# Avoid label overlap
)

# Make plots: PC dimenions 3 & 4
# PCA baseline plot
p2 <- fviz_pca_biplot(lake.pca, title="PCA Dimensions 3 & 4",
                      axes=c(3,4),		# PC dim 3 & 4
                      invisible="ind.sup",	# Hide monitoring data (added below)
                      repel=TRUE,			# Avoid label overlap
                      
                      # Vector Parameters
                      label="",			# No vector labels (too many vars)
                      col.var="red",		# Vector/label colour
                      
                      # Ordination Data Parameters
                      col.ind="darkgrey",	# Point colour
                      pointshape=1,		# Point shape (pch values)
                      
) + theme				# Plot aesthetic parameters (above)

# Add passive data to plots
p2 <- fviz_add(p2, lake.pca$ind.sup$coord,
               axes=c(3,4),		# PC dim 3 & 4
               geom="point",		# Geometry ("point", "text", "arrows", combinations)
               color="#0072B2",		# Colour
               shape=15,			# Point shape (pch values)
               addlabel=FALSE,		# Label points (T/F)
               labelsize=3,		# Label text size
               repel=TRUE,			# Avoid label overlap
)

# Final PCA plot
# Graph results in 1 figure (2 panels)
windows(9,5); grid.arrange(p1, p2, ncol=2)

# --------------------------------------------------------- 

fviz_contrib(lake.pca, choice = "var", axes = 1)
fviz_contrib(lake.pca, choice = "var", axes = 2)
fviz_contrib(lake.pca, choice = "var", axes = 3)
fviz_contrib(lake.pca, choice = "var", axes = 4)

##calculating the WQI

#WQI = (w1 * PC1) + (w2 * PC2) + ... + (wn * PCn)

#WQI: The Water Quality Index (overall score).
#wi: The weight assigned to each principal component.
#PCi: The score for the ith principal component.


# Calculate the scores for the selected principal components
pca_scores <- lake.pca$ind$coord[, c("Dim.1", "Dim.2")]

w1 <- 0.323 ## %variance explained
w2 <- 0.17

# Calculate the Water Quality Index (WQI)
wqi <- w1 * pca_scores[1:15, "Dim.1"] + w2 * pca_scores[1:15, "Dim.2"]
wqi<-as.data.frame(wqi)

wqi$Site<-data$Site[1:841]


## normalize to a scale of 100
# Calculate the minimum and maximum values of the WQI
min_wqi <- min(wqi$wqi)
max_wqi <- max(wqi$wqi)

# Normalize the WQI to a 0-100 scale
normalized_wqi <- 100 * (wqi$wqi - max_wqi) / (min_wqi - max_wqi)
wqi$wqi<-normalized_wqi


## add to the data so I can look

new_data<-left_join(data, wqi, by=join_by("Site"), copy=TRUE)
head(new_data)
view(new_data)

pad_20_dt<-new_data %>%filter(grepl('_20', Site))

pad_30_dt<-new_data %>%filter(grepl('_30', Site))
pad_18_dt<-new_data %>%filter(grepl('_18', Site))

ggplot(data=pad_18_dt, aes(x=fct_inorder(Site), y=wqi, color=wqi))+
  geom_point(size=3)+
  theme_bw()


##looking at PAd 18 PCA score
pad18_pca<-predict.PCA(lake.pca, data%>%filter(grepl('_18', Site)))
# PCA baseline plot

# Add passive data to plots
plot_18 <- fviz_add(p1_r, pad18_pca$coord,
               geom="point",		# Geometry ("point", "text", "arrows", combinations)
               color="pink",		# Colour
               shape=15,			# Point shape (pch values)
               addlabel=TRUE,		# Label points (T/F)
               labelsize=3,		# Label text size
               repel=TRUE,			# Avoid label overlap
)

##looking closely at PAD18

pad18_wqi<-w1 * pad18_pca$coord[, "Dim.1"] + w2 * pad18_pca$coord[, "Dim.2"]

pad18_wqi<-100 * (pad18_wqi - max_wqi) / (min_wqi - max_wqi)


##looking at PAd 30 PCA score
pad30_pca<-predict.PCA(lake.pca, data%>%filter(grepl('_30', Site)))
# PCA baseline plot

# Add passive data to plots
plot_30 <- fviz_add(p1_r, pad30_pca$coord,
                    geom="point",		# Geometry ("point", "text", "arrows", combinations)
                    color="pink",		# Colour
                    shape=15,			# Point shape (pch values)
                    addlabel=TRUE,		# Label points (T/F)
                    labelsize=3,		# Label text size
                    repel=TRUE,			# Avoid label overlap
)

##looking closely at PAD30

pad30_wqi<-w1 * pad30_pca$coord[, "Dim.1"] + w2 * pad30_pca$coord[, "Dim.2"]
pad30_wqi<-100 * (pad30_wqi - max_wqi) / (min_wqi - max_wqi)


##next up is seperating by the years and plotting the data

mayList <- list()
for (i in 15:19) {
  patterns <- c("M", i, "_")
  seas.pca<-ind.coords[grep(paste(patterns, collapse=""), ind.coords$sites),]
  mayList[[i]] <- seas.pca
}
