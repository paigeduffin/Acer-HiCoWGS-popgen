# Author: Paige Duffin
# Isolation by distance using neutral-panel sampling-location FST.
# Inputs (in the working directory):
# - Acer.haplo.locs_ABREV_w.panama.csv: first two columns longitude/latitude;
#   Name contains location labels matching the FST table.
# - STAMPP_fst.bootstrap_samploc.NEU.pan_178k_4.IBD.csv: popA, popB, fst columns.
# The original Mantel matrix construction and test are retained.

library(geosphere)
library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ade4)

coords_file <- "Acer.haplo.locs_ABREV_w.panama.csv"
fst_file <- "STAMPP_fst.bootstrap_samploc.NEU.pan_178k_4.IBD.csv"
out_dir <- "ibd_results"
if (file.exists(out_dir)) stop("Output path already exists: ", out_dir)
dir.create(out_dir)

# Geographic distances and unique location pairs.
coords <- read.csv(coords_file)
rownames(coords) <- coords$Name
pw_dist <- distm(coords[, 1:2], coords[, 1:2], fun = distCosine)
rownames(pw_dist) <- rownames(coords)
colnames(pw_dist) <- rownames(coords)

mat <- as.matrix(pw_dist)
haversine_dist <- as.data.frame(as.table(mat))
haversine_dist_unique <- subset(
    haversine_dist, as.numeric(row(mat)) <= as.numeric(col(mat))
)
colnames(haversine_dist_unique) <- c("popA", "popB", "havers_geo")
haversine_dist_unique <- haversine_dist_unique[
    haversine_dist_unique$popA != haversine_dist_unique$popB, ]

# Match pairwise FST and geographic distances regardless of pair orientation.
weir_fst <- read_csv(fst_file)
haversine_dist_unique$popA <- as.character(haversine_dist_unique$popA)
haversine_dist_unique$popB <- as.character(haversine_dist_unique$popB)
df1_norm <- haversine_dist_unique %>%
    mutate(pop1 = pmin(popA, popB), pop2 = pmax(popA, popB))
df2_norm <- weir_fst %>%
    mutate(pop1 = pmin(popA, popB), pop2 = pmax(popA, popB))
havers.fst_merged <- df1_norm %>%
    inner_join(df2_norm, by = c("pop1", "pop2"), suffix = c(".df1", ".df2"))
havers.fst_merged.clean <- havers.fst_merged %>%
    transmute(popA = pop1, popB = pop2, havers_geo, fst,
              fst.over.1.minus.fst = fst / (1 - fst))

# Original palettes for the two plotting views.
colors_popA_specific <- c(
    "#04A2A0", "#AA36A5", "#699213", "#C6CF18", "#FBCE46",
    "#CF550A", "#7A3DA6", "#4E66C2", "#B42D2D"
)
colors_popB_specific <- c(
    "#AA36A5", "#699213", "#C6CF18", "#FBCE46", "#CF550A",
    "#7A3DA6", "#4E66C2", "#B42D2D", "#EE962F"
)

IBD_plot_popB <- ggplot(havers.fst_merged.clean,
                        aes(x = havers_geo / 1000, y = fst.over.1.minus.fst)) +
    geom_point(aes(color = popB)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey", size = 0.5) +
    xlab("Haversine distance (km)") + ylab("Fst / (1 - Fst)") +
    scale_color_manual(values = colors_popB_specific)

IBD_plot_popA <- ggplot(havers.fst_merged.clean,
                        aes(x = havers_geo / 1000, y = fst.over.1.minus.fst)) +
    geom_point(aes(color = popA)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey", size = 0.5) +
    xlab("Haversine distance (km)") + ylab("Fst / (1 - Fst)") +
    scale_color_manual(values = colors_popA_specific)

pdf(file.path(out_dir, "ibd_plots.pdf"), width = 7, height = 5)
print(IBD_plot_popB)
print(IBD_plot_popA)
dev.off()

# Mantel step: original code, including zero filling and dist(), unchanged.
IBDma=dcast(havers.fst_merged.clean,popA ~ popB, value.var ='fst.over.1.minus.fst')
rownames(IBDma)=IBDma$popA
IBDma$popA=NULL
IBDma[is.na(IBDma)]=0
Distma=pw_dist[match(rownames(IBDma),rownames(pw_dist)),match(colnames(IBDma),colnames(pw_dist))]

mantest=mantel.rtest(dist(IBDma),dist(Distma),nrepet=100000)

print(mantest)
capture.output(print(mantest), file = file.path(out_dir, "mantel_test.txt"))
saveRDS(mantest, file.path(out_dir, "mantel_test.rds"))
write.csv(havers.fst_merged.clean, file.path(out_dir, "ibd_pairs.csv"),
          row.names = FALSE)
