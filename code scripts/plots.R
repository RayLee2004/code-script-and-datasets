# R version: R-4.5.0

########## PART 1: PCoA ##########

rm(list=ls())
library(vegan)
library(ggplot2)
library(patchwork)

groups <- read.table('sample file/pcoa_gorups (pcoa_example).txt', header = FALSE, colClasses=c("character","character"))
colnames(groups) <- c("sample", "group")

length <- length(unique(groups$sample))
length1 <- length(unique(groups$group))
times1 <- length%/%8
res1 <- length%%8
times2 <- length%/%5
res2 <- length%%5

mycol = c("#647ADD","#FFDF67","#F5664D")
col1 = rep(mycol, times1)
col = c(col1, mycol[1:res1])
pich1 <- rep(c(15:18,20,7:14,0:6), times2)
pich <- c(pich1, 15:(15+res2))
n = ifelse(length1 > 30, 2, 1)

# Please replace the datasets here. 
data <- read.table('merged_multigroups (pcoa_example).txt', header = TRUE, 
                   row.names = 1, sep = "\t", comment.char = "", check.names = FALSE)
data <- t(data)

distance <- vegdist(data, method = 'bray')

pcoa <- cmdscale(distance, k = 5, eig = TRUE)
points <- data.frame(scores(pcoa))
point <- data.frame(sample = row.names(points), points) 

pc1 <- round((pcoa$eig/sum(pcoa$eig))*100,2)[1]
pc2 <- round((pcoa$eig/sum(pcoa$eig))*100,2)[2]

colnames(points)[1:2] <- c('dim1', 'dim2')
plotdata = data.frame(rownames(points), points$dim1, points$dim2, groups$group)
colnames(plotdata) = c("sample","dim1","dim2","group")

adonis_result_dis <- adonis2(distance ~ group, data = groups, permutations = 999)
R2 = adonis_result_dis$R2[1]
pvalue = adonis_result_dis$`Pr(>F)`[1]
adonis <- paste("PERMANOVA:\nR2 = ", round(R2,4), "\nP-value = ", pvalue)

plot <- ggplot(plotdata, aes(dim1, dim2)) + 
  geom_point(aes(colour = group), size = 7) + 
  scale_shape_manual(values = pich1) + 
  scale_colour_manual(values = col) +
  xlab(paste("PCoA axis1 ( ",pc1,"%"," )", sep = "")) + 
  ylab(paste("PCoA axis2 ( ",pc2,"%"," )", sep = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        panel.grid = element_blank(), 
        axis.title = element_text(color = 'black', size = 20),
        axis.text = element_text(colour = 'black', size = 16),
        axis.ticks = element_line(color = 'black'), 
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key = element_blank(),
        legend.position = c(0.1, 0.9),
        legend.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.25)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  guides(col = guide_legend(ncol = n), shape = guide_legend(ncol = n)) + 
  stat_ellipse(aes(color = group))

box1 <- ggplot(plotdata, aes(x = groups$group, y = dim1, fill = groups$group)) +
  geom_boxplot(show.legend = FALSE) + 
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
  geom_jitter(show.legend = FALSE) + 
  scale_fill_manual(values = col) +
  coord_flip() + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.25))

box2 <- ggplot(plotdata, aes(x = groups$group, y = dim2, fill = groups$group)) + 
  geom_boxplot(show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
  geom_jitter(show.legend = FALSE) + 
  scale_fill_manual(values = col) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(colour ='black', size = 16, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.25))

box3 <- ggplot(plotdata, aes(dim1, dim2)) +
  geom_text(aes(x = -0.5, y = 0.6, label = adonis), size = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.25))

p <- box1 + box3 + plot + box2 + plot_layout(heights = c(1,4), widths = c(4,1), ncol = 2, nrow = 2)

# Please replace the name of the output image here. 
png("multigroups.png", res = 600, height = 5400, width = 5400, type = "cairo")

dev.off()


########## PART 2 Procrustes Analysis ##########

rm(list = ls())
library(openxlsx)
library(vegan)
library(ggplot2)

# Please replace the datasets of species and environment factors (or another species)
species <- read.csv("sample file/fungi_2024 (procrustes_example).csv", header = T, row.names = 1)
env <- read.csv("sample file/fish_2024 (procrustes_example).csv", header = T, row.names = 1)

spe.dist <- vegdist(species)
env.dist <- vegdist(scale(env), "euclid")

set.seed(123)
mds.s <- monoMDS(spe.dist)
mds.e <- monoMDS(env.dist)
pro.s.e <- procrustes(mds.s, mds.e, symmetric = TRUE)

pro.s.e_t <- protest(mds.s, mds.e, permutations = 999)

Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

ggplot(data = Pro_Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
               arrow = arrow(length = unit(0, 'cm')), color = "#E41A1C", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.2, 'cm')), color = "#377EB8", size = 1) +
  geom_point(aes(X1, X2), color = "#E41A1C", size = 8, shape = 16) +
  geom_point(aes(MDS1, MDS2), color = "#377EB8", size = 8, shape = 16) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[1, 2]/Pro_X[1, 1], size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[2, 2]/Pro_X[2, 1], size = 0.3)

# Please replace the name of the output image here. 
ggsave("fungi-fish.png", height = 5, width = 5)


########## PART 3 Module Compensation Effect Analysis ##########

library(ggplot2)
library(openxlsx)

# Please replace the datasets of mean value and (75% Confidence Interval) Range here.
df_ave <- read.xlsx("sample file/mean_path_length (compensation_example).xlsx", 1)
df_range <- read.xlsx("sample file/mean_path_length (compensation_example).xlsx", 2)

module_levels <- unique(df_ave$module)
df_ave$module_id <- as.numeric(factor(df_ave$module, levels = module_levels))
df_range$module_id <- as.numeric(factor(df_range$module, levels = module_levels))

ggplot() +
  geom_ribbon(data = df_range, aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)), alpha = 0.12) +
  geom_line(data = df_ave, aes(x = module_id, y = value, color = factor(year)), size = 1) +
  scale_x_continuous(breaks = 1:length(module_levels), labels = module_levels) +
  scale_fill_manual(values = c("2022" = "#647ADD", "2023" = "#FFDF67", "2024" = "#F5664D")) +
  scale_color_manual(values = c("2022" = "#647ADD", "2023" = "#FFDF67", "2024" = "#F5664D")) +
  facet_wrap(~year, nrow = 3, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7.5, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 1),
    axis.title = element_blank(),
    plot.title = element_text(size = 34, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "none"
  )

# Please replace the name of the output image here. 
ggsave("mean_path_length_compensation_effect.png", height = 5, width = 5)
