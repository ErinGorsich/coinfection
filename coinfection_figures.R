library(ggplot2)
library(lattice)

###################################################
###################################################
# Figure 1
###################################################
###################################################
df <- data.frame(
    label <- c("Rift Valley fever","parainfluenza", "brucellosis", 
        "bovine respitory \n syncytial virus", "bovine viral \n diarrhea virus", 
        "bovine herpes virus", "A. marginale", "A. centrale"),
    infection <- c("RVF","PI3", "Brucellosis", "BRSV", "BVDV", "BHV", 
        "A. marginale", "A. centrale"),
    change <- c(50, 8, 5, 3, 0, 0, 0, -3), 
    order <- seq(1, 8, 1), 
    incidence <- c("yes", "yes", "sometimes", "no", "unknown", "no", "unknown", 
        "unknown") )

ggplot(data = df, aes(x = label, y = change)) + 
    geom_bar(stat = "identity", colour = "black", cex = 2) + 
    coord_flip() + theme_bw() + 
    ylab("Change in prevalence \n prevalence in bTB+ - prevalence in bTB-") +
    #geom_errorbar(limits, position = position_dodge(width = 0.9), 
    #    width= 0.2, cex = 2 )
    theme(
        panel.border = element_blank(), 
        panel.margin = element_blank(),
        #panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.line.x = element_line(colour= "black", size = 1), 
        axis.line.y = element_line(colour= "black", size = 1),
        axis.title.y = element_blank())

###################################################
###################################################
# FIGURE 2: 
###################################################
###################################################
setwd("~/Documents/collaborations/Bree-Carrie-BTBcoinfection_paper/code/24-May-2017 update/R figures")
source('~/git/coinfection/fixedparams.R', chdir = TRUE)
source('~/git/coinfection/acutemildparams.R', chdir = TRUE)
source('~/git/coinfection/chronicmildparams.R', chdir = TRUE)
add <- read.csv("acute_densitydependent_recov_2019.csv")
add <- add[add$type == "variable", ]
add$recov <- 1/add$change_gamma_tu
cdd <- read.csv("chronic_densitydependent_2019.csv")
cdd <- cdd[cdd$type == "variable", ]
addm <- read.csv("acute_densitydependent_2019.csv")
addm <- addm[addm$type == "variable", ]

# p <- ggplot(data = add, aes(x = gamma_tu, y = beta_tu, fill = EE_R)) +
#     geom_raster(interpolate = TRUE) +
#     theme_bw() + 
#     xlab("Proportional increase in infection duration with co-infection") + 
#     ylab("Proportional increase in transmission with co-infection") +
#     #geom_title("Acute pathogen; density-dependent transmission") +
#     theme(panel.grid.major = element_blank(), 
#         plot.title = element_text(size = 14),  # 14 in pdf format
#         axis.title = element_text(size = 12), 
#         axis.text = element_text(size = 12), # 12 
#         legend.text = element_text(size=12),  # 12 
#         legend.title = element_text(size = 14),
#         legend.text.align = 1, 
#         legend.key.size = unit(1, "cm") )

cols <- colorRampPalette(brewer.pal(9, "Greens"))(50)

ylim = c(min(c(addm$RoAT, add$RoAT, cdd$RoC)),
    max(c(addm$RoAT, add$RoAT, cdd$RoC))+0.1)
#ylim = c(1, 7)

p <- levelplot(RoAT~ recov*change_beta_tu, data= add,
    at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20),
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in infection duration,", 
        gamma[a]/gamma[ta])), 
    col.regions=cols, main="Acute", cex.lab=1.2, 
    colorkey = FALSE)
plot(p, position = c(0, 0, 1/3, 1), more = TRUE) 
p2 <- levelplot(RoAT~ change_alpha_tu*change_beta_tu, 
    data= addm, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Acute", cex.lab = 1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), 
    colorkey = FALSE)
plot(p2, position = c(1/3, 0, 2/3, 1), more = TRUE) 
p3 <- levelplot(RoCT~ change_alpha_tc*change_beta_tc, 
    data= cdd, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Chronic", cex.lab=1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), colorkey = FALSE )
plot(p3, position = c(2/3, 0, 1, 1), more = FALSE) 
legend <- levelplot(RoCT~ change_alpha_tc*change_beta_tc, 
    data= cdd, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Chronic", cex.lab=1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), colorkey = TRUE ) # 1000*500

# EE
ylim = c(min(c(addm$EE_R, add$EE_R, cdd$EE_C)),
    max(c(addm$EE_R, add$EE_R, cdd$EE_C))+0.1)
p <- levelplot(EE_R~ recov*change_beta_tu, data= add,
    at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20),
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in infection duration,", 
        gamma[a]/gamma[ta])), 
    col.regions=cols, main="Acute", cex.lab=1.2, 
    colorkey = FALSE)
plot(p, position = c(0, 0, 1/3, 1), more = TRUE) 
p2 <- levelplot(EE_R~ change_alpha_tu*change_beta_tu, 
    data= addm, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Acute", cex.lab = 1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), 
    colorkey = FALSE)
plot(p2, position = c(1/3, 0, 2/3, 1), more = TRUE) 
p3 <- levelplot(EE_C~ change_alpha_tc*change_beta_tc, 
    data= cdd, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Chronic", cex.lab=1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), colorkey = FALSE )
plot(p3, position = c(2/3, 0, 1, 1), more = FALSE)  # 1000*500

legend <- levelplot(EE_C~ change_alpha_tc*change_beta_tc, 
    data= cdd, at = seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/20), 
    col.regions=cols, main="Chronic", cex.lab=1.2, 
    ylab = expression(paste("Proportional change in transmission, ",
        beta[tc]/beta[c])),
    xlab = expression(paste("Proportional change in mortality,", 
        alpha[ta]/alpha[a])), colorkey = TRUE )
