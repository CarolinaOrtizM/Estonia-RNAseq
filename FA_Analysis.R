####FA Analysis Estonian project ########
#plots and stats

rm(list=ls(all=TRUE)); # graphics.off()

#data 
FA <-read.csv("FA_results.csv", header = T, stringsAsFactors= T)

FA <- FA[-25,] #take the last NA row 

#Check the structure of the data and attach 
head(FA)
str(FA)
attach(FA)

#packages
library("ggplot2")
library(gridExtra)

#plots of the main groups 

############################################ SFA ##########################################################################

FA$SFA <-  FA [,"SFA"]   / 1 


# Calculate mean and standard deviation for each Treatment group
summary_stats1 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(SFA), sd = sd(SFA))

# Calculate mean plus standard deviation
summary_stats1 <- summary_stats1 %>%
  mutate(upper = mean + sd, lower = mean - sd)


SFA <- ggplot(FA, aes(x = TREAT, y = SFA, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats1, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats1, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "SFA [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black"),
    
  )
SFA

#anova
SFAnova <- aov(SFA ~ TREAT) 
print(SFAnova)
summary(SFAnova)

#Tukey post-hoc 
TukeyHSD(SFAnova)

#anovas of SFA´s
C24<-aov(FA$C24.0 ~ FA$TREAT) 
summary(C24)
TukeyHSD(C24)

tapply(FA$SFA, FA$TREAT, FUN = mean)
tapply(FA$SFA, FA$TREAT, FUN = sd)

################################### MUFA ##########################################

# Calculate mean and standard deviation for each Treatment group
summary_stats2 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(MUFA), sd = sd(MUFA))

# Calculate mean plus standard deviation
summary_stats2 <- summary_stats2 %>%
  mutate(upper = mean + sd, lower = mean - sd)


MUFA <- ggplot(FA, aes(x = TREAT, y = MUFA, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats2, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats2, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "MUFA [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")
  )
MUFA

#anova
MUFAnova <- aov(FA$MUFA ~ FA$TREAT) 
print(MUFAnova)
summary(MUFAnova)

#Tukey post-hoc 
TukeyHSD(MUFAnova)

#anovas, means and sd of MUFA´s

tapply(FA$MUFA, FA$TREAT, FUN = mean)
tapply(FA$MUFA, FA$TREAT, FUN = sd)

C22_1_9 <-aov(FA$C22.1n.9 ~ FA$TREAT) 
summary(C22_1_9)
TukeyHSD(C22_1_9)

################################### Omega n-3 #############################################
FA$n.3_PUFA

# Calculate mean and standard deviation for each Treatment group
summary_stats3 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(n.3_PUFA), sd = sd(n.3_PUFA))

# Calculate mean plus standard deviation
summary_stats3 <- summary_stats3 %>%
  mutate(upper = mean + sd, lower = mean - sd)


Omega3 <- ggplot(FA, aes(x = TREAT, y = n.3_PUFA, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats3, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats3, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "Omega 3 PUFA [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
    scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
Omega3

#anova
Omega3Anova <- aov(FA$n.3_PUFA ~ FA$TREAT) 
print(Omega3Anova)
summary(Omega3Anova)

#Tukey post-hoc 
TukeyHSD(Omega3Anova)

#anovas, means and sd of Omega 3 PUFA´s

tapply(FA$n.3_PUFA, FA$TREAT, FUN = mean)
tapply(FA$n.3_PUFA, FA$TREAT, FUN = sd)

C22_6_3 <-aov(FA$C22.6n.3 ~ FA$TREAT) 
summary(C22_6_3)
TukeyHSD(C22_6_3)

FA$C22.6n.3

omega3 <- rowSums(FA[, c("C18.4n.3", "C18.3n.3", "C20.3n.3", "C20.5n.3", "C20.4n.3", "C22.5n.3", "C22.6n.3")])

print(omega3)
mean(omega3)

str(FA)

################################### Omega n-6  #############################################
FA$n.6_PUFA

# Calculate mean and standard deviation for each Treatment group
summary_stats4 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(n.6_PUFA), sd = sd(n.6_PUFA))

# Calculate mean plus standard deviation
summary_stats4 <- summary_stats4 %>%
  mutate(upper = mean + sd, lower = mean - sd)


Omega6 <- ggplot(FA, aes(x = TREAT, y = n.6_PUFA, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats4, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats4, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "Omega 6 PUFA [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
Omega6

#anova
Omega6Anova <- aov(FA$n.6_PUFA ~ FA$TREAT) 
print(Omega6Anova)
summary(Omega6Anova)

#Tukey post-hoc 
TukeyHSD(Omega6Anova)

#anovas, means and sd of Omega 6 PUFA´s

tapply(FA$n.6_PUFA, FA$TREAT, FUN = mean)
tapply(FA$n.6_PUFA, FA$TREAT, FUN = sd)

C22_5_6 <-aov(FA$C22.5n.6 ~ FA$TREAT) 
summary(C22_5_6)
TukeyHSD(C22_5_6)


### plot all together...

firstgrid <- grid.arrange(SFA, MUFA, Omega3, Omega6,  ncol = 2, nrow = 2)

####### Total FA ##########################################################################################
FA$TOTAL

# Calculate mean and standard deviation for each Treatment group
summary_stats5 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(TOTAL), sd = sd(TOTAL))

# Calculate mean plus standard deviation
summary_stats5 <- summary_stats5 %>%
  mutate(upper = mean + sd, lower = mean - sd)

TOT <- ggplot(FA, aes(x = TREAT, y = TOTAL, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats5, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats5, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "Total Fatty acids [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
TOT

#anova
TOTAnova <- aov(FA$TOTAL ~ FA$TREAT) 
print(TOTAnova)
summary(TOTAnova)

#Tukey post-hoc 
TukeyHSD(TOTAnova)

#To fill the table of FA 
tapply(FA$TOTAL, FA$TREAT, FUN = mean)
tapply(FA$TOTAL, FA$TREAT, FUN = sd)

################################Omega 3 / Omega 6 ratio #############################################################################

n3_n6_ratio <- FA$n.3_PUFA /FA$n.6_PUFA

#add it to the data frame 
FA$ratio_n3_n6 <-FA$n.3_PUFA /FA$n.6_PUFA

n3_n6_ratio_an <- aov(n3_n6_ratio ~ FA$TREAT )
summary(n3_n6_ratio_an)
TukeyHSD(n3_n6_ratio_an)

#To fill the table of FA 
tapply(FA$ratio_n3_n6, FA$TREAT, FUN = mean)

#####################################################Unsaturated (MUFAS + PUFAS) / Saturated ratio (SFA) #######################################################
#

U_S_ratio <- (FA$MUFA + FA$n.3_PUFA + FA$n.6_PUFA) /FA$SFA

#add it to the data frame 
FA$U_S_ratio <-(FA$MUFA + FA$n.3_PUFA + FA$n.6_PUFA) /FA$SFA


tapply(FA$U_S_ratio, FA$TREAT, FUN = mean)
tapply(FA$U_S_ratio, FA$TREAT, FUN = sd)


U_S_ratio_an <- aov(U_S_ratio ~ FA$TREAT )
summary(U_S_ratio_an)
TukeyHSD(U_S_ratio_an)

summary_stats7 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(U_S_ratio), sd = sd(U_S_ratio))

# Calculate mean plus standard deviation
summary_stats7 <- summary_stats7 %>%
  mutate(upper = mean + sd, lower = mean - sd)

US <- ggplot(FA, aes(x = TREAT, y = U_S_ratio, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats7, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats7, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "U/S [ng per ind.]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
US






