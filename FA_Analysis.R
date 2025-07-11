####FA Analysis Estonian project ########
#plots and stats

rm(list=ls(all=TRUE)); # graphics.off()

##data with weight correction

FA<-read.csv("FA_results.csv", header = T, stringsAsFactors= T)

#Check the structure of the data and attach 
head(FA)
str(FA)
attach(FA)


############################################ SFA ##########################################################################

FA$Name.FA

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
  labs(x = NULL, y = "SFA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
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
    axis.title.y = element_text(size = 14, color = "black"),
    
  )
SFA

#anova
shapiro.test(FA$SFA)
SFAnova <- aov(FA$SFA ~ FA$TREAT) 
kruskal.test(FA$SFA~ FA$TREAT)
print(SFAnova)
summary(SFAnova)

#Tukey post-hoc 
TukeyHSD(SFAnova)

#anovas of SFA´s

summary(aov(FA$C24.0 ~ FA$TREAT)) 
TukeyHSD(aov(FA$C24.0 ~ FA$TREAT))

tapply(FA$C24.0, FA$TREAT, FUN = mean)
tapply(FA$C24.0, FA$TREAT, FUN = sd)


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
  labs(x = NULL, y = "MUFA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 12, color = "black")
  )
MUFA

#anova
shapiro.test(FA$MUFA)
kruskal.test(FA$MUFA ~ FA$TREAT)
MUFAnova <- aov(FA$MUFA ~ FA$TREAT) 
print(MUFAnova)
summary(MUFAnova)

#Tukey post-hoc 
TukeyHSD(MUFAnova)

#anovas, means and sd of MUFA´s

summary(aov(FA$C24.1 ~ FA$TREAT)) 
TukeyHSD(aov(FA$C24.1 ~ FA$TREAT)) 
 
tapply(FA$MUFA, FA$TREAT, FUN = mean)
tapply(FA$MUFA, FA$TREAT, FUN = sd)


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
  labs(x = NULL, y = "Omega 3 PUFA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
Omega3

#anova
shapiro.test(FA$n.3_PUFA) #not normal 
kruskal.test(FA$n.3_PUFA ~ FA$TREAT)

# Log transformation (add a small constant if there are zeros to avoid -Inf)
FA$log_n.3_PUFA<- log(FA$n.3_PUFA + 1)  # Adding 1 to avoid log(0)
shapiro.test(FA$log_n.3_PUFA ) #normal 

Omega3Anova <- aov(FA$log_n.3_PUFA~ FA$TREAT) 
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

omega3 <- rowSums(FA[, c("C18.4n.3", "C18.3n.3", "C20.3n.3", "C20.5n.3", "C20.4n.3", "C22.5n.3", "C22.6n.3")])

print(omega3)
mean(omega3)

str(FA)

########ALA##########

FA$C18.3n.3
summary_stats7 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C18.3n.3), sd = sd(C18.3n.3))

# Calculate mean plus standard deviation
summary_stats7 <- summary_stats7 %>%
  mutate(upper = mean + sd, lower = mean - sd)

ALA <- ggplot(FA, aes(x = TREAT, y = C18.3n.3, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats7, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats7, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "ALA [ng/µg]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    #axis.text.x =element_text(size = 14, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 12, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
ALA

#normal?
shapiro.test(FA$C18.3n.3 ) #not normal 
kruskal.test(FA$C18.3n.3 ~ FA$TREAT)
#Kruskal-Wallis chi-squared = 9.215, df = 2, p-value = 0.009977
#dunnTest(FA$C18.3n.3 ~ FA$TREAT)  
# C - W 3.005204 0.002654029 0.007962088

#anova
ALAAnova <- aov(FA$C18.3n.3 ~ FA$TREAT) 
#6.814 0.00524**
print(ALAAnova)
summary(ALAAnova)

#Tukey post-hoc 
TukeyHSD(ALAAnova)

tapply(FA$C18.3n.3, FA$TREAT, FUN = mean)
tapply(FA$C18.3n.3, FA$TREAT, FUN = sd)

### log transform

# Log transformation (add a small constant if there are zeros to avoid -Inf)
FA$log_C18.3n.3 <- log(FA$C18.3n.3 + 1)  # Adding 1 to avoid log(0)
shapiro.test(FA$log_C18.3n.3 ) #normal 

# ANOVA with log-transformed data
ALAAnova_log <- aov(log_C18.3n.3 ~ TREAT, data = FA)
print(ALAAnova_log)
summary(ALAAnova_log)

# Tukey post-hoc
TukeyHSD(ALAAnova_log)

ks.test(
  FA$C18.3n.3,  
  "pnorm",  
  mean = mean(FA$C18.3n.3, na.rm = TRUE),  
  sd = sd(FA$C18.3n.3, na.rm = TRUE)
)

######SDA ####

summary_stats16<- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C18.4n.3), sd = sd(C18.4n.3))

# Calculate mean plus standard deviation
summary_stats16 <- summary_stats16 %>%
  mutate(upper = mean + sd, lower = mean - sd)

SDA <- ggplot(FA, aes(x = TREAT, y = C18.4n.3, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats16, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats16, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "SDA [ng/µg]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    #axis.text.x =element_text(size = 14, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 14, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
SDA

#shapiro 
shapiro.test(FA$C18.4n.3 ) # not normal 
kruskal.test(FA$C18.4n.3 ~ FA$TREAT)

# Log transformation (add a small constant if there are zeros to avoid -Inf)
FA$log_C18.4n.3 <- log(FA$C18.4n.3)  # Adding 1 to avoid log(0)
shapiro.test(FA$log_C18.4n.3 ) # normal 

SDAanova <-aov(log_C18.4n.3 ~ TREAT, data = FA)
summary(SDAanova)
TukeyHSD(SDAanova)

####DPA ######

summary_stats17<- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C22.5n.3), sd = sd(C22.5n.3))

# Calculate mean plus standard deviation
summary_stats17 <- summary_stats17 %>%
  mutate(upper = mean + sd, lower = mean - sd)

DPA <- ggplot(FA, aes(x = TREAT, y = C22.5n.3, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1, color = "black") +
  geom_point(data = summary_stats17, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "DPA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 12, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
DPA

shapiro.test(FA$C22.5n.3) #not normal 
kruskal.test(FA$C22.5n.3 ~ FA$TREAT) 

summary(aov(FA$C22.5n.3 ~ FA$TREAT)) 
TukeyHSD(aov(FA$C22.5n.3 ~ FA$TREAT)) 

### log transform

# Log transformation (add a small constant if there are zeros to avoid -Inf)
FA$log_C22.5n.3  <- log(FA$C22.5n.3  + 1)  # Adding 1 to avoid log(0)
shapiro.test(FA$log_C22.5n.3) # not normal 


##### Square root transformation
FA$sqrt_C22.5n.3 <- sqrt(FA$C22.5n.3)
shapiro.test(FA$sqrt_C22.5n.3 ) # not normal 

#### double square root transformation
FA$double_sqrt_C22.5n.3 <- sqrt(sqrt(FA$C22.5n.3))
shapiro.test(FA$double_sqrt_C22.5n.3) # not normal 

FA$INT_C22.5n.3 <- qnorm(
  (rank(FA$C22.5n.3, na.last = "keep") - 0.5) / sum(!is.na(FA$C22.5n.3))
)

shapiro.test(FA$INT_C22.5n.3 )

tapply(FA$C22.5n.3, FA$TREAT, FUN = mean)
tapply(FA$C22.5n.3, FA$TREAT, FUN = sd)

########EPA #######

FA$C20.5n.3
summary_stats11 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C20.5n.3), sd = sd(C20.5n.3))

# Calculate mean plus standard deviation
summary_stats11 <- summary_stats11 %>%
  mutate(upper = mean + sd, lower = mean - sd)

EPA <- ggplot(FA, aes(x = TREAT, y = C20.5n.3, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats11, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats11, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "EPA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
  theme_classic() +
  theme(
    legend.position = "none", 
    #axis.text.x =element_text(size = 14, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 14, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
EPA

shapiro.test(FA$C20.5n.3 ) #normal 
#anova
EPAAnova <- aov(FA$C20.5n.3 ~ FA$TREAT) 
print(EPAAnova)
summary(EPAAnova)

#Tukey post-hoc 
TukeyHSD(EPAAnova)

tapply(FA$C20.5n.3, FA$TREAT, FUN = mean)
tapply(FA$C20.5n.3, FA$TREAT, FUN = sd)

########DHA #######
#FA<-FA[-4, ]

FA$C22.6n.3
summary_stats15 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C22.6n.3), sd = sd(C22.6n.3))

# Calculate mean plus standard deviation
summary_stats15 <- summary_stats15 %>%
  mutate(upper = mean + sd, lower = mean - sd)

DHA <- ggplot(FA, aes(x = TREAT, y = C22.6n.3, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1, color = "black") +
  geom_point(data = summary_stats15, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "DHA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
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
DHA

#anova
shapiro.test(FA$C22.6n.3) #not normal 

DHAAnova <- aov(FA$C22.6n.3 ~ FA$TREAT) 
summary(DHAAnova)

#then a KW test!  
kruskal.test(FA$C22.6n.3 ~ FA$TREAT)
#
dunnTest(FA$C22.6n.3 ~ FA$TREAT)
# M - W -2.967339 0.003003896 0.009011688

#Tukey post-hoc 
TukeyHSD(DHAAnova)

### log transform
# Log transformation (add a small constant if there are zeros to avoid -Inf)
FA$log_C22.6n.3  <- log(FA$C22.6n.3  + 1)  # Adding 1 to avoid log(0)
shapiro.test(FA$log_C22.6n.3  ) # normal 

tapply(FA$C22.6n.3, FA$TREAT, FUN = mean)
tapply(FA$C22.6n.3, FA$TREAT, FUN = sd)

####### DHA/EPA ratio ######
DHA_EPA_ratio <- FA$C22.6n.3 /FA$C20.5n.3

FA$DHA_EPA_ratio <-FA$C22.6n.3 /FA$C20.5n.3

DHA_EPA_ratio_an <- aov(DHA_EPA_ratio ~ FA$TREAT )
summary(DHA_EPA_ratio_an)
TukeyHSD(DHA_EPA_ratio_an)

# Calculate mean and standard deviation for each Treatment group
summary_stats18 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(DHA_EPA_ratio), sd = sd(DHA_EPA_ratio))

# Calculate mean plus standard deviation
summary_stats18 <- summary_stats18 %>%
  mutate(upper = mean + sd, lower = mean - sd)

ratioDHAEPA <- ggplot(FA, aes(x = TREAT, y = DHA_EPA_ratio, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1, color = "black") +
  geom_point(data = summary_stats18, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "Ratio DHA/EPA", title = NULL) +
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
ratioDHAEPA

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
  labs(x = NULL, y = "Omega 6 PUFA [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
Omega6

#anova
shapiro.test(FA$n.6_PUFA) #normal 
Omega6Anova <- aov(FA$n.6_PUFA ~ FA$TREAT) 
print(Omega6Anova)
summary(Omega6Anova)

#Tukey post-hoc 
TukeyHSD(Omega6Anova)

tapply(FA$n.6_PUFA, FA$TREAT, FUN = mean)
tapply(FA$n.6_PUFA, FA$TREAT, FUN = sd)

#anovas, means and sd of Omega 6 PUFA´s
shapiro.test(FA$C22.2n.6)

summary(aov(FA$C22.2n.6 ~ FA$TREAT)) 
TukeyHSD(aov(FA$C22.2n.6 ~ FA$TREAT)) 

tapply(FA$C22.2n.6, FA$TREAT, FUN = mean)
tapply(FA$C22.2n.6, FA$TREAT, FUN = sd)

#####linoleic acid #####

FA$C18.2n.6
summary_stats9 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C18.2n.6), sd = sd(C18.2n.6))

# Calculate mean plus standard deviation
summary_stats9 <- summary_stats9 %>%
  mutate(upper = mean + sd, lower = mean - sd)

LA <- ggplot(FA, aes(x = TREAT, y = C18.2n.6, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats9, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats9, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "LA [ng/µg]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    #axis.text.x =element_text(size = 14, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 14, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
LA


#anova
LA_Anova <- aov(FA$C18.2n.6 ~ FA$TREAT) 
print(LA_Anova)
summary(LA_Anova)

#Tukey post-hoc 
TukeyHSD(LA_Anova)



#####ARA #####
summary_stats10 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(C20.4n.6), sd = sd(C20.4n.6))

# Calculate mean plus standard deviation
summary_stats10 <- summary_stats10 %>%
  mutate(upper = mean + sd, lower = mean - sd)

ARA<- ggplot(FA, aes(x = TREAT, y = C20.4n.6, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats10, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats10, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "ARA [ng/µg]", title = NULL) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x =element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(size = 15, color = "black"),  # Sets y-axis numbers to bigger and black
    axis.title.y = element_text(size = 15, color = "black")) +
  scale_x_discrete(labels = c("C" = "Cold", "M" = "Moderate", "W" = "Warm"))
ARA

#normal?
 shapiro.test(FA$C20.4n.6) #normal 

#anova
ARA_Anova <- aov(FA$C20.4n.6 ~ FA$TREAT) 
print(ARA_Anova)
summary(ARA_Anova)

#Tukey post-hoc 
TukeyHSD(ARA_Anova)

tapply(FA$C20.4n.6, FA$TREAT, FUN = mean)
tapply(FA$C20.4n.6, FA$TREAT, FUN = sd)

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
  labs(x = NULL, y = "Total Fatty acids [ng/µg]", title = NULL) +
  #geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5)+
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
shapiro.test(FA$TOTAL) #normal
TOTAnova <- aov(FA$TOTAL ~ FA$TREAT) 
print(TOTAnova)
summary(TOTAnova)

#Tukey post-hoc 
TukeyHSD(TOTAnova)

#To fill the table of FA 
tapply(FA$TOTAL, FA$TREAT, FUN = mean)
tapply(FA$TOTAL, FA$TREAT, FUN = sd)

######################Omega 3 / Omega 6 ratio ##################################################
n3_n6_ratio <- FA$n.3_PUFA /FA$n.6_PUFA

#add it to the data frame 
FA$ratio_n3_n6 <-FA$n.3_PUFA /FA$n.6_PUFA

n3_n6_ratio_an <- aov(n3_n6_ratio ~ FA$TREAT )
summary(n3_n6_ratio_an)
TukeyHSD(n3_n6_ratio_an)

# Calculate mean and standard deviation for each Treatment group
summary_stats6 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(ratio_n3_n6), sd = sd(ratio_n3_n6))

# Calculate mean plus standard deviation
summary_stats6 <- summary_stats6 %>%
  mutate(upper = mean + sd, lower = mean - sd)

ratio <- ggplot(FA, aes(x = TREAT, y = ratio_n3_n6, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats6, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats6, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "Ratio n3/n6", title = NULL) +
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
ratio

#To fill the table of FA 
tapply(FA$ratio_n3_n6, FA$TREAT, FUN = mean)


PUFAratioAnova <- aov(FA$ratio_n3_n6 ~ FA$TREAT) 
print(PUFAratioAnova)
summary(PUFAratioAnova)
TukeyHSD(PUFAratioAnova)

###########################Unsaturated (MUFAS + PUFAS) / Saturated ratio (SFA) ###########################################

U_S_ratio <- (FA$MUFA + FA$n.3_PUFA + FA$n.6_PUFA) /FA$SFA

#add it to the data frame 
FA$U_S_ratio <-(FA$MUFA + FA$n.3_PUFA + FA$n.6_PUFA) /FA$SFA


tapply(FA$U_S_ratio, FA$TREAT, FUN = mean)
tapply(FA$U_S_ratio, FA$TREAT, FUN = sd)


U_S_ratio_an <- aov(U_S_ratio ~ FA$TREAT )
summary(U_S_ratio_an)
TukeyHSD(U_S_ratio_an)

summary_stats8 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(U_S_ratio), sd = sd(U_S_ratio))

# Calculate mean plus standard deviation
summary_stats8 <- summary_stats8 %>%
  mutate(upper = mean + sd, lower = mean - sd)


US <- ggplot(FA, aes(x = TREAT, y = U_S_ratio, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats8, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats8, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "U/S ratio ", title = NULL) +
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


###########plot all together ######

final_grid <- grid.arrange(SFA, MUFA, ALA, LA, EPA, ARA, 
                           ncol = 2, nrow = 3,
                           widths = c(1, 1),  # Equal column widths
                           heights = c(1, 1,1),     # Equal row heights
                           padding = unit(0.5, "line"))  # Reduce space between plots

suplement_grid <- grid.arrange(Omega3, Omega6, 
                           ncol = 2, nrow = 1,
                           widths = c(1, 1),  # Equal column widths
                           #heights = c(1, 1),     # Equal row heights
                           padding = unit(0.5, "line"))  # Reduce space between plots

ratios_grid<- grid.arrange(US, ratio,  
                        ncol = 2, nrow = 1,
                        widths = c(1, 1),  # Equal column widths
                        #heights = c(1, 1),     # Equal row heights
                        padding = unit(0.5, "line"))  # Reduce space between plots

omega3s_grid <- grid.arrange(SDA, DPA, DHA, 
                             ncol = 3, nrow = 1,
                             #widths = c(1, 1, 1),  # Equal column widths
                             #heights = c(1, 1, 1),     # Equal row heights
                             padding = unit(0.5, "line"))  # Reduce space between plots

grid8 <- grid.arrange(SFA, ALA, MUFA, SDA, LA, EPA, ARA, DPA, 
                      ncol = 2, nrow = 4,
                      widths = c(1, 1),  # Equal column widths
                      heights = c(1, 1,1,1))     # Equal row heights
                      #padding = unit(0.5, "line"))  # Reduce space between plots

### dry weight 

FA$dryweight1mg

# Calculate mean and standard deviation for each Treatment group
summary_stats12 <- FA %>%
  group_by(TREAT) %>%
  summarise(mean = mean(dryweight1mg), sd = sd(dryweight1mg))

# Calculate mean plus standard deviation
summary_stats12 <- summary_stats12 %>%
  mutate(upper = mean + sd, lower = mean - sd)

dry <- ggplot(FA, aes(x = TREAT, y = dryweight1mg, fill= TREAT)) +
  geom_violin(aes(fill = TREAT), alpha = 0.5) +
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats12, aes(x = TREAT, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats12, aes(x = TREAT, y = mean, fill = TREAT), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "dry weight [µg]", title = NULL) +
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
dry

#anova
dryAnova <- aov(FA$dryweight1mg ~ FA$TREAT) 
print(dryAnova)
summary(dryAnova)

#Tukey post-hoc 
TukeyHSD(dryAnova)









