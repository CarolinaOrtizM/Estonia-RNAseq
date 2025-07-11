##### Fat content model###########
######Estonian winters 

rm(list=ls(all=TRUE)); # graphics.off()

#working directory
setwd("~/brain/Lipids")

#packages
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)

#data 
Fatt <-read.csv("FatContent.csv", header = T, stringsAsFactors= T)
###subset egg sacs ####
Egg1 = subset(Fatt, Fatt$Eggsac == "1") #there is no difference in the treatment regardless 1st or 2nd egg sac
Egg2 =subset(Fatt, Fatt$Eggsac == "2")

Eggsacs2 <-aov(Egg2$Fatpercentage ~ Egg2$Treatment) 
summary(Eggsacs2)
TukeyHSD(Eggsacs2)

###subset openings 

Onemonth = subset(Fatt, Fatt$Opening == "Dec")
ThreeMonth = subset(Fatt, Fatt$Opening == "Feb") 


###### check the levels ####
levels(Fatt$Treatment)

#rescale the response variable, because it needs to be 0 to 1 
Fatt$Fatpercentage <- Fatt$Fatpercentage/ 100 
Fatt$dryweight
Fatt$Age

#GLMM 
model <- glmmTMB(Fatpercentage ~ Treatment + Opening + dryweight + scale(Age) + (1|ID),
                  data = Fatt,
                  tweedie())


ModelDry <-glmmTMB(dryweight ~ Treatment + Opening + scale(Age) + (1|ID),
                   data = Fatt,
                   gaussian())

# Display the summary of the model
summary(model)

#check assumtions 
#DHARMA test
obj <- simulateResiduals(model, plot = F)
plot(obj, quantreg = F)
Anova(model, type="III")

#post hoc posthoc with Tukey correction for the chosen model 
emm <- emmeans(model_, specs = pairwise ~ Opening | Treatment)
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")

emm <- emmeans(model, specs = pairwise ~  Opening)
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")

###plot the model ###############################################
#22.04.24

model_ <- glmmTMB(Fatpercentage ~ Treatment + Opening,
                 data = Fatt,
                 tweedie())

# Extract unique treatment levels
treatment_levels <- unique(Fatt$Treatment)

# Create a grid of factors for prediction
pred_data <- expand.grid(
  Opening = factor(c("Dec", "Feb"), levels = c("Dec", "Feb")),
  Treatment = factor(c("Cold", "Moderate", "Warm"), levels = c("Cold", "Moderate", "Warm"))
)

# Generate predictions with confidence intervals
predictions <- predict(model_, newdata = pred_data, type = "response")
conf_int <- predict(model_, newdata = pred_data, type = "response", se.fit = TRUE)

# Include the predictions and confidence intervals in pred_data
pred_data$Fatpercentage <- predictions
pred_data$lower <- predictions - 1.96 * conf_int$se.fit  # 95% CI lower bound
pred_data$upper <- predictions + 1.96 * conf_int$se.fit  # 95% CI upper bound
pred_data$DataType <- 'Predicted'

# Transform 'Opening' to numeric in pred_data for consistency
pred_data$Opening_numeric <- as.numeric(pred_data$Opening == "Feb") + 1

# Prepare observed data with necessary columns
observed_data <- Fatt %>%
  mutate(Opening_numeric = case_when(Opening == "Dec" ~ 1, Opening == "Feb" ~ 2),
         DataType = 'Observed',
         lower = NA,  # Not applicable for observed data
         upper = NA)  # Not applicable for observed data

# Ensure both datasets have the same columns, in the same order
if (!"Treatment" %in% names(observed_data)) {
  observed_data$Treatment <- NA  # Assign NA or an appropriate value
}

# Ensure both datasets have identical columns for combining
columns_to_use <- c("Opening_numeric", "Fatpercentage", "DataType", "lower", "upper", "Treatment")

# Add missing columns with NA values where necessary
all_columns <- unique(c(names(observed_data), columns_to_use))
for (col in all_columns) {
  if (!col %in% names(observed_data)) {
    observed_data[[col]] <- NA
  }
  if (!col %in% names(pred_data)) {
    pred_data[[col]] <- NA
  }
}

# Now use rbind with the columns in order
combined_data <- rbind(
  observed_data[columns_to_use],
  pred_data[columns_to_use]
)


# Custom colors
custom_colors <- c("#0072b2", "#009e73", "#d55e00") 

# Plotting the combined data with confidence interval ribbons
Fatt<- ggplot(combined_data, aes(x = Opening_numeric, y = Fatpercentage, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = subset(combined_data, DataType == 'Predicted'),
              aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = subset(combined_data, DataType == 'Observed'),
              aes(x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), # Shift points
              width = 0.05, height = 0, shape = 19, alpha = 0.8) +
  geom_line(data = subset(combined_data, DataType == 'Predicted'), linetype = "solid", size = 1) +
  stat_summary(data = subset(combined_data, DataType == 'Observed'), 
               aes(group = Treatment, x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), 
               fun = median, 
               geom = "point", 
               shape = 22,  # Square shape
               size = 4, 
               position = position_dodge(width = 0.1)) +  # Adjust dodge width if needed
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = c(1, 2), labels = c("1 month", "3 months"))+
  labs(x = NULL, y = "Fat content proportion", color = "Treatment", fill = "Treatment") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black")) +   
  ylim(0.03, 0.16)  # Set y-axis limits from 0 to 0.17
Fatt
ggsave("Fat_content.tiff", width=4, height=6)


#### dry weight only model #####
Fatt$dryweight 
Fatt$dryweight <- Fatt$dryweight/30

# Calculate mean and standard deviation for each Treatment group
summary_stats13 <- Fatt %>%
  group_by(Treatment) %>%
  summarise(mean = mean(dryweight), sd = sd(dryweight))

# Calculate mean plus standard deviation
summary_stats13 <- summary_stats13 %>%
  mutate(upper = mean + sd, lower = mean - sd)

dryF <- ggplot(Onemonth, aes(x = Treatment, y = dryweight, fill= Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.5) +
  geom_point(aes(fill = Treatment), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +
  geom_errorbar(data = summary_stats13, aes(x = Treatment, ymin = lower, ymax = upper),
                width = 0.2, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = summary_stats13, aes(x = Treatment, y = mean, fill = Treatment), shape = 21, size = 5, color = "black") +
  labs(x = NULL, y = "dry weight [mg]", title = NULL) +
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
dryF


Fatt$Fatpercentage <- Fatt$Fatpercentage/30
Fatt$Fatpercentage <- Fatt$Fatpercentage/ 100


model_ <- glmmTMB(dryweight ~ Treatment + Opening,
                  data = Fatt,
                  tweedie())

# Extract unique treatment levels
treatment_levels <- unique(Fatt$Treatment)

# Create a grid of factors for prediction
pred_data <- expand.grid(
  Opening = factor(c("Dec", "Feb"), levels = c("Dec", "Feb")),
  Treatment = factor(c("Cold", "Moderate", "Warm"), levels = c("Cold", "Moderate", "Warm"))
)

# Generate predictions with confidence intervals
predictions <- predict(model_, newdata = pred_data, type = "response")
conf_int <- predict(model_, newdata = pred_data, type = "response", se.fit = TRUE)

# Include the predictions and confidence intervals in pred_data
pred_data$dryweight <- predictions
pred_data$lower <- predictions - 1.96 * conf_int$se.fit  # 95% CI lower bound
pred_data$upper <- predictions + 1.96 * conf_int$se.fit  # 95% CI upper bound
pred_data$DataType <- 'Predicted'

# Transform 'Opening' to numeric in pred_data for consistency
pred_data$Opening_numeric <- as.numeric(pred_data$Opening == "Feb") + 1

# Prepare observed data with necessary columns
observed_data <- Fatt %>%
  mutate(Opening_numeric = case_when(Opening == "Dec" ~ 1, Opening == "Feb" ~ 2),
         DataType = 'Observed',
         lower = NA,  # Not applicable for observed data
         upper = NA)  # Not applicable for observed data

# Ensure both datasets have the same columns, in the same order
if (!"Treatment" %in% names(observed_data)) {
  observed_data$Treatment <- NA  # Assign NA or an appropriate value
}

# Ensure both datasets have identical columns for combining
columns_to_use <- c("Opening_numeric", "dryweight", "DataType", "lower", "upper", "Treatment")

# Add missing columns with NA values where necessary
all_columns <- unique(c(names(observed_data), columns_to_use))
for (col in all_columns) {
  if (!col %in% names(observed_data)) {
    observed_data[[col]] <- NA
  }
  if (!col %in% names(pred_data)) {
    pred_data[[col]] <- NA
  }
}

# Now use rbind with the columns in order
combined_data <- rbind(
  observed_data[columns_to_use],
  pred_data[columns_to_use]
)


# Custom colors
custom_colors <- c("#0072b2", "#009e73", "#d55e00") 

# Plotting the combined data with confidence interval ribbons
dryW<- ggplot(combined_data, aes(x = Opening_numeric, y = dryweight, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = subset(combined_data, DataType == 'Predicted'),
              aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = subset(combined_data, DataType == 'Observed'),
              aes(x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), # Shift points
              width = 0.05, height = 0, shape = 19, alpha = 0.8) +
  geom_line(data = subset(combined_data, DataType == 'Predicted'), linetype = "solid", size = 1) +
  stat_summary(data = subset(combined_data, DataType == 'Observed'), 
               aes(group = Treatment, x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), 
               fun = median, 
               geom = "point", 
               shape = 22,  # Square shape
               size = 4, 
               position = position_dodge(width = 0.1)) +  # Adjust dodge width if needed
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = c(1, 2), labels = c("1 month", "3 months"))+
  labs(x = NULL, y = "dryweight [mg]", color = "Treatment", fill = "Treatment") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))
  #ylim(0.03, 0.16)  # Set y-axis limits from 0 to 0.17
dryW



###### ploting Fat proportion vs FA ####

FattFA <- read.csv("FA_results5.csv", header = T, stringsAsFactors= T)

FattFA$Fatcontent1

ggplot(FattFA, aes(x = TOTAL, y = Fatcontent1)) +  
  geom_point(aes(fill = TREAT), shape = 21, position = position_jitter(width = 0.08), color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = c("#0072b2", "#009e73", "#d55e00")) +  
  geom_smooth(method = "lm") + 
  geom_text(aes(label=Name.FA), vjust=1.5, hjust=1.5) +
  labs(x = "Total FA", y = "Fat content proportion per ind.", title = NULL) +
  theme_minimal() 
  

str(Fatt)
