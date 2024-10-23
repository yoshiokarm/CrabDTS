#  Crab DTS ANALYSES--------------------------------------------------------------
#  code by Reyn Yoshioka
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# Analyses for Yoshioka and Groner (YEAR)
# Ambient-temperature preservation of fatty acids in marine crabs for remote 
# field collections

# 1: Getting Started -----------------------------------------------------------
# > 1.1: Clear everything ----- 
# clear workspace
# rm(list=ls()) # Currently commented out to prevent folks from clearing things
# clear console
cat("\014")

# > 1.2: Load packages -----
library(here) # working directory management
library(tidyverse) # all the tidy things
library(grid)
library(scales)
library(ggrepel)
library(gridExtra)
library(viridis) # colorblind-friendly and perceptually uniform plotting
library(glmmTMB) # mixed effects glms 
library(DHARMa) # GLM residual diagnostics
library(easyCODA) # Compositional data, needed for G&G logratio approach
library(patchwork) # combining plots
library(stringr) # working with strings

# > 1.3: Universal settings -----
# Plot resolution
res = 300

# > 1.3: Here ----
here::i_am("CrabDTS_Analysis.R")

# > 1.4: Read in data -----
# Sample log and data related to the crab
df_samplog = read.csv("CrabDTS_SampleLog.csv")
df_samplog$Row = NULL
df_samplog$X = NULL

# Fatty acid data
df_full = read.csv("CrabDTS_FA_data.csv")
df_full$X = NULL
df_Copeman = read.csv("CrabDTS_Copeman_etal_2021.csv")

# > 1.5: Add columns -----
# > > 1.5.1: Quantity Notes -----
df_full$notes_quant = ifelse(df_full$ug == "<MRL",
                             "<MRL",
                             ifelse(df_full$ug == "<LOQ",
                                    "<LOQ",
                                    ifelse(df_full$ug == "<LOD",
                                           "<LOD",
                                           NA)))
df_full$ug = as.numeric(df_full$ug)

# > > 1.5.2: Sample Information -----
df_full = left_join(df_full,
                    df_samplog,
                    by = join_by(SampleID))

# > > 1.5.3: Zeroes to <LOD, <MRL, <LOQ -----
df_full$ug0 = df_full$ug
df_full[is.na(df_full$ug0), ]$ug0 = 0

# > 1.6: Data organization and resplitting -----
# > > 1.6.1: Separate Phases -----
# Select lab validation data 
df_LV = df_full[df_full$Set =="Lab Validation",]
df_LV = df_LV[!is.na(df_LV$FA),]

# Select field validation data
df_FV = df_full[df_full$Set == "Field Validation",]
df_FV = df_FV[!is.na(df_FV$FA),]

# > > 1.6.2:  Reattach total FA rows as columns -----
# for lab validation
banana = df_LV[df_LV$FA == "TFA_ug",]
banana$TFA_ug = banana$ug0
df_LV = left_join(df_LV,
                  banana[,c("SampleID", "TFA_ug")])
rm(banana)
df_LV = df_LV[df_LV$FA != "TFA_ug", ]

# for field validation
banana = df_FV[df_FV$FA == "TFA_ug",]
banana$TFA_ug = banana$ug0
df_FV = left_join(df_FV,
                  banana[,c("SampleID", "TFA_ug")])
rm(banana)
df_FV = df_FV[df_FV$FA != "TFA_ug", ]

# > > 1.6.3:  Calculate percentages -----
df_LV$perc = df_LV$ug0 / df_LV$TFA_ug * 100
df_FV$perc = df_FV$ug0 / df_FV$TFA_ug * 100

# > >  1.6.5: Separate Tissues -----
df_LV_HP = df_LV[df_LV$Tissue == "HP",]
df_LV_HL = df_LV[df_LV$Tissue == "HL",]

df_FV_HP = df_FV[df_FV$Tissue == "HP",]
df_FV_HL = df_FV[df_FV$Tissue == "HL",]

# > 1.7: Lookups and Labels ----- 
labels_logratios = data.frame(y = c(0.005,
                                    0.01,
                                    0.02,
                                    0.05,
                                    0.1,
                                    0.2,
                                    0.5,
                                    1,
                                    2,
                                    5,
                                    10,
                                    20,
                                    50,
                                    100,
                                    200),
                              label = c("1:200",
                                        "1:100",
                                        "1:50",
                                        "1:20",
                                        "1:10",
                                        "1:5",
                                        "1:2",
                                        "1:1",
                                        "2:1",
                                        "5:1",
                                        "10:1",
                                        "20:1",
                                        "50:1",
                                        "100:1",
                                        "200:1"))
labels_logratios$y = log(labels_logratios$y)

FA_labeller = c(
  C16.1n7c = "Palmitoleic",
  C18.1n9c = "Oleic",
  C18.2n6c = "LIN",
  C18.3n3c = "ALA",
  C20.4n6c = "ARA",
  C20.5n3 = "EPA",
  C22.6n3 = "DHA",
  C22.4n6c = "Adrenic",
  C14.0 = "Myristic",
  C16.0 = "Palmitic",
  C18.0 = "Stearic",
  C18.4n3 = "SDA",
  C22.5n3 = "DPA",
  C20.1n9c = "20:1n9",
  C22.1n9c = "22:1n9",
  C17.0 = "17:0"
)

TypeMo_labeller = c(
  DTS_2 = "DTS, 2 mo",
  DTS_3 = "DTS, 3 mo",
  DTS_5 = "DTS, 5 mo",
  Tissue_5 = "Frozen, 5 mo"
)

sort(unique(df_full$FA))
OBCFA_list = c("C11.0",
               "C13.0",
               "C15.0",
               "C17.0",
               "C17.1n7c")
sort(unique(df_full$FA))[sort(unique(df_full$FA)) %in% OBCFA_list]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# 2: Hemolymph vs. Hepatopancreas-----

# > 2.1: Lab Validation -----
# > > 2.1.1: Quantifiable Fatty Acids -----
length(unique(df_LV[df_LV$ug0 > 0,]$FA))
df_LV_quantFA = df_LV |>
  filter(ug0 > 0) |>
  group_by(Tissue, Type, Tag_no_Crab_ID, Time_storage_mo) |>
  summarize(count = length(FA))

banana1 = 
  df_LV_quantFA |>
  filter(Time_storage_mo == 0) |>
  ungroup() |>
  select(!c(Time_storage_mo, Type))
names(banana1)[names(banana1) == "count"] = "baseline"

df_LV_quantFA = left_join(df_LV_quantFA,
                          banana1,
                          by = join_by(Tag_no_Crab_ID, Tissue))
rm(banana1)

df_LV_quantFA$prop = df_LV_quantFA$count / df_LV_quantFA$baseline

ggplot(aes(x = Time_storage_mo,
           y = count,
           color = Type,
           group = interaction(Type, Tissue, Time_storage_mo)),
       data = df_LV_quantFA) +
  geom_boxplot() +
  scale_color_manual(values = c("firebrick",
                                "dodgerblue"),
                     labels = c("DTS",
                                "FRZ")) +
  labs(x = "Storage time, mo",
       y = "# quantifiable fatty acids") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "grey30", fill = NA)) +
  facet_grid(.~ Tissue,
             labeller = as_labeller(c(HL = "hemolymph",
                                      HP = "hepatopancreas")))
ggsave("CrabDTS_Fig01.tiff",
       device = "tiff",
       height = 3,
       width = 5,
       dpi = res)

glm_LV_quantFA = 
  glmmTMB(prop ~ Type + Tissue * Time_storage_mo + (1|Tag_no_Crab_ID),
          family = "gaussian",
          data = df_LV_quantFA[df_LV_quantFA$Time_storage_mo != 0,])
summary(glm_LV_quantFA)
simulateResiduals(glm_LV_quantFA,
                  plot = TRUE)

# > > 2.1.2: Hemolymph lipid content -----
df_LV_HL_quant = 
  df_LV_HL |>
  select(Tag_no_Crab_ID,
         Time_storage_mo,
         sample_vol_uL,
         Type,
         TFA_ug) |>
  unique()
df_LV_HL_quant$lipid_ugpuL = 
  df_LV_HL_quant$TFA_ug /
  df_LV_HL_quant$sample_vol_uL 

plot_LV_HL_quant =
  ggplot(aes(x = Time_storage_mo,
             y = lipid_ugpuL,
             color = Type,
             group = interaction(Type, Time_storage_mo)),
         data = df_LV_HL_quant) +
  geom_boxplot() +
  # geom_text(aes(label = Tag_no_Crab_ID),
  #           size = 3) +
  scale_color_manual(values = c("firebrick",
                                "dodgerblue"),
                     labels = c("DTS",
                                "FRZ")) +
  labs(x = "Storage time, mo",
       y = "fatty acid content, μg/μL") +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification.inside = c(1, 1))
plot_LV_HL_quant

ggsave("CrabDTS_Fig02.tiff",
       device = "tiff",
       height = 3,
       width = 4,
       dpi = res)

boxplot(lipid_ugpuL ~ Type,
        data = df_LV_HL_quant[df_LV_HL_quant$Time_storage_mo == 5,])
wilcox.test(df_LV_HL_quant[df_LV_HL_quant$Time_storage_mo == 5 &
                             df_LV_HL_quant$Type == "DTS",]$lipid_ugpuL,
            df_LV_HL_quant[df_LV_HL_quant$Time_storage_mo == 5 &
                             df_LV_HL_quant$Type == "Tissue",]$lipid_ugpuL,
            paired = TRUE)


# > 2.1: Field Validation -----
# > > 2.2.1: Quantifiable Fatty Acids -----
length(unique(df_FV[df_FV$ug0 > 0,]$FA))
unique(df_FV[df_FV$ug0 > 0,]$FA)

df_FV_quantFA = 
  df_FV |>
  filter(Type == "DTS" &
           (Dev != "mat" |
              Sex != "M") &
           ug0 > 0) |>
  group_by(Tissue, Tag_no_Crab_ID) |>
  summarize(count = length(FA))

df_FV_quantFA |>
  group_by(Tissue) |>
  summarize(M = mean(count),
            SD = sd(count))

boxplot(count ~ Tissue, data = df_FV_quantFA)
wilcox.test(df_FV_quantFA[df_FV_quantFA$Tissue == "HP",]$count,
            df_FV_quantFA[df_FV_quantFA$Tissue == "HL",]$count,
            paired = TRUE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# 3: Logratio Approach -----
# > 3.1: Lab Validation Data -----
# > > 3.1.1: Organize Data -----

# Combine Odd and Branched Fatty Acids as Bacterial Markers
# Note: no branched FA, but lumping all the same
# Sum OBCFAs
banana = df_LV_HP |>
  filter(FA %in% OBCFA_list) |>
  group_by(SampleID,
           Set,
           Species,
           Date_sampling,
           Date_extraction,
           Type,
           Tissue,
           CW_mm,
           Tag_no_Crab_ID,
           Sex,
           Shell,
           Time_storage_mo,
           TFA_ug) |>
  summarize(perc = sum(perc),
            FA = "OBCFA")

# Attach summed OBCFAs and remove individual FAs
df_LV_HP = 
  bind_rows(df_LV_HP,
            banana) |>
  filter(! FA %in% OBCFA_list)

# Check
sort(unique(df_LV_HP$FA))

# Summarize
df_LV_HP_sum = df_LV_HP |>
  group_by(Type, Time_storage_mo, FA) |>
  summarise(M = mean(perc),
            SD = sd(perc))

# Get FAs that are at least 0.5% of any type/time combination
df_LV_HP_list_0.5 = unique(df_LV_HP_sum[df_LV_HP_sum$M >= 0.5,]$FA)
df_LV_HP_list_0.5
# Get FAs that are at least 5% of any type/time combination
df_LV_HP_list_5 = unique(df_LV_HP_sum[df_LV_HP_sum$M >= 5,]$FA)
df_LV_HP_list_5

# Wide format data frame with reduced set
df_LV_HPw = df_LV_HP[df_LV_HP$FA %in% df_LV_HP_list_0.5, ] |>
  pivot_wider(names_from = FA,
              values_from = perc,
              id_cols = c(Tag_no_Crab_ID, Type, Time_storage_mo))

ma_LV_HPw = df_LV_HPw[ , 4:length(df_LV_HPw) ]

# Replace 0s with minimum values
ma_LV_HPw_0NA = ma_LV_HPw
ma_LV_HPw_0NA[ma_LV_HPw_0NA == 0] = NA
FAmin = apply(ma_LV_HPw_0NA, 2, min, na.rm = T)

for(j in 1:ncol(ma_LV_HPw)) {
  for(i in 1:nrow(ma_LV_HPw)) {
    if(ma_LV_HPw[i, j] == 0) ma_LV_HPw[i, j] <- 0.5 * FAmin[j]
  }
}

# renormalize to 100% ("re-close")
ma_LV_HPw = ma_LV_HPw / rowSums(ma_LV_HPw) * 100

#  G&G2020's output function (verbatim):
print.ratios <- function(rationames, R2, procr=NA, N=10) {
  # function prints ratios and the corresponding R2, optionally Procrustes correlations
  # prints first 10 ratios by default  
  # split names of parts in ratio into FA1 and FA2
  # notice that the ratios can be reported as FA1/FA2 or FA2/FA1  
  foo    <- strsplit(rationames, "/")
  parts  <- matrix(unlist(foo), ncol=2, byrow=TRUE)
  df   <- as.data.frame(parts)[1:N,]
  if(is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2))
    colnames(df) <- c("FA1", "FA2","R2")
  }
  if(!is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2), round(procr[1:N], 3))
    colnames(df) <- c("FA1", "FA2", "R2","Procrustes")
  }
  print(df[1:N,])
}

# > > 3.1.2: Stepwise LR selection -----
# > > > Step 1 ----
LV_HP_step1 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10)
print.ratios(LV_HP_step1$names.top, LV_HP_step1$R2.top)
# FA1     FA2    R2
# 1  C20.2n6c C22.6n3 39.14
# 2   C22.6n3   OBCFA 38.41 <- DHA is LCPUFA, OBCFA for bacteria
# 3   C22.5n3 C22.6n3 37.25
# 4   C20.5n3 C22.6n3 36.77
# 5  C22.4n6c C22.6n3 36.07
# 6  C20.1n9c C22.6n3 35.67
# 7  C20.4n6c C22.6n3 34.70
# 8  C18.2n6c C22.6n3 33.85
# 9  C18.3n3c C22.6n3 32.27
# 10 C18.1n9t C22.6n3 31.99
LV_HP_LR1 = LV_HP_step1$logratios.top[,2]
LV_HP_NR1 = LV_HP_step1$ratios.top[2,]
LV_HP_R21 = LV_HP_step1$R2.top[2]

# > > > Step 2 -----
LV_HP_step2 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR1)
print.ratios(LV_HP_step2$names.top, LV_HP_step2$R2.top)
# FA1      FA2    R2
# 1  C20.1n9c C20.4n6c 56.78 <- 20:1n9 for copes/zoops, ARA is LCPUFA
# 2   C18.4n3 C20.4n6c 56.45
# 3  C16.1n7c  C20.5n3 55.47
# 4  C16.1n7c C20.4n6c 55.32
# 5     C14.0 C20.4n6c 55.20
# 6   C18.4n3  C20.5n3 55.09
# 7  C16.1n7c    OBCFA 54.92
# 8  C16.1n7c  C22.6n3 54.92
# 9   C18.4n3    OBCFA 54.68
# 10  C18.4n3  C22.6n3 54.68
LV_HP_LR2 = cbind(LV_HP_LR1, LV_HP_step2$logratios.top[,1])
LV_HP_NR2 = rbind(LV_HP_NR1, LV_HP_step2$ratios.top[1,])
LV_HP_R22 = c(LV_HP_R21, LV_HP_step2$R2.top[1])

# > > > Step 3 -----
LV_HP_step3 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR2)
print.ratios(LV_HP_step3$names.top, LV_HP_step3$R2.top)
# FA1      FA2    R2
# 1  C16.0  C20.5n3 68.51 <- palmitic common SAFA, EPA is and EFA/diatom marker
# 2  C16.0  C22.6n3 68.20
# 3  C16.0    OBCFA 68.20
# 4  C16.0  C22.5n3 67.91
# 5  C16.0 C20.2n6c 67.73
# 6  C16.0 C20.4n6c 67.53
# 7  C16.0 C20.1n9c 67.53
# 8  C16.0 C22.4n6c 66.78
# 9  C16.0 C18.1n9c 66.59
# 10 C16.0 C18.3n3c 66.14
LV_HP_LR3 = cbind(LV_HP_LR2, LV_HP_step3$logratios.top[,1])
LV_HP_NR3 = rbind(LV_HP_NR2, LV_HP_step3$ratios.top[1,])
LV_HP_R23 = c(LV_HP_R22, LV_HP_step3$R2.top[1])

# > > > Step 4 -----
LV_HP_step4 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR3)
print.ratios(LV_HP_step4$names.top, LV_HP_step4$R2.top)
# FA1      FA2    R2
# 1  C20.4n6c  C20.5n3 75.39 <- both eicosanoid LCPUFAs, browns and reds, resp.
# 2  C20.1n9c  C20.5n3 75.39
# 3     C16.0 C20.1n9c 75.39
# 4     C16.0 C20.4n6c 75.39
# 5     C16.0    C18.0 74.87
# 6     C18.0  C20.5n3 74.87
# 7     C16.0 C18.1n9c 74.71
# 8  C18.1n9c  C20.5n3 74.71
# 9     C18.0  C18.4n3 74.45
# 10    C16.0 C18.2n6c 74.33
LV_HP_LR4 = cbind(LV_HP_LR3, LV_HP_step4$logratios.top[,1])
LV_HP_NR4 = rbind(LV_HP_NR3, LV_HP_step4$ratios.top[1,])
LV_HP_R24 = c(LV_HP_R23, LV_HP_step4$R2.top[1])

# > > > Step 5 -----
LV_HP_step5 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR4)
print.ratios(LV_HP_step5$names.top, LV_HP_step5$R2.top)
# FA1      FA2    R2
# 1     C18.0  C18.4n3 81.32 <- stearic common SAFA, stearidonic for browns
# 2     C18.0 C20.4n6c 81.23
# 3     C18.0 C20.1n9c 81.23
# 4     C16.0    C18.0 81.23
# 5     C18.0  C20.5n3 81.23
# 6     C14.0    C18.0 81.13
# 7     C18.0 C18.3n3c 80.99
# 8  C16.1n7c    C18.0 80.82
# 9     C18.0 C18.2n6c 80.81
# 10    C14.0 C18.1n9c 80.61
LV_HP_LR5 = cbind(LV_HP_LR4, LV_HP_step5$logratios.top[,1])
LV_HP_NR5 = rbind(LV_HP_NR4, LV_HP_step5$ratios.top[1,])
LV_HP_R25 = c(LV_HP_R24, LV_HP_step5$R2.top[1])

# > > > Step 6 -----
LV_HP_step6 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR5)
print.ratios(LV_HP_step6$names.top, LV_HP_step6$R2.top)
# FA1      FA2    R2
# 1  C18.1n9c C22.1n9c 85.68 <- oleic for carnivory, detritus, 22:1 for zoops
# 2  C18.1n9c C20.1n9c 85.53
# 3  C18.1n9c  C20.5n3 85.53
# 4  C18.1n9c C20.4n6c 85.53
# 5     C16.0 C18.1n9c 85.53
# 6  C20.4n6c C22.1n9c 85.41
# 7   C20.5n3 C22.1n9c 85.41
# 8     C16.0 C22.1n9c 85.41
# 9  C20.1n9c C22.1n9c 85.41
# 10 C22.1n9c    OBCFA 85.37
LV_HP_LR6 = cbind(LV_HP_LR5, LV_HP_step6$logratios.top[,1])
LV_HP_NR6 = rbind(LV_HP_NR5, LV_HP_step6$ratios.top[1,])
LV_HP_R26 = c(LV_HP_R25, LV_HP_step6$R2.top[1])

# > > > Step 7 -----
LV_HP_step7 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR6)
print.ratios(LV_HP_step7$names.top, LV_HP_step7$R2.top)
# FA1      FA2    R2
# 1  C20.4n6c C22.1n9c 89.54
# 2  C20.1n9c C22.1n9c 89.54
# 3  C18.1n9c C20.1n9c 89.54
# 4  C18.1n9c C20.4n6c 89.54
# 5   C20.5n3 C22.1n9c 89.54
# 6  C18.1n9c  C20.5n3 89.54 <- oleic more abundant in browns, EPA in reds
# 7     C16.0 C18.1n9c 89.54
# 8     C16.0 C22.1n9c 89.54
# 9   C18.4n3  C20.5n3 89.10
# 10    C18.0 C20.1n9c 89.10
LV_HP_LR7 = cbind(LV_HP_LR6, LV_HP_step7$logratios.top[,6])
LV_HP_NR7 = rbind(LV_HP_NR6, LV_HP_step7$ratios.top[6,])
LV_HP_R27 = c(LV_HP_R26, LV_HP_step7$R2.top[6])

# > > > Step 8 -----
LV_HP_step8 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR7)
print.ratios(LV_HP_step8$names.top, LV_HP_step8$R2.top)
# FA1      FA2    R2
# 1   C18.4n3  C20.5n3 92.44
# 2   C18.4n3 C20.4n6c 92.44 <- both are brown indicators, but ARA is LCPUFA
# 3     C18.0 C20.4n6c 92.44
# 4  C18.1n9c  C18.4n3 92.44
# 5     C18.0 C20.1n9c 92.44
# 6     C18.0 C18.1n9c 92.44
# 7   C18.4n3 C20.1n9c 92.44
# 8   C18.4n3 C22.1n9c 92.44
# 9     C16.0  C18.4n3 92.44
# 10    C16.0    C18.0 92.44
LV_HP_LR8 = cbind(LV_HP_LR7, LV_HP_step8$logratios.top[,2])
LV_HP_NR8 = rbind(LV_HP_NR7, LV_HP_step8$ratios.top[2,])
LV_HP_R28 = c(LV_HP_R27, LV_HP_step8$R2.top[2])

# > > > Step 9
LV_HP_step9 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR8)
print.ratios(LV_HP_step9$names.top, LV_HP_step9$R2.top)
# FA1      FA2    R2
# 1  C22.4n6c  C22.6n3 94.40 <- both adrenic and DHA are LCPUFA, but adrenic not straightforward
# 2  C22.4n6c    OBCFA 94.40
# 3  C20.4n6c C22.4n6c 94.38
# 4     C16.0 C22.4n6c 94.38
# 5  C20.1n9c C22.4n6c 94.38
# 6     C18.0 C22.4n6c 94.38
# 7  C18.1n9c C22.4n6c 94.38
# 8  C22.1n9c C22.4n6c 94.38
# 9   C18.4n3 C22.4n6c 94.38
# 10  C20.5n3 C22.4n6c 94.38
LV_HP_LR9 = cbind(LV_HP_LR8, LV_HP_step9$logratios.top[,1])
LV_HP_NR9 = rbind(LV_HP_NR8, LV_HP_step9$ratios.top[1,])
LV_HP_R29 = c(LV_HP_R28, LV_HP_step9$R2.top[1])

# > > > Step 10 -----
LV_HP_step10 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR9)
print.ratios(LV_HP_step10$names.top, LV_HP_step10$R2.top)
# FA1      FA2    R2
# 1  C16.1n7c  C20.5n3 96.09
# 2  C16.1n7c    C18.0 96.09
# 3  C16.1n7c C20.4n6c 96.09
# 4  C16.1n7c  C18.4n3 96.09
# 5     C16.0 C16.1n7c 96.09 <- diatom marker, sensu Copeman et al.
# 6  C16.1n7c C20.1n9c 96.09
# 7  C16.1n7c C22.1n9c 96.09
# 8  C16.1n7c C18.1n9c 96.09
# 9  C22.1n9c  C22.6n3 96.06
# 10    C18.0  C22.6n3 96.06
LV_HP_LR10 = cbind(LV_HP_LR9, LV_HP_step10$logratios.top[,5])
LV_HP_NR10 = rbind(LV_HP_NR9, LV_HP_step10$ratios.top[5,])
LV_HP_R210 = c(LV_HP_R29, LV_HP_step10$R2.top[5])

# > > > Step 11 -----
LV_HP_step11 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR10)
print.ratios(LV_HP_step11$names.top, LV_HP_step11$R2.top)
# FA1      FA2    R2
# 1  C22.1n9c  C22.6n3 97.75
# 2  C20.4n6c  C22.6n3 97.75 <- both arguable LCEFA, ARA common in browns
# 3   C18.4n3 C22.4n6c 97.75
# 4  C22.1n9c    OBCFA 97.75
# 5  C20.1n9c    OBCFA 97.75
# 6  C22.1n9c C22.4n6c 97.75
# 7  C18.1n9c  C22.6n3 97.75
# 8   C20.5n3    OBCFA 97.75
# 9  C18.1n9c C22.4n6c 97.75
# 10    C16.0    OBCFA 97.75
LV_HP_LR11 = cbind(LV_HP_LR10, LV_HP_step11$logratios.top[,2])
LV_HP_NR11 = rbind(LV_HP_NR10, LV_HP_step11$ratios.top[2,])
LV_HP_R211 = c(LV_HP_R210, LV_HP_step11$R2.top[2])

# > > > Step 12 -----
LV_HP_step12 = STEP(data = ma_LV_HPw, nsteps = 1, top = 10, previous = LV_HP_LR11)
print.ratios(LV_HP_step12$names.top, LV_HP_step12$R2.top)
# FA1      FA2    R2
# 1     C14.0  C22.5n3 98.55
# 2     C14.0 C20.2n6c 98.53
# 3   C20.5n3  C22.5n3 98.48
# 4  C18.1n9c  C22.5n3 98.48
# 5     C16.0  C22.5n3 98.48
# 6   C18.4n3  C22.5n3 98.48
# 7  C22.4n6c  C22.5n3 98.48
# 8   C22.5n3    OBCFA 98.48
# 9  C20.4n6c  C22.5n3 98.48
# 10  C22.5n3  C22.6n3 98.48

# Stop here. DPA is only somwhat informative and diminishing returns

# > > 3.1.3: Organize step 11 LRs -----
rownames(LV_HP_NR11) = paste("Step", 1:11, sep="")
colnames(LV_HP_NR11) = c("FA1","FA2")
LV_HP_FR = as.data.frame(cbind(LV_HP_NR11,
                               Ratio = paste(colnames(ma_LV_HPw)[LV_HP_NR11[,1]],
                                             "/",
                                             colnames(ma_LV_HPw)[LV_HP_NR11[,2]],
                                             sep="")))
LV_HP_FR$FA1_lab = colnames(ma_LV_HPw)[LV_HP_NR11[,1]]
LV_HP_FR$FA2_lab = colnames(ma_LV_HPw)[LV_HP_NR11[,2]]
LV_HP_FR$step = 1:11

LV_HP_FR$R2_cum = LV_HP_R211
LV_HP_FR = LV_HP_FR |>
  arrange(step) |>
  mutate(R2 = R2_cum - lag(R2_cum,
                           default = 0))

LV_HP_FR

LV_HP_FR_export = LV_HP_FR
LV_HP_FR_export$R2_cum = round(LV_HP_FR_export$R2_cum, 3)
LV_HP_FR_export$R2 = round(LV_HP_FR_export$R2, 3)
LV_HP_FR_export
write.csv(LV_HP_FR_export, "CrabDTS_Table02_partial.csv")

colnames(LV_HP_LR11) = LV_HP_FR[,3]

# LV_HP_PiR = sort(unique(as.numeric(LV_HP_NR11))) # "parts in ratios"
LV_HP_PiR = unique(as.numeric(t(LV_HP_NR11))) # "parts in ratios"
colnames(ma_LV_HPw)[LV_HP_PiR]

# > > 3.1.4: Prep acyclic graph -----
LV_HP_PiR_dim = data.frame(FA = colnames(ma_LV_HPw)[LV_HP_PiR],
                           dim1 = c(2, # C22.6n3
                                    1, # OBCFA
                                    1, # C20.1n9c
                                    2, # C20.4n6c
                                    1, # C16.0
                                    2, # C20.5n3
                                    4, # C18.0
                                    3, # C18.4n3
                                    3, # C18.1n9c
                                    4, # C22.1n9c
                                    3, # C22.4n6c
                                    1),# C16.1n7c
                           dim2 = c(4, # C22.6n3
                                    4, # OBCFA
                                    3, # C20.1n9c
                                    3, # C20.4n6c
                                    2, # C16.0
                                    2, # C20.5n3
                                    3, # C18.0
                                    3, # C18.4n3
                                    2, # C18.1n9c
                                    2, # C22.1n9c
                                    4, # C22.4n6c
                                    1))# C16.1n7c
LV_HP_PiR_dim$FA_labs = str_replace_all(LV_HP_PiR_dim$FA,
                                        FA_labeller)

LV_HP_FR = left_join(LV_HP_FR,
                     LV_HP_PiR_dim,
                     by = c("FA1_lab" = "FA"))
LV_HP_FR = left_join(LV_HP_FR,
                     LV_HP_PiR_dim,
                     by = c("FA2_lab" = "FA"))

# > > 3.1.5: Prep PCA ----
LV_HP_PCA = PCA(LV_HP_LR11, weight = FALSE)

PLOT.PCA(LV_HP_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

# Extract and calculate contribution vectors
# This bit of code is derived from the PLOT.PCA function in easyCODA
axes_inv = c(1, 1)
dim = c(1, 2)
LV_HP_RPC = LV_HP_PCA$rowcoord[, dim] %*% diag(LV_HP_PCA$sv[dim] * axes_inv)
LV_HP_CSC = LV_HP_PCA$colcoord[, dim] %*% diag(axes_inv)
LV_HP_CPC = LV_HP_CSC %*% diag(LV_HP_PCA$sv[dim])
LV_HP_CCC = LV_HP_CSC * sqrt(LV_HP_PCA$colmass)
LV_HP_CRD = LV_HP_CCC

# rescale contribution vectors for visibility
vectorscale = 1.25
LV_HP_CRD_scale = LV_HP_CRD * vectorscale

# Extract points
df_LV_HP_RPC = cbind(df_LV_HPw[, 1:3], LV_HP_RPC)
colnames(df_LV_HP_RPC) = c(colnames(df_LV_HPw[, 1:3]),'dim1','dim2')

# Aggregate by type and storage time
df_LV_HP_RPC_MSD = df_LV_HP_RPC |>
  group_by(Type, Time_storage_mo) |>
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SD = sd(dim1),
            dim2SD = sd(dim2))

# Calculate variance explained by PCA dimensions
LV_HP_PCA_perc_1 = 100 * LV_HP_PCA$sv[dim[1]]^2 / sum(LV_HP_PCA$sv^2)
LV_HP_PCA_perc_2 = 100 * LV_HP_PCA$sv[dim[2]]^2 / sum(LV_HP_PCA$sv^2)

df_LV_HP_RPC$ID2 = paste(df_LV_HP_RPC$Type,
                         df_LV_HP_RPC$Time_storage_mo,
                         df_LV_HP_RPC$Tag_no_Crab_ID,
                         sep = ".")

# Format arrows
df_LV_HP_arrows = data.frame(dim1 = LV_HP_CRD_scale[ ,1],
                             dim2 = LV_HP_CRD_scale[, 2],
                             LR = LV_HP_PCA$colnames)

df_LV_HP_arrows$step = 1:11

df_LV_HP_arrows$angle = 90 - ((atan2(df_LV_HP_arrows$dim1,
                                     df_LV_HP_arrows$dim2) *
                                 180) /
                                pi)
df_LV_HP_arrows$angle = ifelse(df_LV_HP_arrows$angle > 90 & df_LV_HP_arrows$angle < 270,
                               df_LV_HP_arrows$angle + 180,
                               df_LV_HP_arrows$angle)
df_LV_HP_arrows$hjust = ifelse(df_LV_HP_arrows$angle > 90,
                               1,
                               0)

# Format arrow labels
df_LV_HP_arrows$labs = df_LV_HP_arrows$LR
df_LV_HP_arrows$labs = str_replace_all(df_LV_HP_arrows$labs,
                                       FA_labeller)
df_LV_HP_arrows$labs = str_replace_all(df_LV_HP_arrows$labs,
                                       "/",
                                       " / ")

# Create hulls
df_LV_HP_RPC_hull = df_LV_HP_RPC |>
  group_by(Type, Time_storage_mo) |>
  slice(chull(dim1, dim2))

df_LV_HP_RPC_hull2 = df_LV_HP_RPC |>
  group_by(Tag_no_Crab_ID) |>
  slice(chull(dim1, dim2))

# > > 3.1.6: Plot PCA -----
# > > > 3.1.6.1: By Type, Storage Time -----
plot_LV_HP_PCA = ggplot() +
  # geom_polygon(aes(x = dim1,
  #                  y = dim2,
  #                  color = Time_storage_mo,
  #                  linetype = Type,
  #                  group = interaction(Time_storage_mo, Type)),
  #              linewidth = 1,
  #              fill = NA,
  #              data = df_LV_HP_RPC_hull) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = Time_storage_mo,
                 shape = Type),
             alpha = 0.5,
             size = 2,
             data = df_LV_HP_RPC) +
  # geom_point(aes(x = ifelse(dim1 > 0,
  #                           dim1 + 0.075 * cos(angle * pi/180),
  #                           dim1 - 0.075 * cos(angle * pi/180)),
  #                y = ifelse(dim1 > 0,
  #                           dim2 + 0.075 * sin(angle * pi/180),
  #                           dim2 - 0.075 * sin(angle * pi/180))),
  #            colour = "grey50",
  #            size = 4,
  #            data = df_LV_HP_arrows)+
  # geom_text(aes(x = ifelse(dim1 > 0,
  #                          dim1 + 0.075 * cos(angle * pi/180),
#                          dim1 - 0.075 * cos(angle * pi/180)),
#               y = ifelse(dim1 > 0,
#                          dim2 + 0.075 * sin(angle * pi/180),
#                          dim2 - 0.075 * sin(angle * pi/180)),
#               label = step),
#           colour = "white",
#           size = unit(2, "pt"),
#           data = df_LV_HP_arrows)+
# geom_segment(aes(x = 0,
#                  y = 0,
#                  xend = LV_HP_CRD_scale[ ,1],
#                  yend = LV_HP_CRD_scale[ ,2]),
#              color = "grey50",
#              arrow = arrow(length = unit(0.2, "cm"))) +
geom_errorbar(aes(x = dim1M,
                  ymin = dim2M - dim2SD,
                  ymax = dim2M + dim2SD,
                  color = Time_storage_mo),
              width = 0,
              linewidth = 0.5,
              data = df_LV_HP_RPC_MSD) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - dim1SD,
                     xmax = dim1M + dim1SD,
                     color = Time_storage_mo),
                 height = 0,
                 linewidth = 0.5,
                 data = df_LV_HP_RPC_MSD) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = Time_storage_mo,
                 shape = Type),
             size = 3,
             data = df_LV_HP_RPC_MSD) +
  scale_color_viridis_c(option = "magma",
                        begin = 0.2,
                        end = 0.9,
                        direction = -1) +
  # scale_linetype_manual(values = c(1, 2),
  #                       labels = c("DTS", "FRZ"),
  #                       name = "Type") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = NULL,
       # x = paste("PC 1 (",
       #           round(LV_HP_PCA_perc_1, 1),
       #           "%)",
       #           sep=""),
       y = paste("PC 2 (",
                 round(LV_HP_PCA_perc_2, 1),
                 "%)",
                 sep=""),
       color = "Storage \nTime") +
  coord_fixed(xlim = c(-1, 0.9),
              ylim = c(-0.9, 0.9)) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        legend.box = "vertical",
        axis.text.y.right = element_blank(),
        # axis.text.y.right = element_text(colour="grey50",
        #                                  size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))
plot_LV_HP_PCA

# > > > 3.1.6.2: By Individual -----
plot_LV_HP_PCA2 = ggplot() +
  geom_polygon(aes(x = dim1,
                   y = dim2,
                   color = Tag_no_Crab_ID),
               linewidth = 1,
               fill = NA,
               data = df_LV_HP_RPC_hull2) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = Tag_no_Crab_ID,
                 shape = Type),
             alpha = 0.5,
             size = 2,
             data = df_LV_HP_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = Time_storage_mo,
  #               label = ID2),
  #           size = 2,
  #           data = df_LV_HP_RPC) +
  # geom_point(aes(x = ifelse(dim1 > 0,
  #                           dim1 + 0.075 * cos(angle * pi/180),
  #                           dim1 - 0.075 * cos(angle * pi/180)),
  #                y = ifelse(dim1 > 0,
  #                           dim2 + 0.075 * sin(angle * pi/180),
#                           dim2 - 0.075 * sin(angle * pi/180))),
#            colour = "grey50",
#            size = 4,
#            data = df_LV_HP_arrows)+
# geom_text(aes(x = ifelse(dim1 > 0,
#                          dim1 + 0.075 * cos(angle * pi/180),
#                          dim1 - 0.075 * cos(angle * pi/180)),
#               y = ifelse(dim1 > 0,
#                          dim2 + 0.075 * sin(angle * pi/180),
#                          dim2 - 0.075 * sin(angle * pi/180)),
#               label = step),
#           colour = "white",
#           size = unit(2, "pt"),
#           data = df_LV_HP_arrows)+
# geom_segment(aes(x = 0,
#                  y = 0,
#                  xend = LV_HP_CRD_scale[ ,1],
#                  yend = LV_HP_CRD_scale[ ,2]),
#              color = "grey50",
#              arrow = arrow(length = unit(0.2, "cm"))) +
# geom_text(aes(x = dim1 * 1.1,
#               y = dim2 * 1.1,
#               label = labs),
#           colour = "grey50",
#           size = 2,
#           hjust = df_LV_HP_arrows$hjust,
#           angle = df_LV_HP_arrows$angle,
#           data = df_LV_HP_arrows)+
# geom_errorbar(aes(x = dim1M,
#                   ymin = dim2M - dim2SD,
#                   ymax = dim2M + dim2SD),
#               color = "grey70",
#               width = 0,
#               linewidth = 1,
#               data = df_LV_HP_RPC_MSD) +
# geom_errorbarh(aes(y = dim2M,
#                    xmin = dim1M - dim1SD,
#                    xmax = dim1M + dim1SD),
#                color = "grey70",
#                height = 0,
#                linewidth = 1,
#                data = df_LV_HP_RPC_MSD) +
# geom_point(aes(x = dim1M,
#                y = dim2M,
#                shape = Type),
#            color = "grey70",
#            size = 3,
#            data = df_LV_HP_RPC_MSD) +
scale_color_viridis_d(option = "mako",
                      begin = 0.2,
                      end = 0.9,
                      direction = -1,
                      guide = "none") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     guide = "none") +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(LV_HP_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       y = NULL,
       # y = paste("PC 2 (",
       #           round(LV_HP_PCA_perc_2, 1),
       #           "%)",
       #           sep=""),
       color = "Crab ID",
       shape = "Type") +
  coord_fixed(xlim = c(-1, 0.9),
              ylim = c(-0.9, 0.9)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        axis.text.y = element_blank(),
        legend.box = "vertical",
        axis.text.y.right = element_blank(),
        # axis.text.y.right = element_text(colour="grey50",
        #                                  size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))
plot_LV_HP_PCA2

# > > > 3.1.6.3: Baseplot with Vectors -----
plot_LV_HP_PCA_base =
  ggplot() +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = LV_HP_CRD_scale[ ,1],
                   yend = LV_HP_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(aes(x = ifelse(dim1 > 0,
                            dim1 + 0.085 * cos(angle * pi/180),
                            dim1 - 0.085 * cos(angle * pi/180)),
                 y = ifelse(dim1 > 0,
                            dim2 + 0.085 * sin(angle * pi/180),
                            dim2 - 0.085 * sin(angle * pi/180))),
             colour = "grey50",
             size = 3,
             data = df_LV_HP_arrows)+
  geom_text(aes(x = ifelse(dim1 > 0,
                           dim1 + 0.085 * cos(angle * pi/180),
                           dim1 - 0.085 * cos(angle * pi/180)),
                y = ifelse(dim1 > 0,
                           dim2 + 0.085 * sin(angle * pi/180),
                           dim2 - 0.085 * sin(angle * pi/180)),
                label = step),
            colour = "white",
            size = unit(2, "pt"),
            data = df_LV_HP_arrows)+
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = NULL,
       # x = paste("PC 1 (",
       #           round(LV_HP_PCA_perc_1, 1),
       #           "%)",
       #           sep=""),
       # y = paste("PC 2 (",
       #           round(LV_HP_PCA_perc_2, 1),
       #           "%)",
       #           sep=""),
       y = NULL) +
  coord_fixed(xlim = c(-1, 0.9),
              ylim = c(-0.9, 0.9)) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        axis.text.y = element_blank(),
        legend.box = "vertical",
        axis.text.y.right = element_text(colour="grey50",
                                         size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))
plot_LV_HP_PCA_base

# > > 3.1.7: Plot acyclic graph -----
plot_LV_HP_LR_graph = ggplot() +
  geom_segment(aes(x = dim1.x,
                   y = dim2.x,
                   xend = dim1.y,
                   yend = dim2.y,
                   color = step),
               linewidth = 1,
               data = LV_HP_FR) +
  geom_label(aes(x = dim1,
                 y = dim2,
                 label = FA_labs),
             size = unit(3, "pt"),
             data = LV_HP_PiR_dim) +
  geom_point(aes(x = (dim1.x + dim1.y)/2,
                 y = (dim2.x + dim2.y)/2,
                 color = step),
             size = 6,
             data = LV_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2,
                y = (dim2.x + dim2.y)/2,
                label = step),
            color = "white",
            size = unit(3, "pt"),
            data = LV_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2 + 0.125,
                y = (dim2.x + dim2.y)/2 - 0.125,
                label = round(R2*100, 0),
                color = step),
            size = unit(2, "pt"),
            data = LV_HP_FR) +
  scale_color_gradient(high = "dodgerblue4",
                       low = "dodgerblue1") +
  coord_cartesian(xlim = c(0.5, 4.5),
                  clip = "off") +
  theme_void() +
  theme(legend.position = "none")
plot_LV_HP_LR_graph 

plot_LV_HP_LR_graph/ (plot_LV_HP_PCA | plot_LV_HP_PCA2) +
  plot_annotation(tag_levels = "A") 

# > > 3.1.8: Export acylic graph and PCAs ----- 
# tiff("plot_C_maenas_LV_HP_LR_doc.tiff",
#      height = res * 6.5,
#      width = res * 7.5,
#      res = res)
# plot_LV_HP_LR_graph / (plot_LV_HP_PCA | plot_LV_HP_PCA2) +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(heights = c(1, 2))
# dev.off()
# 
# tiff("plot_C_maenas_LV_HP_LR_wide.tiff",
#      height = res * 5,
#      width = res * 12,
#      res = res)
# plot_LV_HP_LR_graph + plot_LV_HP_PCA + plot_LV_HP_PCA2 +
#   plot_annotation(tag_levels = "A")
# dev.off()

# > > 3.1.9: Univariate LR plot -----
df_LV_HP_LR11w = df_LV_HPw[,1:3]
df_LV_HP_LR11w$TypeMo = paste(df_LV_HP_LR11w$Type,
                              df_LV_HP_LR11w$Time_storage_mo,
                              sep = "_")

df_LV_HP_LR11w$TypeMo = factor(df_LV_HP_LR11w$TypeMo,
                               levels = c("Tissue_0",
                                          "DTS_2",
                                          "DTS_3",
                                          "DTS_5",
                                          "Tissue_5"))

df_LV_HP_LR11w = cbind(df_LV_HP_LR11w,
                       LV_HP_LR11)

df_LV_HP_LR11 = df_LV_HP_LR11w |>
  pivot_longer(cols = 5:length(df_LV_HP_LR11w),
               names_to = "FAs",
               values_to = "LR")

df_LV_HP_LR11$FAs = str_replace_all(df_LV_HP_LR11$FAs,
                                    FA_labeller)

df_LV_HP_LR11$FAs = factor(df_LV_HP_LR11$FAs,
                           levels = str_replace_all(LV_HP_FR$Ratio,
                                                    FA_labeller))
df_LV_HP_LR11_sum = df_LV_HP_LR11 |>
  group_by(TypeMo, FAs) |>
  summarise(M = mean(LR),
            SD = sd(LR))

# > > > 3.1.9.1: Plot -----
plot_LV_HP_LR_cat = ggplot() +
  geom_hline(aes(yintercept = y),
             linetype = 2,
             color = "grey75",
             data = labels_logratios) +
  geom_point(aes(x = TypeMo,
                 y = LR,
                 color = Tag_no_Crab_ID),
             data = df_LV_HP_LR11) +
  geom_errorbar(aes(x = TypeMo,
                    ymin = M - SD,
                    ymax = M + SD),
                width = 0,
                color = "firebrick",
                data = df_LV_HP_LR11_sum) +
  geom_point(aes(x = TypeMo,
                 y = M),
             size = 3,
             shape = 9,
             color = "firebrick",
             data = df_LV_HP_LR11_sum) +
  scale_color_viridis_d(option = "mako",
                        begin = 0.2,
                        end = 0.9) +
  scale_x_discrete(labels = c("FRZ, 0 mo",
                              "DTS, 2 mo",
                              "DTS, 3 mo",
                              "DTS, 5 mo",
                              "FRZ, 5 mo")) +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),
                                         breaks = exp(labels_logratios$y),
                                         labels = labels_logratios$label,
                                         name = "Ratio"),
                     limits = c(-5.5, 4)) +
  labs(x = "Sample type, storage time",
       y = "Logratio",
       color = "Crab ID") +
  theme_classic() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text = element_text(size = unit(4, "pt")),
        axis.title = element_text(size = unit(6, "pt")),
        axis.text.y.right = element_text(color="grey70"),
        axis.ticks.y.right = element_line(color="grey70"),
        axis.line.y.right = element_line(color="grey70"),
        axis.title.y.right = element_text(color = "grey70"),
        legend.background = element_blank(),
        legend.text = element_text(size = unit(4, "pt")),
        legend.title = element_text(size = unit(6, "pt")),
        legend.spacing = unit(0, "pt"),
        strip.text = element_text(size = unit(4, "pt"))) +
  facet_grid(.~ FAs)
plot_LV_HP_LR_cat

# > > > 3.1.9.2: Export -----
# tiff("plot_C_maenas_LV_HP_LR_cat.tiff",
#      height = res * 4,
#      width = res * 7.5,
#      res = res)
# plot_LV_HP_LR_cat
# dev.off()

# > > 3.1.10: PERMANOVA on LRs -----
# Sample type and storage time are severely unbalanced; this strategy was 
# chosen for practical sample numbers and costs. To analyze the data, we 
# separate the data into two sets: one to assess Sample type at month 5 and 
# the other to assess Time for DTS samples. 

# Isolate month 5
df_LV_HP_LR11w_Type = 
  df_LV_HP_LR11w |>
  filter(Time_storage_mo == 5)

# Isolate month 3 DTS and month 0 baseline
df_LV_HP_LR11w_Type2 = 
  df_LV_HP_LR11w |>
  filter(Time_storage_mo == 2 |
           (Time_storage_mo == 5 &
              Type == "Tissue"))

# Isolate DTS samples, and month 0 frozen as baseline
df_LV_HP_LR11w_Time = 
  df_LV_HP_LR11w |>
  filter(Time_storage_mo == 0 |
           Type == "DTS")

# Isolate frozen samples
df_LV_HP_LR11w_Time2 = 
  df_LV_HP_LR11w |>
  filter(Type == "Tissue")

# To control for the repeated-measures nature of the experiment, we constrain 
# permutations based on the individual. 
LV_LR11_Type_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HP_LR11w_Type$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_LR11_Type2_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HP_LR11w_Type2$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_LR11_Time_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HP_LR11w_Time$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_LR11_Time2_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HP_LR11w_Time2$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

# We run them once to check our general structure. We do not use this output for
# inference because we know the pseudo-Fs are not calculated how we want.
LV_LR11_Type_PERMANOVA_base = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 5 ,] ~ 
                                        Tag_no_Crab_ID +
                                        Type,
                                      data = df_LV_HP_LR11w_Type,
                                      method = "euclidean",
                                      permutations = LV_LR11_Type_perm, 
                                      by = "terms")

LV_LR11_Type_PERMANOVA_base

LV_LR11_Type2_PERMANOVA_base = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 2 |
                                                    (df_LV_HP_LR11w$Time_storage_mo == 5 &
                                                       df_LV_HP_LR11w$Type == "Tissue"),] ~ 
                                         Tag_no_Crab_ID +
                                         Type,
                                       data = df_LV_HP_LR11w_Type2,
                                       method = "euclidean",
                                       permutations = LV_LR11_Type2_perm, 
                                       by = "terms")

LV_LR11_Type2_PERMANOVA_base

LV_LR11_Time_PERMANOVA_base = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 0 |
                                                   df_LV_HP_LR11w$Type == "DTS" ,] ~ 
                                        Tag_no_Crab_ID +
                                        Time_storage_mo,
                                      data = df_LV_HP_LR11w_Time,
                                      method = "euclidean",
                                      permutations = LV_LR11_Time_perm, 
                                      by = "terms")

LV_LR11_Time_PERMANOVA_base

LV_LR11_Time2_PERMANOVA_base = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Type == "Tissue" ,] ~ 
                                         Tag_no_Crab_ID +
                                         Time_storage_mo,
                                       data = df_LV_HP_LR11w_Time2,
                                       method = "euclidean",
                                       permutations = LV_LR11_Time2_perm, 
                                       by = "terms")

LV_LR11_Time2_PERMANOVA_base

# Note that for both, we include Tag_no_Crab_ID to allow proper estimation of
# the pseudo-F statistics. These rely on proper denominators relevant to each 
# term. Using the residual as done with adonis2() default is note (always) appropriate 
# (Anderson and Ter Braak 2003 Journal of Statistical Computation and 
# Simulation). Should be based on exchangeable units, which are the individual
# crab. The method used below is from Bakker 2024, Applied Multivariate 
# Statistics in R.

# Check number of possible permutations
numPerms(df_LV_HP_LR11w_Type,
         LV_LR11_Type_perm)

numPerms(df_LV_HP_LR11w_Type2,
         LV_LR11_Type2_perm)

numPerms(df_LV_HP_LR11w_Time,
         LV_LR11_Time_perm)

numPerms(df_LV_HP_LR11w_Time2,
         LV_LR11_Time2_perm)

# Run PERMANOVA

# Type: 

# F_Type = MS_Type / MS_Residual

# Type is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_LR11_Type_perms = rbind(1:nrow(df_LV_HP_LR11w_Type),
                           shuffleSet(n = nrow(df_LV_HP_LR11w_Type),
                                      control = LV_LR11_Type_perm,
                                      nset = 9999))

LV_LR11_Type_PERMANOVA = matrix(nrow = nrow(LV_LR11_Type_perms),
                                ncol = 4)

colnames(LV_LR11_Type_PERMANOVA) = c("Tag_no_Crab_ID",
                                     "Type",
                                     "Residual",
                                     "Total")

for (i in 1:nrow(LV_LR11_Type_perms)) {
  df_temp = df_LV_HP_LR11w_Type[LV_LR11_Type_perms[i, ], ]
  res_temp = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 5 ,] ~
                       Tag_no_Crab_ID +
                       Type,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_LR11_Type_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_LR11_Type_PERMANOVA = 
  LV_LR11_Type_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 9))

# Calculate p-value
LV_LR11_Type_PERMANOVA$F_Type[1]
with(LV_LR11_Type_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

set.seed(NULL)

# Type2: 

# F_Type = MS_Type / MS_Residual

# Type is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_LR11_Type2_perms = rbind(1:nrow(df_LV_HP_LR11w_Type2),
                            shuffleSet(n = nrow(df_LV_HP_LR11w_Type2),
                                       control = LV_LR11_Type2_perm,
                                       nset = 9999))

LV_LR11_Type2_PERMANOVA = matrix(nrow = nrow(LV_LR11_Type2_perms),
                                 ncol = 4)

colnames(LV_LR11_Type2_PERMANOVA) = c("Tag_no_Crab_ID",
                                      "Type",
                                      "Residual",
                                      "Total")

for (i in 1:nrow(LV_LR11_Type2_perms)) {
  df_temp = df_LV_HP_LR11w_Type2[LV_LR11_Type2_perms[i, ], ]
  res_temp = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 2 |
                                  (df_LV_HP_LR11w$Time_storage_mo == 5 &
                                     df_LV_HP_LR11w$Type == "Tissue"),] ~
                       Tag_no_Crab_ID +
                       Type,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_LR11_Type2_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_LR11_Type2_PERMANOVA = 
  LV_LR11_Type2_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 9))

# Calculate p-value
LV_LR11_Type2_PERMANOVA$F_Type[1]
with(LV_LR11_Type2_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

set.seed(NULL)

# Time: 

# F_Time_storage_mo = MS_Time_storage_mo / MS_Residual

# Time is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_LR11_Time_perms = rbind(1:nrow(df_LV_HP_LR11w_Time),
                           shuffleSet(n = nrow(df_LV_HP_LR11w_Time),
                                      control = LV_LR11_Time_perm,
                                      nset = 9999))

LV_LR11_Time_PERMANOVA = matrix(nrow = nrow(LV_LR11_Time_perms),
                                ncol = 4)

colnames(LV_LR11_Time_PERMANOVA) = c("Tag_no_Crab_ID",
                                     "Time_storage_mo",
                                     "Residual",
                                     "Total")

for (i in 1:nrow(LV_LR11_Time_perms)) {
  df_temp = df_LV_HP_LR11w_Time[LV_LR11_Time_perms[i, ], ]
  res_temp = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Time_storage_mo == 0 |
                                  df_LV_HP_LR11w$Type == "DTS" ,] ~
                       Tag_no_Crab_ID +
                       Time_storage_mo,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_LR11_Time_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_LR11_Time_PERMANOVA = 
  LV_LR11_Time_PERMANOVA |>
  data.frame() |>
  mutate(F_Time_storage_mo = (Time_storage_mo / 1) / (Residual / 9))

# Calculate p-value
LV_LR11_Time_PERMANOVA$F_Time_storage_mo[1]
with(LV_LR11_Time_PERMANOVA, 
     sum(F_Time_storage_mo >= F_Time_storage_mo[1]) / 
       length(F_Time_storage_mo))

set.seed(NULL)

# Time2: 

# F_Time_storage_mo = MS_Time_storage_mo / MS_Residual

# Time is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_LR11_Time2_perms = rbind(1:nrow(df_LV_HP_LR11w_Time2),
                            shuffleSet(n = nrow(df_LV_HP_LR11w_Time2),
                                       control = LV_LR11_Time2_perm,
                                       nset = 9999))

LV_LR11_Time2_PERMANOVA = matrix(nrow = nrow(LV_LR11_Time2_perms),
                                 ncol = 4)

colnames(LV_LR11_Time2_PERMANOVA) = c("Tag_no_Crab_ID",
                                      "Time_storage_mo",
                                      "Residual",
                                      "Total")

for (i in 1:nrow(LV_LR11_Time2_perms)) {
  df_temp = df_LV_HP_LR11w_Time2[LV_LR11_Time2_perms[i, ], ]
  res_temp = adonis2(LV_HP_LR11[df_LV_HP_LR11w$Type == "Tissue" ,] ~
                       Tag_no_Crab_ID +
                       Time_storage_mo,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_LR11_Time2_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_LR11_Time2_PERMANOVA = 
  LV_LR11_Time2_PERMANOVA |>
  data.frame() |>
  mutate(F_Time_storage_mo = (Time_storage_mo / 1) / (Residual / 9))

# Calculate p-value
LV_LR11_Time2_PERMANOVA$F_Time_storage_mo[1]
with(LV_LR11_Time2_PERMANOVA, 
     sum(F_Time_storage_mo >= F_Time_storage_mo[1]) / 
       length(F_Time_storage_mo))

set.seed(NULL)

# > > 3.1.11: nMDS and PERMANOVA on raw percentages -----
# Do not include replacements for 0 values.

# > > > 3.1.11.1: nMDS -----
# > > > > 3.1.11.1.1: Fit NMDS ----
LV_nMDS = metaMDS(ma_LV_HPw,
                  k = 2,
                  distance = "euclidean",
                  wascores = F,
                  noshare = F,
                  autotransform = F)
NMDS_stress = LV_nMDS$stress

# View Stressplot
stressplot(LV_nMDS)

# Extract scores from fit
df_LV_nMDS_plot = scores(LV_nMDS, display = c("sites")) # pull scores from MDS fit 
df_LV_nMDS_plot = data.frame(df_LV_HPw[,1:3], df_LV_nMDS_plot)

# get arrows of FAs correlating to first two dimensions
LV_nMDS_envfit = envfit(LV_nMDS, ma_LV_HPw, perm = 99999) 
# fit vector arrows of FAs correlating with dimensions

# pull scores vector arrows of FAs correlating with dimensions
FA_arrows = scores(LV_nMDS_envfit, display = c("vectors")) 

# calculate scaling multiplier based on min range of points in NMDS
scalefactor = 0.35 * min(max(df_LV_nMDS_plot$NMDS1) - min(df_LV_nMDS_plot$NMDS1), 
                         max(df_LV_nMDS_plot$NMDS2) - min(df_LV_nMDS_plot$NMDS2))
# makes scaling multiplier based on min range of points

# scaling of arrows to fit plot area (with vegan's ordiArrowMul function)
#FA_arrows=FA_arrows*ordiArrowMul(LV_nMDS_plot) 

# scaling of arrows to fit plot area based on scalefactor above
FA_arrows = FA_arrows * scalefactor 

# save those as a dataframe
FA_arrows = as.data.frame(FA_arrows) 

# Use onlyFA with over 5% representation for arrows:
FA_arrows = FA_arrows[rownames(FA_arrows) %in% df_LV_HP_list_5,]

# FA arrow text formatting
arrow_angle = 90 - ((atan2(FA_arrows$NMDS1, FA_arrows$NMDS2) * 180) / pi)
arrow_angle = ifelse(arrow_angle > 90 & arrow_angle < 270,
                     arrow_angle + 180,
                     arrow_angle)
arrow_hjust = ifelse(arrow_angle > 90,
                     1,
                     0)
FA_arrows$labs = as.character(rownames(FA_arrows))
FA_arrows$labs = gsub("C", "", FA_arrows$labs)
FA_arrows$labs = gsub("B", "", FA_arrows$labs)
FA_arrows$labs = sub("c", "*c", FA_arrows$labs)
FA_arrows$labs = sub("\\.", ":", FA_arrows$labs)
FA_arrows$labs = sub("n", "*omega*-", FA_arrows$labs)

# Stress annotation
fontfamily = "sans"

stressanot = grobTree(textGrob(paste("stress = ",
                                     round(NMDS_stress, 3),
                                     sep = ''),
                               x = 0.03,
                               y = 0.975,
                               hjust = 0,
                               gp = gpar(fontsize = 6,
                                         fontfamily = fontfamily)))

# Aggregate by type and storage time
df_LV_nMDS_MSD = df_LV_nMDS_plot |>
  group_by(Type, Time_storage_mo) |>
  summarize(NMDS1M = mean(NMDS1),
            NMDS2M = mean(NMDS2),
            NMDS1SD = sd(NMDS1),
            NMDS2SD = sd(NMDS2))

# Create hulls
df_LV_nMDS_hull = df_LV_nMDS_plot |>
  group_by(Type, Time_storage_mo) |>
  slice(chull(NMDS1, NMDS2))

df_LV_nMDS_hull2 = df_LV_nMDS_plot |>
  group_by(Tag_no_Crab_ID) |>
  slice(chull(NMDS1, NMDS2))

# > > > > 3.1.11.1.1: Plot NMDS ----
# > > > > > 3.1.11.1.1.1: By Type, Storage Time ----
plot_LV_HP_nMDS = ggplot()+
  # geom_polygon(aes(x = NMDS1,
  #                  y = NMDS2,
  #                  color = Time_storage_mo,
  #                  linetype = Type,
  #                  group = interaction(Time_storage_mo, Type)),
  #              linewidth = 1,
  #              fill = NA,
  #              data = df_LV_nMDS_hull) +
  geom_point(data = df_LV_nMDS_plot,
             aes(NMDS1,
                 NMDS2,
                 colour = Time_storage_mo,
                 shape = Type),
             size = 2,
             alpha = 0.5)+
  # geom_segment(data = FA_arrows,
  #              aes(x = 0,
  #                  y = 0,
  #                  xend = NMDS1,
  #                  yend = NMDS2),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              colour = "gray50",
  #              alpha = 1)+
  # geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
  #           aes(NMDS1, NMDS2, label =  FA_arrows$labs),
  #           color = "grey50",
#           size = unit(2, "pt"),
#           parse = T,
#           angle = arrow_angle,
#           hjust = arrow_hjust) +
geom_errorbar(aes(x = NMDS1M,
                  ymin = NMDS2M - NMDS2SD,
                  ymax = NMDS2M + NMDS2SD,
                  color = Time_storage_mo),
              width = 0,
              linewidth = 0.5,
              data = df_LV_nMDS_MSD) +
  geom_errorbarh(aes(y = NMDS2M,
                     xmin = NMDS1M - NMDS1SD,
                     xmax = NMDS1M + NMDS1SD,
                     color = Time_storage_mo),
                 height = 0,
                 linewidth = 0.5,
                 data = df_LV_nMDS_MSD) +
  geom_point(aes(x = NMDS1M,
                 y = NMDS2M,
                 color = Time_storage_mo,
                 shape = Type),
             size = 3,
             data = df_LV_nMDS_MSD) +
  scale_color_viridis_c(option = "magma",
                        begin = 0.2,
                        end = 0.9,
                        direction = -1) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("DTS", "FRZ"),
                        name = "Type") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +  
  theme_classic() +
  coord_fixed(xlim = c(-19, 12),
              ylim = c(-16, 15)) +
  labs(x = NULL,
       # x = "nMDS dim 1",
       y = "nMDS dim 2",
       color = "Storage \nTime") +
  annotation_custom(stressanot) +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical")
plot_LV_HP_nMDS

# > > > > > 3.1.11.1.1.2: By Individual ----
plot_LV_HP_nMDS2 = ggplot()+
  geom_polygon(aes(x = NMDS1,
                   y = NMDS2,
                   color = Tag_no_Crab_ID),
               linewidth = 1,
               fill = NA,
               data = df_LV_nMDS_hull2) +
  geom_point(data = df_LV_nMDS_plot,
             aes(NMDS1,
                 NMDS2,
                 colour = Tag_no_Crab_ID,
                 shape = Type),
             size = 2,
             alpha = 0.5)+
  # geom_segment(data = FA_arrows,
  #              aes(x = 0,
  #                  y = 0,
  #                  xend = NMDS1,
  #                  yend = NMDS2),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              colour = "gray50",
  #              alpha = 1)+
  # geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
  #           aes(NMDS1, NMDS2, label =  FA_arrows$labs),
  #           color = "grey50",
#           size = unit(2, "pt"),
#           parse = T,
#           angle = arrow_angle,
#           hjust = arrow_hjust) +
scale_color_viridis_d(option = "mako",
                      begin = 0.2,
                      end = 0.9,
                      direction = -1,
                      guide = "none") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type",
                     guide = "none") +  
  theme_classic() +
  coord_fixed(xlim = c(-19, 12),
              ylim = c(-16, 15)) +
  labs(x = "nMDS dim 1",
       # y = "nMDS dim 2",
       y = NULL,
       color = "Crab ID") +
  annotation_custom(stressanot)+
  theme(legend.position = "none",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical")
plot_LV_HP_nMDS2

# > > > > > 3.1.11.1.1.3: Baseplot with Vectors ----
plot_LV_HP_nMDS_base = 
  ggplot()+
  geom_segment(data = FA_arrows,
               aes(x = 0,
                   y = 0,
                   xend = NMDS1,
                   yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "gray50",
               alpha = 1)+
  geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
            aes(NMDS1, NMDS2, label =  FA_arrows$labs),
            color = "grey50",
            size = unit(2, "pt"),
            parse = T,
            angle = arrow_angle,
            hjust = arrow_hjust) +
  theme_classic() +
  coord_fixed(xlim = c(-19, 12),
              ylim = c(-16, 15)) +
  labs(x = NULL,
       # x = "nMDS dim 1",
       # y = "nMDS dim 2",
       y = NULL) +
  annotation_custom(stressanot) +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical")
plot_LV_HP_nMDS_base

# > > > > > 3.1.11.1.1.4: Combine and Export ----
plot_LV_HP_PCA + 
  plot_LV_HP_PCA2 + 
  plot_LV_HP_PCA_base +
  plot_LV_HP_nMDS +
  plot_LV_HP_nMDS2 +
  plot_LV_HP_nMDS_base +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(20, "pt"))

ggsave("CrabDTS_Fig03.tiff",
       device = "tiff",
       height = 5.75,
       width = 7.5,
       dpi = res)

# > > > 3.1.11.2: PERMANOVA -----

# Use the same analysis strategy as for logratios:
# Sample type and storage time are severely unbalanced; this strategy was 
# chosen for practical sample numbers and costs. To analyze the data, we 
# separate the data into two sets: one to assess Sample type at month 5 and 
# the other to assess Time for DTS samples. 

# Isolate month 5
df_LV_HPw_Type = 
  df_LV_HPw |>
  filter(Time_storage_mo == 5)

# Isolate month 3 DTS and month 0 baseline
df_LV_HPw_Type2 = 
  df_LV_HPw |>
  filter(Time_storage_mo == 2 |
           (Time_storage_mo == 5 &
              Type == "Tissue"))

# Isolate DTS samples, and month 0 frozen as baseline
df_LV_HPw_Time = 
  df_LV_HPw |>
  filter(Time_storage_mo == 0 |
           Type == "DTS")

# Isolate frozen samples
df_LV_HPw_Time2 = 
  df_LV_HPw |>
  filter(Type == "Tissue")

# To control for the repeated-measures nature of the experiment, we constrain 
# permutations based on the individual. 
LV_perc_Type_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HPw_Type$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_perc_Type2_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HPw_Type2$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_perc_Time_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HPw_Time$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

LV_perc_Time2_perm = 
  how(within = Within(type = "free"), # samples are free to shuffle within individual
      plots = Plots(strata = df_LV_HPw_Time2$Tag_no_Crab_ID, 
                    type = "free"), # samples cannot shuffle across individual
      nperm = 9999) 

# We run them once to check our general structure. We do not use this output for
# inference because we know the pseudo-Fs are not calculated how we want.
LV_perc_Type_PERMANOVA_base = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 5 ,] ~ 
                                        Tag_no_Crab_ID +
                                        Type,
                                      data = df_LV_HPw_Type,
                                      method = "euclidean",
                                      permutations = LV_perc_Type_perm, 
                                      by = "terms")

LV_perc_Type_PERMANOVA_base

LV_perc_Type2_PERMANOVA_base = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 2 |
                                                   (df_LV_HPw$Time_storage_mo == 5 &
                                                      df_LV_HPw$Type == "Tissue"),] ~ 
                                         Tag_no_Crab_ID +
                                         Type,
                                       data = df_LV_HPw_Type2,
                                       method = "euclidean",
                                       permutations = LV_perc_Type2_perm, 
                                       by = "terms")

LV_perc_Type2_PERMANOVA_base

LV_perc_Time_PERMANOVA_base = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 0 |
                                                  df_LV_HPw$Type == "DTS" ,] ~ 
                                        Tag_no_Crab_ID +
                                        Time_storage_mo,
                                      data = df_LV_HPw_Time,
                                      method = "euclidean",
                                      permutations = LV_perc_Time_perm, 
                                      by = "terms")

LV_perc_Time_PERMANOVA_base

LV_perc_Time2_PERMANOVA_base = adonis2(ma_LV_HPw[df_LV_HPw$Type == "Tissue" ,] ~ 
                                         Tag_no_Crab_ID +
                                         Time_storage_mo,
                                       data = df_LV_HPw_Time2,
                                       method = "euclidean",
                                       permutations = LV_perc_Time2_perm, 
                                       by = "terms")

LV_perc_Time2_PERMANOVA_base

# Note that for both, we include Tag_no_Crab_ID to allow proper estimation of
# the pseudo-F statistics. These rely on proper denominators relevant to each 
# term. Using the residual as done with adonis2() default is not (always) appropriate 
# (Anderson and Ter Braak 2003 Journal of Statistical Computation and 
# Simulation). Should be based on exchangeable units, which are the individual
# crab. The method used below is from Bakker 2024, Applied Multivariate 
# Statistics in R.

# Check number of possible permutations
numPerms(df_LV_HPw_Type,
         LV_perc_Type_perm)

numPerms(df_LV_HPw_Type2,
         LV_perc_Type2_perm)

numPerms(df_LV_HPw_Time,
         LV_perc_Time_perm)

numPerms(df_LV_HPw_Time2,
         LV_perc_Time2_perm)

# Run PERMANOVA

# Type: 

# F_Type = MS_Type / MS_Residual

# Type is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_perc_Type_perms = rbind(1:nrow(df_LV_HPw_Type),
                           shuffleSet(n = nrow(df_LV_HPw_Type),
                                      control = LV_perc_Type_perm,
                                      nset = 9999))

LV_perc_Type_PERMANOVA = matrix(nrow = nrow(LV_perc_Type_perms),
                                ncol = 4)

colnames(LV_perc_Type_PERMANOVA) = c("Tag_no_Crab_ID",
                                     "Type",
                                     "Residual",
                                     "Total")

for (i in 1:nrow(LV_perc_Type_perms)) {
  df_temp = df_LV_HPw_Type[LV_perc_Type_perms[i, ], ]
  res_temp = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 5 ,] ~
                       Tag_no_Crab_ID +
                       Type,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_perc_Type_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_perc_Type_PERMANOVA = 
  LV_perc_Type_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 9))

# Calculate p-value
LV_perc_Type_PERMANOVA$F_Type[1]
with(LV_perc_Type_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

set.seed(NULL)

# Type2: 

# F_Type = MS_Type / MS_Residual

# Type is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_perc_Type2_perms = rbind(1:nrow(df_LV_HPw_Type2),
                            shuffleSet(n = nrow(df_LV_HPw_Type2),
                                       control = LV_perc_Type2_perm,
                                       nset = 9999))

LV_perc_Type2_PERMANOVA = matrix(nrow = nrow(LV_perc_Type2_perms),
                                 ncol = 4)

colnames(LV_perc_Type2_PERMANOVA) = c("Tag_no_Crab_ID",
                                      "Type",
                                      "Residual",
                                      "Total")

for (i in 1:nrow(LV_perc_Type2_perms)) {
  df_temp = df_LV_HPw_Type2[LV_perc_Type2_perms[i, ], ]
  res_temp = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 2 |
                                 (df_LV_HPw$Time_storage_mo == 5 &
                                    df_LV_HPw$Type == "Tissue"),] ~
                       Tag_no_Crab_ID +
                       Type,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_perc_Type2_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_perc_Type2_PERMANOVA = 
  LV_perc_Type2_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 9))

# Calculate p-value
LV_perc_Type2_PERMANOVA$F_Type[1]
with(LV_perc_Type2_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

set.seed(NULL)

# Time: 

# F_Time_storage_mo = MS_Time_storage_mo / MS_Residual

# Time is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_perc_Time_perms = rbind(1:nrow(df_LV_HPw_Time),
                           shuffleSet(n = nrow(df_LV_HPw_Time),
                                      control = LV_perc_Time_perm,
                                      nset = 9999))

LV_perc_Time_PERMANOVA = matrix(nrow = nrow(LV_perc_Time_perms),
                                ncol = 4)

colnames(LV_perc_Time_PERMANOVA) = c("Tag_no_Crab_ID",
                                     "Time_storage_mo",
                                     "Residual",
                                     "Total")

for (i in 1:nrow(LV_perc_Time_perms)) {
  df_temp = df_LV_HPw_Time[LV_perc_Time_perms[i, ], ]
  res_temp = adonis2(ma_LV_HPw[df_LV_HPw$Time_storage_mo == 0 |
                                 df_LV_HPw$Type == "DTS" ,] ~
                       Tag_no_Crab_ID +
                       Time_storage_mo,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_perc_Time_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_perc_Time_PERMANOVA = 
  LV_perc_Time_PERMANOVA |>
  data.frame() |>
  mutate(F_Time_storage_mo = (Time_storage_mo / 1) / (Residual / 9))

# Calculate p-value
LV_perc_Time_PERMANOVA$F_Time_storage_mo[1]
with(LV_perc_Time_PERMANOVA, 
     sum(F_Time_storage_mo >= F_Time_storage_mo[1]) / 
       length(F_Time_storage_mo))

set.seed(NULL)

# Time2: 

# F_Time_storage_mo = MS_Time_storage_mo / MS_Residual

# Time is at the level of exchangeable units

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a nesting 
# factor/random-ish effect.

set.seed(207)

LV_perc_Time2_perms = rbind(1:nrow(df_LV_HPw_Time2),
                            shuffleSet(n = nrow(df_LV_HPw_Time2),
                                       control = LV_perc_Time2_perm,
                                       nset = 9999))

LV_perc_Time2_PERMANOVA = matrix(nrow = nrow(LV_perc_Time2_perms),
                                 ncol = 4)

colnames(LV_perc_Time2_PERMANOVA) = c("Tag_no_Crab_ID",
                                      "Time_storage_mo",
                                      "Residual",
                                      "Total")

for (i in 1:nrow(LV_perc_Time2_perms)) {
  df_temp = df_LV_HPw_Time2[LV_perc_Time2_perms[i, ], ]
  res_temp = adonis2(ma_LV_HPw[df_LV_HPw$Type == "Tissue" ,] ~
                       Tag_no_Crab_ID +
                       Time_storage_mo,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  LV_perc_Time2_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-F
LV_perc_Time2_PERMANOVA = 
  LV_perc_Time2_PERMANOVA |>
  data.frame() |>
  mutate(F_Time_storage_mo = (Time_storage_mo / 1) / (Residual / 9))

# Calculate p-value
LV_perc_Time2_PERMANOVA$F_Time_storage_mo[1]
with(LV_perc_Time2_PERMANOVA, 
     sum(F_Time_storage_mo >= F_Time_storage_mo[1]) / 
       length(F_Time_storage_mo))

set.seed(NULL)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# > 3.2: Field Validation Data -----
# > > 3.2.1: W/ mat M ----
# Including mature males

# > > > 3.2.1.1: Organize data -----

# Combine Odd and Branched Fatty Acids as Bacterial Markers
# Note: no branched FA, but lumping all the same
# Sum OBCFAs
banana = df_FV_HP |>
  filter(FA %in% OBCFA_list) |>
  group_by(SampleID,
           Set,
           Species,
           Date_sampling,
           Date_extraction,
           Type,
           Tissue,
           CW_mm,
           Tag_no_Crab_ID,
           Sex,
           Dev,
           Shell,
           Chela_Clutch,
           Time_storage_mo,
           TFA_ug) |>
  summarize(perc = sum(perc),
            FA = "OBCFA")

# Attach summed OBCFAs and remove individual FAs
df_FV_HP = 
  bind_rows(df_FV_HP,
            banana) |>
  filter(! FA %in% OBCFA_list)

# Check
sort(unique(df_FV_HP$FA))

# Summarize
df_FV_HP_sum = df_FV_HP |>
  group_by(Type, Sex, Dev, FA) |>
  summarise(M = mean(perc),
            SD = sd(perc))

# Get FAs that are at least 0.5% of any type/sex/maturity combination
df_FV_HP_list_0.5 = unique(df_FV_HP_sum[df_FV_HP_sum$M >= 0.5,]$FA)
df_FV_HP_list_0.5

# Wide format data frame with reduced set
df_FV_HPw = df_FV_HP[df_FV_HP$FA %in% df_FV_HP_list_0.5, ] |>
  pivot_wider(names_from = FA,
              values_from = perc,
              id_cols = c(Tag_no_Crab_ID, Type, Sex, Dev))

ma_FV_HPw = df_FV_HPw[ , 5:length(df_FV_HPw) ]

# Replace 0s with minimum values
ma_FV_HPw_0NA = ma_FV_HPw
ma_FV_HPw_0NA[ma_FV_HPw_0NA == 0] = NA
FAmin = apply(ma_FV_HPw_0NA, 2, min, na.rm = T)

for(j in 1:ncol(ma_FV_HPw)) {
  for(i in 1:nrow(ma_FV_HPw)) {
    if(ma_FV_HPw[i, j] == 0) ma_FV_HPw[i, j] <- 0.5 * FAmin[j]
  }
}

# renormalize to 100% ("re-close")
ma_FV_HPw = ma_FV_HPw / rowSums(ma_FV_HPw) * 100

# > > > 3.2.1.2: Stepwise LR selection -----
# > > > > Step 1 ----
FV_HP_step1 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10)
print.ratios(FV_HP_step1$names.top, FV_HP_step1$R2.top)
# FA1        FA2    R2
# 1  C16.0    C22.6n3 53.73 <- Palmitic vs. DHA
# 2  C16.0   C20.4n6c 53.51
# 3  C16.0    C22.5n3 53.16
# 4  C16.0    C20.5n3 52.47
# 5  C16.0 C20.1n9.7c 51.61
# 6  C14.0      C16.0 50.73
# 7  C16.0      OBCFA 50.31
# 8  C16.0   C22.4n6c 48.02
# 9  C16.0   C18.2n6c 48.01
# 10 C16.0   C20.2n6c 47.99
FV_HP_LR1 = FV_HP_step1$logratios.top[,1]
FV_HP_NR1 = FV_HP_step1$ratios.top[1,]
FV_HP_R21 = FV_HP_step1$R2.top[1]

# > > > > Step 2 -----
FV_HP_step2 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR1)
print.ratios(FV_HP_step2$names.top, FV_HP_step2$R2.top)
# FA1        FA2    R2
# 1     C16.0   C16.1n7c 76.28 <- diatom indicator, sensu Copeman et al.
# 2  C16.1n7c    C22.6n3 76.28
# 3   C20.5n3    C22.6n3 76.27
# 4     C16.0    C20.5n3 76.27
# 5  C16.1n7c C20.1n9.7c 76.27
# 6  C16.1n7c   C22.1n9c 76.19
# 7     C16.0    C18.4n3 76.18
# 8   C18.4n3    C22.6n3 76.18
# 9   C18.4n3   C22.1n9c 76.04
# 10  C16.4n1   C22.1n9c 75.91
FV_HP_LR2 = cbind(FV_HP_LR1, FV_HP_step2$logratios.top[,1])
FV_HP_NR2 = rbind(FV_HP_NR1, FV_HP_step2$ratios.top[1,])
FV_HP_R22 = c(FV_HP_R21, FV_HP_step2$R2.top[1])

# > > > > Step 3 -----
FV_HP_step3 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR2)
print.ratios(FV_HP_step3$names.top, FV_HP_step3$R2.top)
# FA1        FA2    R2
# 1  C18.1n9c    C20.5n3 86.48 <- sep browns, processed material and diatoms, reds
# 2     C16.0   C18.1n9c 86.17
# 3  C16.1n7c   C18.1n9c 86.17
# 4  C18.1n9c    C22.6n3 86.17
# 5  C18.1n9c    C18.4n3 85.86
# 6  C18.1n9c   C20.4n6c 85.57
# 7  C18.1n9c   C18.1n9t 85.56
# 8  C18.1n9c      OBCFA 85.45
# 9  C18.1n9c    C22.5n3 85.29
# 10 C18.1n9c C20.1n9.7c 85.18
FV_HP_LR3 = cbind(FV_HP_LR2, FV_HP_step3$logratios.top[,1])
FV_HP_NR3 = rbind(FV_HP_NR2, FV_HP_step3$ratios.top[1,])
FV_HP_R23 = c(FV_HP_R22, FV_HP_step3$R2.top[1])

# > > > > Step 4 -----
FV_HP_step4 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR3)
print.ratios(FV_HP_step4$names.top, FV_HP_step4$R2.top)
# FA1      FA2    R2
# 1     C18.0 C18.1n9c 91.25 <- SAFA and desat product
# 2     C18.0  C20.5n3 91.25
# 3   C16.4n1    C18.0 91.13
# 4     C18.0  C18.4n3 91.08
# 5     C16.0    C18.0 91.02
# 6     C18.0  C22.6n3 91.02
# 7  C16.1n7c    C18.0 91.02
# 8     C18.0 C20.4n6c 90.94
# 9     C18.0    OBCFA 90.84
# 10    C14.0    C18.0 90.77
FV_HP_LR4 = cbind(FV_HP_LR3, FV_HP_step4$logratios.top[,1])
FV_HP_NR4 = rbind(FV_HP_NR3, FV_HP_step4$ratios.top[1,])
FV_HP_R24 = c(FV_HP_R23, FV_HP_step4$R2.top[1])

# > > > > Step 5 -----
FV_HP_step5 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR4)
print.ratios(FV_HP_step5$names.top, FV_HP_step5$R2.top)
# FA1        FA2    R2
# 1  C20.1n9.7c    C20.5n3 94.58 <- 20:1n9 for copes/zoops and 20:5n3 for diatoms
# 2    C18.1n9c C20.1n9.7c 94.58
# 3       C18.0 C20.1n9.7c 94.58
# 4     C16.4n1   C22.1n9c 94.51
# 5     C18.4n3   C20.2n6c 94.44
# 6     C16.4n1   C20.2n6c 94.43
# 7     C20.5n3   C22.1n9c 94.39
# 8       C18.0   C22.1n9c 94.39
# 9    C18.1n9c   C22.1n9c 94.39
# 10   C18.1n9c   C20.2n6c 94.39
FV_HP_LR5 = cbind(FV_HP_LR4, FV_HP_step5$logratios.top[,1])
FV_HP_NR5 = rbind(FV_HP_NR4, FV_HP_step5$ratios.top[1,])
FV_HP_R25 = c(FV_HP_R24, FV_HP_step5$R2.top[1])

# > > > > Step 6 -----
FV_HP_step6 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR5)
print.ratios(FV_HP_step6$names.top, FV_HP_step6$R2.top)
# FA1      FA2    R2
# 1     C16.4n1 C22.1n9c 96.50
# 2     C18.4n3 C22.1n9c 96.40 <- <- SDA and 22:1n9 for copes/zoops
# 3    C18.1n9t C22.1n9c 96.38
# 4     C20.5n3 C22.1n9c 96.33
# 5       C18.0 C22.1n9c 96.33
# 6  C20.1n9.7c C22.1n9c 96.33
# 7    C18.1n9c C22.1n9c 96.33
# 8    C22.1n9c  C22.6n3 96.31
# 9    C16.1n7c C22.1n9c 96.31
# 10      C16.0 C22.1n9c 96.31
FV_HP_LR6 = cbind(FV_HP_LR5, FV_HP_step6$logratios.top[,2])
FV_HP_NR6 = rbind(FV_HP_NR5, FV_HP_step6$ratios.top[2,])
FV_HP_R26 = c(FV_HP_R25, FV_HP_step6$R2.top[2])

# > > > > Step 7 -----
FV_HP_step7 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR6)
print.ratios(FV_HP_step7$names.top, FV_HP_step7$R2.top)
# FA1        FA2    R2
# 1    C14.0    C22.5n3 97.48
# 2    C14.0    C22.6n3 97.47 <- myristic common SAFA, DHA more interpretable than DPA
# 3    C14.0      C16.0 97.47
# 4    C14.0   C16.1n7c 97.47
# 5    C14.0   C22.4n6c 97.46
# 6    C14.0   C18.1n9t 97.44
# 7  C16.4n1      OBCFA 97.44
# 8    C14.0   C22.1n9c 97.42
# 9    C14.0    C18.4n3 97.42
# 10   C14.0 C20.1n9.7c 97.42
FV_HP_LR7 = cbind(FV_HP_LR6, FV_HP_step7$logratios.top[,2])
FV_HP_NR7 = rbind(FV_HP_NR6, FV_HP_step7$ratios.top[2,])
FV_HP_R27 = c(FV_HP_R26, FV_HP_step7$R2.top[2])

# > > > > Step 8 -----
FV_HP_step8 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR7)
print.ratios(FV_HP_step8$names.top, FV_HP_step8$R2.top)
# FA1        FA2    R2
# 1   C16.4n1      OBCFA 98.50 <- 16:4n1 unclear, OBCFA for bacteria
# 2  C20.4n6c   C22.1n9c 98.47
# 3   C18.4n3   C20.4n6c 98.47
# 4   C16.4n1   C20.4n6c 98.47
# 5   C16.4n1 C20.1n9.7c 98.45
# 6   C16.4n1   C18.1n9c 98.45
# 7   C16.4n1      C18.0 98.45
# 8   C16.4n1    C20.5n3 98.45
# 9   C16.4n1    C22.5n3 98.42
# 10 C22.1n9c    C22.5n3 98.42
FV_HP_LR8 = cbind(FV_HP_LR7, FV_HP_step8$logratios.top[,1])
FV_HP_NR8 = rbind(FV_HP_NR7, FV_HP_step8$ratios.top[1,])
FV_HP_R28 = c(FV_HP_R27, FV_HP_step8$R2.top[1])

# > > > > Step 9 -----
FV_HP_step9 = STEP(data = ma_FV_HPw, nsteps = 1, top = 10, previous = FV_HP_LR8)
print.ratios(FV_HP_step9$names.top, FV_HP_step9$R2.top)
# FA1     FA2    R2
# 1     C16.0   OBCFA 98.98
# 2     C16.0 C16.4n1 98.98
# 3  C16.1n7c C16.4n1 98.98
# 4     C14.0 C16.4n1 98.98
# 5   C22.6n3   OBCFA 98.98
# 6     C14.0   OBCFA 98.98
# 7  C16.1n7c   OBCFA 98.98
# 8   C16.4n1 C22.6n3 98.98
# 9     C16.0 C22.5n3 98.93
# 10    C14.0 C22.5n3 98.93

# Stop here. Recovering less than 1% variation and prominent FA (16:4n1) has 
# limited interpretability.

# > > > 3.2.1.3: Organize step 8 LRs -----
rownames(FV_HP_NR8) = paste("Step", 1:8, sep="")
colnames(FV_HP_NR8) = c("FA1","FA2")
FV_HP_FR = as.data.frame(cbind(FV_HP_NR8,
                               Ratio = paste(colnames(ma_FV_HPw)[FV_HP_NR8[,1]],
                                             "/",
                                             colnames(ma_FV_HPw)[FV_HP_NR8[,2]],
                                             sep="")))
FV_HP_FR$FA1_lab = colnames(ma_FV_HPw)[FV_HP_NR8[,1]]
FV_HP_FR$FA2_lab = colnames(ma_FV_HPw)[FV_HP_NR8[,2]]
FV_HP_FR$step = 1:8

FV_HP_FR$R2_cum = FV_HP_R28
FV_HP_FR = FV_HP_FR |>
  arrange(step) |>
  mutate(R2 = R2_cum - lag(R2_cum,
                           default = 0))

FV_HP_FR

# FV_HP_FR_export = FV_HP_FR
# FV_HP_FR_export$R2_cum = round(FV_HP_FR_export$R2_cum, 3)
# FV_HP_FR_export$R2 = round(FV_HP_FR_export$R2, 3)
# FV_HP_FR_export
# write.csv(FV_HP_FR_export, "FV_HP_FR_export.csv")

colnames(FV_HP_LR8) = FV_HP_FR[,3]

# FV_HP_PiR = sort(unique(as.numeric(FV_HP_NR8))) # "parts in ratios"
FV_HP_PiR = unique(as.numeric(t(FV_HP_NR8))) # "parts in ratios"
colnames(ma_FV_HPw)[FV_HP_PiR]

# > > > 3.2.1.4: Prep acyclic graph -----
FV_HP_PiR_dim = data.frame(FA = colnames(ma_FV_HPw)[FV_HP_PiR],
                           dim1 = c(2, # C16.0
                                    3, # C22.6n3
                                    1, # C16.1n7c
                                    2, # C18.1n9c
                                    3, # C20.5n3
                                    1, # C18.0
                                    4, # C20.1n9.7c
                                    1, # C18.4n3
                                    2, # C22.1n9c
                                    4, # C14.0
                                    3, # C16.4n1
                                    4),# OBCFA
                           dim2 = c(3, # C16.0
                                    3, # C22:6n3
                                    3, # C16.1n7c
                                    2, # C18.1n9c
                                    2, # C20.5n3
                                    2, # C18.0
                                    2, # C20.1n9.7c
                                    1, # C18.4n3
                                    1, # C22.1n9c
                                    3, # C14.0
                                    1, # C16.4n1
                                    1)) # OBCFA

FV_HP_PiR_dim$FA_labs = str_replace_all(FV_HP_PiR_dim$FA,
                                        FA_labeller)

FV_HP_FR = left_join(FV_HP_FR,
                     FV_HP_PiR_dim,
                     by = c("FA1_lab" = "FA"))
FV_HP_FR = left_join(FV_HP_FR,
                     FV_HP_PiR_dim,
                     by = c("FA2_lab" = "FA"))

# > > > 3.2.1.5: Prep PCA -----
FV_HP_PCA = PCA(FV_HP_LR8, weight = FALSE)

PLOT.PCA(FV_HP_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

# Extract and calculate contribution vectors
axes_inv = c(1, 1)
dim = c(1, 2)
FV_HP_RPC = FV_HP_PCA$rowcoord[, dim] %*% diag(FV_HP_PCA$sv[dim] * axes_inv)
FV_HP_CSC = FV_HP_PCA$colcoord[, dim] %*% diag(axes_inv)
FV_HP_CPC = FV_HP_CSC %*% diag(FV_HP_PCA$sv[dim])
FV_HP_CCC = FV_HP_CSC * sqrt(FV_HP_PCA$colmass)
FV_HP_CRD = FV_HP_CCC

# rescale contribution vectors for visibility
vectorscale = 1
FV_HP_CRD_scale = FV_HP_CRD * vectorscale

# Extract points
df_FV_HP_RPC = cbind(df_FV_HPw[, 1:4], FV_HP_RPC)
colnames(df_FV_HP_RPC) = c(colnames(df_FV_HPw[, 1:4]),'dim1','dim2')

# Create SexDev
df_FV_HP_RPC$SexDev = paste(df_FV_HP_RPC$Sex,
                            df_FV_HP_RPC$Dev,
                            sep = "_")

df_FV_HP_RPC$SexDev = factor(df_FV_HP_RPC$SexDev,
                             levels = c("F_imm",
                                        "F_mat",
                                        "M_imm",
                                        "M_mat"))

# Aggregate by type, sex, and maturity
df_FV_HP_RPC_MSD = df_FV_HP_RPC |>
  group_by(Type, SexDev) |>
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SD = sd(dim1),
            dim2SD = sd(dim2))

# Calculate variance explained by PCA dimensions
FV_HP_PCA_perc_1 = 100 * FV_HP_PCA$sv[dim[1]]^2 / sum(FV_HP_PCA$sv^2)
FV_HP_PCA_perc_2 = 100 * FV_HP_PCA$sv[dim[2]]^2 / sum(FV_HP_PCA$sv^2)

df_FV_HP_RPC$ID2 = paste(df_FV_HP_RPC$Type,
                         df_FV_HP_RPC$Sex,
                         df_FV_HP_RPC$Mat,
                         df_FV_HP_RPC$Tag_no_Crab_ID,
                         sep = ".")

# Format arrows
df_FV_HP_arrows = data.frame(dim1 = FV_HP_CRD_scale[ ,1],
                             dim2 = FV_HP_CRD_scale[, 2],
                             LR = FV_HP_PCA$colnames)

df_FV_HP_arrows$step = 1:8

df_FV_HP_arrows$angle = 90 - ((atan2(df_FV_HP_arrows$dim1,
                                     df_FV_HP_arrows$dim2) *
                                 180) /
                                pi)
df_FV_HP_arrows$angle = ifelse(df_FV_HP_arrows$angle > 90 & df_FV_HP_arrows$angle < 270,
                               df_FV_HP_arrows$angle + 180,
                               df_FV_HP_arrows$angle)
df_FV_HP_arrows$hjust = ifelse(df_FV_HP_arrows$angle > 90,
                               1,
                               0)

# Format arrow labels
df_FV_HP_arrows$labs = df_FV_HP_arrows$LR
df_FV_HP_arrows$labs = str_replace_all(df_FV_HP_arrows$labs,
                                       FA_labeller)
df_FV_HP_arrows$labs = str_replace_all(df_FV_HP_arrows$labs,
                                       "/",
                                       " / ")

# establish levels for IDs to keep colors consistent
df_FV_HP_RPC$Tag_no_Crab_ID = factor(df_FV_HP_RPC$Tag_no_Crab_ID,
                                     levels = paste0("FV", sprintf("%02.0f", 1:20)))

# Create hulls
df_FV_HP_RPC_hull = df_FV_HP_RPC |>
  group_by(Type, SexDev) |>
  slice(chull(dim1, dim2))

df_FV_HP_RPC_hull2 = df_FV_HP_RPC |>
  group_by(Tag_no_Crab_ID) |>
  slice(chull(dim1, dim2))

# > > > 3.2.1.6: Plot PCA -----
plot_FV_HP_PCA = ggplot() +
  geom_polygon(aes(x = dim1,
                   y = dim2,
                   color = SexDev,
                   linetype = Type,
                   group = interaction(SexDev, Type)),
               linewidth = 1,
               fill = NA,
               data = df_FV_HP_RPC_hull) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = SexDev,
                 shape = Type),
             alpha = 0.5,
             size = 2,
             data = df_FV_HP_RPC) +
  geom_point(aes(x = ifelse(dim1 > 0,
                            dim1 + 0.075 * cos(angle * pi/180),
                            dim1 - 0.075 * cos(angle * pi/180)),
                 y = ifelse(dim1 > 0,
                            dim2 + 0.075 * sin(angle * pi/180),
                            dim2 - 0.075 * sin(angle * pi/180))),
             colour = "grey50",
             size = 4,
             data = df_FV_HP_arrows)+
  geom_text(aes(x = ifelse(dim1 > 0,
                           dim1 + 0.075 * cos(angle * pi/180),
                           dim1 - 0.075 * cos(angle * pi/180)),
                y = ifelse(dim1 > 0,
                           dim2 + 0.075 * sin(angle * pi/180),
                           dim2 - 0.075 * sin(angle * pi/180)),
                label = step),
            colour = "white",
            size = unit(2, "pt"),
            data = df_FV_HP_arrows)+
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FV_HP_CRD_scale[ ,1],
                   yend = FV_HP_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_errorbar(aes(x = dim1M,
                    ymin = dim2M - dim2SD,
                    ymax = dim2M + dim2SD,
                    color = SexDev),
                width = 0,
                linewidth = 1,
                data = df_FV_HP_RPC_MSD) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - dim1SD,
                     xmax = dim1M + dim1SD,
                     color = SexDev),
                 height = 0,
                 linewidth = 1,
                 data = df_FV_HP_RPC_MSD) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = SexDev,
                 shape = Type),
             size = 5,
             data = df_FV_HP_RPC_MSD) +
  scale_color_viridis_d(option = "magma",
                        begin = 0.2,
                        end = 0.9,
                        direction = -1,
                        labels = c("Female, imm.",
                                   "Female, mat.",
                                   "Male, imm.",
                                   "Male, mat.")) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("DTS", "FRZ"),
                        name = "Type") +
  scale_shape_manual(values = c(18, 20),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FV_HP_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FV_HP_PCA_perc_2, 1),
                 "%)",
                 sep=""),
       color = "Sex, maturity") +
  coord_fixed(
    # xlim = c(-1, 0.8),
    # ylim = c(-0.8, 1)
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(8, "pt")),
        legend.text = element_text(size = unit(8, "pt")),
        axis.title = element_text(size = unit(10, "pt")),
        axis.text = element_text(size = unit(8, "pt")),
        legend.box = "vertical",
        axis.text.y.right = element_text(colour="grey50",
                                         size = unit(8, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(8, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50")) +
  guides(color = guide_legend(override.aes = list(size = 2),
                              nrow = 2),
         shape = guide_legend(override.aes = list(size = 2)))

plot_FV_HP_PCA

plot_FV_HP_PCA2 = ggplot() +
  geom_polygon(aes(x = dim1,
                   y = dim2,
                   color = Tag_no_Crab_ID,
                   group = Tag_no_Crab_ID),
               linewidth = 1,
               fill = NA,
               data = df_FV_HP_RPC_hull2) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = Tag_no_Crab_ID,
                 shape = Type),
             alpha = 0.5,
             size = 3,
             data = df_FV_HP_RPC) +
  # geom_text(aes(x = dim1,
  #                y = dim2,
  #                color = Tag_no_Crab_ID,
  #               label = Tag_no_Crab_ID),
  #            alpha = 0.5,
  #            size = 3,
  #            data = df_FV_HP_RPC) +
  geom_point(aes(x = ifelse(dim1 > 0,
                            dim1 + 0.075 * cos(angle * pi/180),
                            dim1 - 0.075 * cos(angle * pi/180)),
                 y = ifelse(dim1 > 0,
                            dim2 + 0.075 * sin(angle * pi/180),
                            dim2 - 0.075 * sin(angle * pi/180))),
             colour = "grey50",
             size = 4,
             data = df_FV_HP_arrows)+
  geom_text(aes(x = ifelse(dim1 > 0,
                           dim1 + 0.075 * cos(angle * pi/180),
                           dim1 - 0.075 * cos(angle * pi/180)),
                y = ifelse(dim1 > 0,
                           dim2 + 0.075 * sin(angle * pi/180),
                           dim2 - 0.075 * sin(angle * pi/180)),
                label = step),
            colour = "white",
            size = unit(2, "pt"),
            data = df_FV_HP_arrows)+
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FV_HP_CRD_scale[ ,1],
                   yend = FV_HP_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_viridis_d(option = "turbo",
                        begin = 0.2,
                        end = 0.9,
                        direction = -1,
                        guide = "none") +
  scale_shape_manual(values = c(18, 20),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FV_HP_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FV_HP_PCA_perc_2, 1),
                 "%)",
                 sep=""),
       color = "Crab ID") +
  coord_fixed(
    # xlim = c(-1, 0.8),
    # ylim = c(-0.8, 1)
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(8, "pt")),
        legend.text = element_text(size = unit(8, "pt")),
        axis.title = element_text(size = unit(10, "pt")),
        axis.text = element_text(size = unit(8, "pt")),
        legend.box = "vertical",
        axis.text.y.right = element_text(colour="grey50",
                                         size = unit(8, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(8, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50")) +
  guides(shape = guide_legend(override.aes = list(size = 2)))

plot_FV_HP_PCA2

# > > > 3.2.1.7: Plot acyclic graph -----
plot_FV_HP_LR_graph = ggplot() +
  geom_segment(aes(x = dim1.x,
                   y = dim2.x,
                   xend = dim1.y,
                   yend = dim2.y,
                   color = step),
               linewidth = 1,
               data = FV_HP_FR) +
  geom_label(aes(x = dim1,
                 y = dim2,
                 label = FA_labs),
             size = unit(3, "pt"),
             data = FV_HP_PiR_dim) +
  geom_point(aes(x = (dim1.x + dim1.y)/2,
                 y = (dim2.x + dim2.y)/2,
                 color = step),
             size = 6,
             data = FV_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2,
                y = (dim2.x + dim2.y)/2,
                label = step),
            color = "white",
            size = unit(3, "pt"),
            data = FV_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2 + 0.125,
                y = (dim2.x + dim2.y)/2 - 0.125,
                label = round(R2*100, 0),
                color = step),
            size = unit(2, "pt"),
            data = FV_HP_FR) +
  scale_color_gradient(high = "dodgerblue4",
                       low = "dodgerblue1") +
  # coord_fixed(xlim = c(0.25, 4.75),
  #             clip = "off") +
  coord_cartesian(xlim = c(0.5, 4.5),
                  clip = "off") +
  theme_void() +
  theme(legend.position = "none")
plot_FV_HP_LR_graph 

plot_FV_HP_LR_graph | plot_FV_HP_PCA | plot_FV_HP_PCA2 +
  plot_annotation(tag_levels = "A") 

# > > > 3.2.1.8: Export acyclic graph and PCA -----
# tiff("plot_C_opilio_FV_HP_LR_doc.tiff",
#      height = res * 6.5,
#      width = res * 7.5,
#      res = res)
# plot_FV_HP_LR_graph / (plot_FV_HP_PCA | plot_FV_HP_PCA2) +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(heights = c(1, 2))
# dev.off()

# > > > 3.2.1.9: Univariate LR plot -----
df_FV_HP_LR8w = df_FV_HPw[,1:4]
df_FV_HP_LR8w$SexDev = paste(df_FV_HP_LR8w$Sex,
                             df_FV_HP_LR8w$Dev,
                             sep = "_")

df_FV_HP_LR8w$SexDev = factor(df_FV_HP_LR8w$SexDev,
                              levels = c("F_imm",
                                         "F_mat",
                                         "M_imm",
                                         "M_mat"))

df_FV_HP_LR8w = cbind(df_FV_HP_LR8w,
                      FV_HP_LR8)

df_FV_HP_LR8 = df_FV_HP_LR8w |>
  pivot_longer(cols = 6:length(df_FV_HP_LR8w),
               names_to = "FAs",
               values_to = "LR")

df_FV_HP_LR8$FAs = str_replace_all(df_FV_HP_LR8$FAs,
                                   FA_labeller)

df_FV_HP_LR8$FAs = factor(df_FV_HP_LR8$FAs,
                          levels = str_replace_all(FV_HP_FR$Ratio,
                                                   FA_labeller))
df_FV_HP_LR8_sum = df_FV_HP_LR8 |>
  group_by(Type, SexDev, FAs) |>
  summarise(M = mean(LR),
            SD = sd(LR))

# > > > > 3.2.1.9.1: Plot -----
plot_FV_HP_LR_cat = ggplot() +
  geom_hline(aes(yintercept = y),
             linetype = 2,
             color = "grey75",
             data = labels_logratios) +
  geom_errorbar(aes(x = SexDev,
                    color = Type,
                    ymin = M - SD,
                    ymax = M + SD),
                width = 0,
                alpha = 0.5,
                position = position_dodge(0.5),
                data = df_FV_HP_LR8_sum) +
  geom_point(aes(x = SexDev,
                 color = Type,
                 y = M),
             size = 2,
             alpha = 0.5,
             shape = 9,
             position = position_dodge(0.5),
             data = df_FV_HP_LR8_sum) +
  geom_point(aes(x = SexDev,
                 y = LR,
                 color = Type),
             size = 0.5,
             position = position_dodge(0.5),
             data = df_FV_HP_LR8) +
  scale_color_manual(values = c("firebrick",
                                "dodgerblue")) +
  scale_x_discrete(labels = c("Female, imm.",
                              "Female, mat.",
                              "Male, imm.",
                              "Male, mat.")) +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),
                                         breaks = exp(labels_logratios$y),
                                         labels = labels_logratios$label,
                                         name = "Ratio"),
                     limits = c(floor(min(df_FV_HP_LR8$LR) / 0.5) * 0.5,
                                ceiling(max(df_FV_HP_LR8$LR) / 0.5) * 0.5)) +
  labs(x = "Sex, Maturity",
       y = "Logratio",
       color = "Sample Type") +
  theme_classic() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text = element_text(size = unit(4, "pt")),
        axis.title = element_text(size = unit(6, "pt")),
        axis.text.y.right = element_text(color="grey70"),
        axis.ticks.y.right = element_line(color="grey70"),
        axis.line.y.right = element_line(color="grey70"),
        axis.title.y.right = element_text(color = "grey70"),
        legend.background = element_blank(),
        legend.text = element_text(size = unit(4, "pt")),
        legend.title = element_text(size = unit(6, "pt")),
        legend.spacing = unit(0, "pt"),
        strip.text = element_text(size = unit(4, "pt"))) +
  facet_grid(.~ FAs)
plot_FV_HP_LR_cat

# > > > > 3.2.1.9.2: Export -----
# tiff("plot_C_opilio_FV_HP_LR_cat.tiff",
#      height = res * 4,
#      width = res * 7.5,
#      res = res)
# plot_FV_HP_LR_cat
# dev.off()

# We do not formally analyze the data including mature males because mature
# males lack frozen samples for hepatopancreas.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# > > 3.2.2: W/o mat M ----
# Including mature males
# > > > 3.2.2.1 : Organize data ----

# Subset to exclude mature males and include only HP
df_FVnMM_HP = df_FV[df_FV$Tissue == "HP" &
                      !(df_FV$Sex == "M" &
                          df_FV$Dev == "mat"),]
summary(interaction(df_FVnMM_HP$Sex, df_FVnMM_HP$Dev))

# Combine Odd and Branched Fatty Acids as Bacterial Markers
# Note: no branched FA, but lumping all the same
# Sum OBCFAs
banana = df_FVnMM_HP |>
  filter(FA %in% OBCFA_list) |>
  group_by(SampleID,
           Set,
           Species,
           Date_sampling,
           Date_extraction,
           Type,
           Tissue,
           CW_mm,
           Tag_no_Crab_ID,
           Sex,
           Dev,
           Shell,
           Chela_Clutch,
           Time_storage_mo,
           TFA_ug) |>
  summarize(perc = sum(perc),
            FA = "OBCFA")

# Attach summed OBCFAs and remove individual FAs
df_FVnMM_HP = 
  bind_rows(df_FVnMM_HP,
            banana) |>
  filter(! FA %in% OBCFA_list)

# Check
sort(unique(df_FVnMM_HP$FA))

# Summarize
df_FVnMM_HP_sum = df_FVnMM_HP |>
  group_by(Type, Sex, Dev, FA) |>
  summarise(M = mean(perc),
            SD = sd(perc))

# Get FAs that are at least 0.5% of any type/sex/maturity combination
df_FVnMM_HP_list_0.5 = unique(df_FVnMM_HP_sum[df_FVnMM_HP_sum$M >= 0.5,]$FA)
df_FVnMM_HP_list_0.5
# Get FAs that are at least 5% of any type/sex/maturity combination
df_FVnMM_HP_list_5 = unique(df_FVnMM_HP_sum[df_FVnMM_HP_sum$M >= 5,]$FA)
df_FVnMM_HP_list_5

# Wide format data frame with reduced set
df_FVnMM_HPw = df_FVnMM_HP[df_FVnMM_HP$FA %in% df_FVnMM_HP_list_0.5, ] |>
  pivot_wider(names_from = FA,
              values_from = perc,
              id_cols = c(Tag_no_Crab_ID, Type, Sex, Dev))

ma_FVnMM_HPw = df_FVnMM_HPw[ , 5:length(df_FVnMM_HPw) ]

# Replace 0s with minimum values
ma_FVnMM_HPw_0NA = ma_FVnMM_HPw
ma_FVnMM_HPw_0NA[ma_FVnMM_HPw_0NA == 0] = NA
FAmin = apply(ma_FVnMM_HPw_0NA, 2, min, na.rm = T)

for(j in 1:ncol(ma_FVnMM_HPw)) {
  for(i in 1:nrow(ma_FVnMM_HPw)) {
    if(ma_FVnMM_HPw[i, j] == 0) ma_FVnMM_HPw[i, j] <- 0.5 * FAmin[j]
  }
}

# renormalize to 100% ("re-close")
ma_FVnMM_HPw = ma_FVnMM_HPw / rowSums(ma_FVnMM_HPw) * 100

# > > > 3.2.2.2 : Stepwise LR selection ----
# > > > > Step 1 -----
FVnMM_HP_step1 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10)
print.ratios(FVnMM_HP_step1$names.top, FVnMM_HP_step1$R2.top)
# FA1        FA2    R2
# 1  C18.1n9c    C20.5n3 47.45 <- Oleic vs EPA
# 2  C16.1n7c   C18.1n9c 47.34
# 3  C18.1n9c    C22.5n3 46.57
# 4  C18.1n9c   C20.4n6c 45.48
# 5  C18.1n9c    C22.6n3 45.43
# 6  C18.1n9c      OBCFA 45.39
# 7  C18.1n9c    C18.4n3 44.82
# 8  C18.1n9c   C18.1n9t 44.76
# 9     C16.0   C18.1n9c 44.57
# 10 C18.1n9c C20.1n9.7c 43.46
FVnMM_HP_LR1 = FVnMM_HP_step1$logratios.top[,1]
FVnMM_HP_NR1 = FVnMM_HP_step1$ratios.top[1,]
FVnMM_HP_R21 = FVnMM_HP_step1$R2.top[1]

# > > > > Step 2 -----
FVnMM_HP_step2 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR1)
print.ratios(FVnMM_HP_step2$names.top, FVnMM_HP_step2$R2.top)
# FA1      FA2    R2
# 1     C18.0 C18.1n9c 66.10 <- stearic and oleic, desat product
# 2     C18.0  C20.5n3 66.10
# 3     C18.0 C20.4n6c 65.74
# 4  C16.1n7c    C18.0 65.65
# 5     C18.0  C22.5n3 65.54
# 6     C18.0  C22.6n3 65.45
# 7     C18.0  C18.4n3 65.37
# 8     C18.0    OBCFA 65.35
# 9     C18.0 C18.1n9t 65.13
# 10    C16.0    C18.0 65.11
FVnMM_HP_LR2 = cbind(FVnMM_HP_LR1, FVnMM_HP_step2$logratios.top[,1])
FVnMM_HP_NR2 = rbind(FVnMM_HP_NR1, FVnMM_HP_step2$ratios.top[1,])
FVnMM_HP_R22 = c(FVnMM_HP_R21, FVnMM_HP_step2$R2.top[1])

# > > > > Step 3 -----
FVnMM_HP_step3 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR2)
print.ratios(FVnMM_HP_step3$names.top, FVnMM_HP_step3$R2.top)
# FA1        FA2    R2
# 1       C18.0 C20.1n9.7c 77.97
# 2  C20.1n9.7c    C20.5n3 77.97 <- 20:1n9/7 for copes/zoops, EPA for diatoms
# 3    C18.1n9c C20.1n9.7c 77.97
# 4       C18.0   C18.1n9t 77.84
# 5    C18.1n9c   C18.1n9t 77.84
# 6    C18.1n9t    C20.5n3 77.84
# 7  C20.1n9.7c   C20.4n6c 76.50
# 8  C20.1n9.7c    C22.6n3 75.92
# 9    C16.1n7c      C18.0 75.40
# 10   C16.1n7c   C18.1n9c 75.40
FVnMM_HP_LR3 = cbind(FVnMM_HP_LR2, FVnMM_HP_step3$logratios.top[,2])
FVnMM_HP_NR3 = rbind(FVnMM_HP_NR2, FVnMM_HP_step3$ratios.top[2,])
FVnMM_HP_R23 = c(FVnMM_HP_R22, FVnMM_HP_step3$R2.top[2])

# > > > > Step 4 -----
FVnMM_HP_step4 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR3)
print.ratios(FVnMM_HP_step4$names.top, FVnMM_HP_step4$R2.top)
# FA1      FA2    R2
# 1     C18.4n3 C22.1n9c 83.81 <- SDA for cryptophytes, 22:1n9 for copes/zoops
# 2     C16.4n1 C22.1n9c 83.76
# 3    C16.1n7c C22.1n9c 83.74
# 4       C16.0 C22.1n9c 83.69
# 5    C18.1n9t C22.1n9c 83.59
# 6     C20.5n3 C22.1n9c 83.32
# 7       C18.0 C22.1n9c 83.32
# 8    C18.1n9c C22.1n9c 83.32
# 9  C20.1n9.7c C22.1n9c 83.32
# 10   C16.1n7c  C22.6n3 83.16
FVnMM_HP_LR4 = cbind(FVnMM_HP_LR3, FVnMM_HP_step4$logratios.top[,1])
FVnMM_HP_NR4 = rbind(FVnMM_HP_NR3, FVnMM_HP_step4$ratios.top[1,])
FVnMM_HP_R24 = c(FVnMM_HP_R23, FVnMM_HP_step4$R2.top[1])

# > > > > Step 5 -----
FVnMM_HP_step5 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR4)
print.ratios(FVnMM_HP_step5$names.top, FVnMM_HP_step5$R2.top)
# FA1      FA2    R2
# 1  C16.1n7c C18.1n9t 88.39
# 2     C14.0  C22.6n3 88.35
# 3     C16.0  C16.4n1 88.25
# 4     C16.0 C22.1n9c 88.25
# 5     C16.0  C18.4n3 88.25
# 6   C22.6n3    OBCFA 88.24 <- DHA for dinos, OBCFA for bacteria
# 7     C14.0  C22.5n3 88.18
# 8  C20.4n6c  C22.6n3 88.15
# 9   C22.5n3    OBCFA 88.13
# 10    C14.0    C18.0 88.08
FVnMM_HP_LR5 = cbind(FVnMM_HP_LR4, FVnMM_HP_step5$logratios.top[,6])
FVnMM_HP_NR5 = rbind(FVnMM_HP_NR4, FVnMM_HP_step5$ratios.top[6,])
FVnMM_HP_R25 = c(FVnMM_HP_R24, FVnMM_HP_step5$R2.top[6])

# > > > > Step 6 -----
FVnMM_HP_step6 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR5)
print.ratios(FVnMM_HP_step6$names.top, FVnMM_HP_step6$R2.top)
# FA1        FA2    R2
# 1     C16.0   C22.1n9c 92.55
# 2     C16.0    C18.4n3 92.55 <- palmitic common SAFA, SDA for browns, cryptos
# 3     C16.0    C16.4n1 92.47
# 4   C20.5n3   C22.1n9c 91.90
# 5   C18.4n3 C20.1n9.7c 91.90
# 6  C18.1n9c    C18.4n3 91.90
# 7     C18.0   C22.1n9c 91.90
# 8     C18.0    C18.4n3 91.90
# 9  C18.1n9c   C22.1n9c 91.90
# 10  C18.4n3    C20.5n3 91.90
FVnMM_HP_LR6 = cbind(FVnMM_HP_LR5, FVnMM_HP_step6$logratios.top[,2])
FVnMM_HP_NR6 = rbind(FVnMM_HP_NR5, FVnMM_HP_step6$ratios.top[2,])
FVnMM_HP_R26 = c(FVnMM_HP_R25, FVnMM_HP_step6$R2.top[2])

# > > > > Step 7 -----
FVnMM_HP_step7 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR6)
print.ratios(FVnMM_HP_step7$names.top, FVnMM_HP_step7$R2.top)
# FA1      FA2    R2
# 1  C22.1n9c  C22.6n3 95.29
# 2     C16.0    OBCFA 95.29
# 3   C18.4n3    OBCFA 95.29
# 4  C22.1n9c    OBCFA 95.29
# 5   C18.4n3  C22.6n3 95.29 <- cryptos vs dinos
# 6     C16.0  C22.6n3 95.29
# 7     C16.0 C20.4n6c 95.14
# 8  C20.4n6c C22.1n9c 95.14
# 9   C18.4n3 C20.4n6c 95.14
# 10 C16.1n7c C20.4n6c 95.11
FVnMM_HP_LR7 = cbind(FVnMM_HP_LR6, FVnMM_HP_step7$logratios.top[,5])
FVnMM_HP_NR7 = rbind(FVnMM_HP_NR6, FVnMM_HP_step7$ratios.top[5,])
FVnMM_HP_R27 = c(FVnMM_HP_R26, FVnMM_HP_step7$R2.top[5])

# > > > > Step 8 -----
FVnMM_HP_step8 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR7)
print.ratios(FVnMM_HP_step8$names.top, FVnMM_HP_step8$R2.top)
# FA1        FA2    R2
# 1  C16.1n7c   C22.1n9c 97.26
# 2  C16.1n7c    C18.4n3 97.26
# 3  C16.1n7c    C22.6n3 97.26
# 4  C16.1n7c      OBCFA 97.26
# 5     C16.0   C16.1n7c 97.26 <- diatom marker sensu Copeman et al. 
# 6  C16.1n7c    C16.4n1 97.07
# 7  C16.1n7c C20.1n9.7c 97.01
# 8  C16.1n7c      C18.0 97.01
# 9  C16.1n7c    C20.5n3 97.01
# 10 C16.1n7c   C18.1n9c 97.01
FVnMM_HP_LR8 = cbind(FVnMM_HP_LR7, FVnMM_HP_step8$logratios.top[,5])
FVnMM_HP_NR8 = rbind(FVnMM_HP_NR7, FVnMM_HP_step8$ratios.top[5,])
FVnMM_HP_R28 = c(FVnMM_HP_R27, FVnMM_HP_step8$R2.top[5])

# > > > > Step 9 -----
FVnMM_HP_step9 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR8)
print.ratios(FVnMM_HP_step9$names.top, FVnMM_HP_step9$R2.top)
# FA1      FA2    R2
# 1     C20.5n3 C22.4n6c 98.23
# 2       C18.0 C22.4n6c 98.23
# 3    C18.1n9c C22.4n6c 98.23
# 4  C20.1n9.7c C22.4n6c 98.23
# 5    C20.4n6c C22.4n6c 98.19 <- ARA and Adrenic, elongation product
# 6    C18.1n9t C22.4n6c 98.18
# 7       C14.0 C22.4n6c 98.18
# 8    C16.1n7c C22.4n6c 98.16
# 9     C18.4n3 C22.4n6c 98.16
# 10   C22.4n6c    OBCFA 98.16
FVnMM_HP_LR9 = cbind(FVnMM_HP_LR8, FVnMM_HP_step9$logratios.top[,5])
FVnMM_HP_NR9 = rbind(FVnMM_HP_NR8, FVnMM_HP_step9$ratios.top[5,])
FVnMM_HP_R29 = c(FVnMM_HP_R28, FVnMM_HP_step9$R2.top[5])

# > > > > Step 10 -----
FVnMM_HP_step10 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR9)
print.ratios(FVnMM_HP_step10$names.top, FVnMM_HP_step10$R2.top)
# FA1      FA2    R2
# 1  C18.2n6c C20.4n6c 98.81
# 2  C18.2n6c C22.4n6c 98.81
# 3  C16.1n7c C18.2n6c 98.80 <- diatoms vs. greens
# 4     C16.0 C18.2n6c 98.80
# 5  C18.2n6c  C22.6n3 98.80
# 6  C18.2n6c  C18.4n3 98.80
# 7  C18.2n6c C22.1n9c 98.80
# 8  C18.2n6c    OBCFA 98.80
# 9  C18.2n6c  C22.5n3 98.80
# 10    C14.0 C18.2n6c 98.80
FVnMM_HP_LR10 = cbind(FVnMM_HP_LR9, FVnMM_HP_step10$logratios.top[,3])
FVnMM_HP_NR10 = rbind(FVnMM_HP_NR9, FVnMM_HP_step10$ratios.top[3,])
FVnMM_HP_R210 = c(FVnMM_HP_R29, FVnMM_HP_step10$R2.top[3])

# > > > > Step 11 -----
FVnMM_HP_step11 = STEP(data = ma_FVnMM_HPw, nsteps = 1, top = 10, previous = FVnMM_HP_LR10)
print.ratios(FVnMM_HP_step11$names.top, FVnMM_HP_step11$R2.top)
# FA1        FA2    R2
# 1   C20.5n3   C22.1n9c 99.22
# 2  C18.1n9c   C18.2n6c 99.22
# 3     C18.0      OBCFA 99.22
# 4     C18.0   C22.1n9c 99.22
# 5  C16.1n7c   C18.1n9c 99.22
# 6   C18.4n3 C20.1n9.7c 99.22
# 7   C20.5n3      OBCFA 99.22
# 8  C18.1n9c    C18.4n3 99.22
# 9  C18.1n9c    C22.6n3 99.22
# 10 C18.2n6c C20.1n9.7c 99.22

# Stop here. Very little additional variation accounted.

# > > > 3.2.2.3 : Organize step 10 LRs ----
rownames(FVnMM_HP_NR10) = paste("Step", 1:10, sep="")
colnames(FVnMM_HP_NR10) = c("FA1","FA2")
FVnMM_HP_FR = as.data.frame(cbind(FVnMM_HP_NR10,
                                  Ratio = paste(colnames(ma_FVnMM_HPw)[FVnMM_HP_NR10[,1]],
                                                "/",
                                                colnames(ma_FVnMM_HPw)[FVnMM_HP_NR10[,2]],
                                                sep="")))
FVnMM_HP_FR$FA1_lab = colnames(ma_FVnMM_HPw)[FVnMM_HP_NR10[,1]]
FVnMM_HP_FR$FA2_lab = colnames(ma_FVnMM_HPw)[FVnMM_HP_NR10[,2]]
FVnMM_HP_FR$step = 1:10

FVnMM_HP_FR$R2_cum = FVnMM_HP_R210
FVnMM_HP_FR = FVnMM_HP_FR |>
  arrange(step) |>
  mutate(R2 = R2_cum - lag(R2_cum,
                           default = 0))

FVnMM_HP_FR

# FVnMM_HP_FR_export = FVnMM_HP_FR
# FVnMM_HP_FR_export$R2_cum = round(FVnMM_HP_FR_export$R2_cum, 3)
# FVnMM_HP_FR_export$R2 = round(FVnMM_HP_FR_export$R2, 3)
# FVnMM_HP_FR_export
# write.csv(FVnMM_HP_FR_export, "CrabDTS_Table05_partial.csv")

colnames(FVnMM_HP_LR10) = FVnMM_HP_FR[,3]

# FVnMM_HP_PiR = sort(unique(as.numeric(FVnMM_HP_NR10))) # "parts in ratios"
FVnMM_HP_PiR = unique(as.numeric(t(FVnMM_HP_NR10))) # "parts in ratios"
colnames(ma_FVnMM_HPw)[FVnMM_HP_PiR]

# > > > 3.2.2.4 : Prep acyclic graph ----
FVnMM_HP_PiR_dim = data.frame(FA = colnames(ma_FVnMM_HPw)[FVnMM_HP_PiR],
                              dim1 = c(2, # C18.1n9c
                                       3, # C20.5n3
                                       1, # C18.0
                                       4, # C20.1n9.7c
                                       2, # C18.4n3
                                       1, # C22.1n9c
                                       3, # C22.6n3
                                       4, # OBCFA
                                       2, # C16.0
                                       3, # C16.1n7c
                                       1, # C20.4n6c
                                       1, # C22.4n6c
                                       4),# C18.2n6c
                              dim2 = c(4, # C18.1n9c
                                       4, # C20.5n3
                                       4, # C18.0
                                       4, # C20.1n9.7c
                                       3, # C18.4n3
                                       3, # C22.1n9c
                                       3, # C22.6n3
                                       3, # OBCFA
                                       2, # C16.0
                                       2, # C16.1n7c
                                       2, # C20.4n6c
                                       1, # C22.4n6c
                                       2)) # C18.2n6c

FVnMM_HP_PiR_dim$FA_labs = str_replace_all(FVnMM_HP_PiR_dim$FA,
                                           FA_labeller)

FVnMM_HP_FR = left_join(FVnMM_HP_FR,
                        FVnMM_HP_PiR_dim,
                        by = c("FA1_lab" = "FA"))
FVnMM_HP_FR = left_join(FVnMM_HP_FR,
                        FVnMM_HP_PiR_dim,
                        by = c("FA2_lab" = "FA"))

# > > > 3.2.2.5 : Prep PCA ----
FVnMM_HP_PCA = PCA(FVnMM_HP_LR10, weight = FALSE)

PLOT.PCA(FVnMM_HP_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 0.25)
points(FVnMM_HP_PCA$rowpcoord)

# Extract and calculate contribution vectors
axes_inv = c(1, 1)
dim = c(1, 2)
FVnMM_HP_RPC = FVnMM_HP_PCA$rowcoord[, dim] %*% diag(FVnMM_HP_PCA$sv[dim] * axes_inv)
FVnMM_HP_CSC = FVnMM_HP_PCA$colcoord[, dim] %*% diag(axes_inv)
FVnMM_HP_CPC = FVnMM_HP_CSC %*% diag(FVnMM_HP_PCA$sv[dim])
FVnMM_HP_CCC = FVnMM_HP_CSC * sqrt(FVnMM_HP_PCA$colmass)
FVnMM_HP_CRD = FVnMM_HP_CCC

# rescale contribution vectors for visibility
vectorscale = 0.33
FVnMM_HP_CRD_scale = FVnMM_HP_CRD * vectorscale

# Extract points
df_FVnMM_HP_RPC = cbind(df_FVnMM_HPw[, 1:4], FVnMM_HP_RPC)
colnames(df_FVnMM_HP_RPC) = c(colnames(df_FVnMM_HPw[, 1:4]),'dim1','dim2')

# create SexDev
df_FVnMM_HP_RPC$SexDev = paste(df_FVnMM_HP_RPC$Sex,
                               df_FVnMM_HP_RPC$Dev,
                               sep = "_")

df_FVnMM_HP_RPC$SexDev = factor(df_FVnMM_HP_RPC$SexDev,
                                levels = c("F_imm",
                                           "F_mat",
                                           "M_imm",
                                           "M_mat"))

# Aggregate by type, sex, and maturity
df_FVnMM_HP_RPC_MSD = df_FVnMM_HP_RPC |>
  group_by(Type, SexDev) |>
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SD = sd(dim1),
            dim2SD = sd(dim2))

# Calculate variance explained by PCA dimensions
FVnMM_HP_PCA_perc_1 = 100 * FVnMM_HP_PCA$sv[dim[1]]^2 / sum(FVnMM_HP_PCA$sv^2)
FVnMM_HP_PCA_perc_2 = 100 * FVnMM_HP_PCA$sv[dim[2]]^2 / sum(FVnMM_HP_PCA$sv^2)

df_FVnMM_HP_RPC$ID2 = paste(df_FVnMM_HP_RPC$Type,
                            df_FVnMM_HP_RPC$Sex,
                            df_FVnMM_HP_RPC$Mat,
                            df_FVnMM_HP_RPC$Tag_no_Crab_ID,
                            sep = ".")

# Format arrows
df_FVnMM_HP_arrows = data.frame(dim1 = FVnMM_HP_CRD_scale[ ,1],
                                dim2 = FVnMM_HP_CRD_scale[, 2],
                                LR = FVnMM_HP_PCA$colnames)

df_FVnMM_HP_arrows$step = 1:10

df_FVnMM_HP_arrows$angle = 90 - ((atan2(df_FVnMM_HP_arrows$dim1,
                                        df_FVnMM_HP_arrows$dim2) *
                                    180) /
                                   pi)
df_FVnMM_HP_arrows$angle = ifelse(df_FVnMM_HP_arrows$angle > 90 & df_FVnMM_HP_arrows$angle < 270,
                                  df_FVnMM_HP_arrows$angle + 180,
                                  df_FVnMM_HP_arrows$angle)
df_FVnMM_HP_arrows$hjust = ifelse(df_FVnMM_HP_arrows$angle > 90,
                                  1,
                                  0)

# Format arrow labels
df_FVnMM_HP_arrows$labs = df_FVnMM_HP_arrows$LR
df_FVnMM_HP_arrows$labs = str_replace_all(df_FVnMM_HP_arrows$labs,
                                          FA_labeller)
df_FVnMM_HP_arrows$labs = str_replace_all(df_FVnMM_HP_arrows$labs,
                                          "/",
                                          " / ")

# establish levels for IDs to keep colors consistent
df_FVnMM_HP_RPC$Tag_no_Crab_ID = factor(df_FVnMM_HP_RPC$Tag_no_Crab_ID,
                                        levels = paste0("FV", sprintf("%02.0f", 1:20)))

# Create hulls
df_FVnMM_HP_RPC_hull = df_FVnMM_HP_RPC |>
  group_by(Type, SexDev) |>
  slice(chull(dim1, dim2))

df_FVnMM_HP_RPC_hull2 = df_FVnMM_HP_RPC |>
  group_by(Tag_no_Crab_ID) |>
  slice(chull(dim1, dim2))

# > > > 3.2.2.6 : Plot PCA ----
# > > > > 3.2.2.6.1 : By Type, SexDev ----
plot_FVnMM_HP_PCA = ggplot() +
  # geom_polygon(aes(x = dim1,
  #                  y = dim2,
  #                  color = SexDev,
  #                  linetype = Type,
  #                  group = interaction(SexDev, Type)),
  #              linewidth = 1,
  #              fill = NA,
  #              data = df_FVnMM_HP_RPC_hull) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = SexDev,
                 shape = Type),
             alpha = 0.5,
             size = 2,
             data = df_FVnMM_HP_RPC) +
  # geom_segment(aes(x = 0,
  #                  y = 0,
  #                  xend = FVnMM_HP_CRD_scale[ ,1],
  #                  yend = FVnMM_HP_CRD_scale[ ,2]),
  #              color = "grey50",
  #              arrow = arrow(length = unit(0.2, "cm"))) +
  geom_errorbar(aes(x = dim1M,
                    ymin = dim2M - dim2SD,
                    ymax = dim2M + dim2SD,
                    color = SexDev),
                width = 0,
                linewidth = 0.5,
                data = df_FVnMM_HP_RPC_MSD) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - dim1SD,
                     xmax = dim1M + dim1SD,
                     color = SexDev),
                 height = 0,
                 linewidth = 0.5,
                 data = df_FVnMM_HP_RPC_MSD) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = SexDev,
                 shape = Type),
             size = 3,
             data = df_FVnMM_HP_RPC_MSD) +
  # geom_point(aes(x = ifelse(dim1 > 0,
  #                           dim1 + 0.075 * cos(angle * pi/180),
  #                           dim1 - 0.075 * cos(angle * pi/180)),
  #                y = ifelse(dim1 > 0,
  #                           dim2 + 0.075 * sin(angle * pi/180),
  #                           dim2 - 0.075 * sin(angle * pi/180))),
  #            colour = "grey50",
  #            size = 4,
  #            data = df_FVnMM_HP_arrows)+
  # geom_text(aes(x = ifelse(dim1 > 0,
  #                          dim1 + 0.075 * cos(angle * pi/180),
#                          dim1 - 0.075 * cos(angle * pi/180)),
#               y = ifelse(dim1 > 0,
#                          dim2 + 0.075 * sin(angle * pi/180),
#                          dim2 - 0.075 * sin(angle * pi/180)),
#               label = step),
#           colour = "white",
#           size = unit(2, "pt"),
#           data = df_FVnMM_HP_arrows)+
scale_color_viridis_d(option = "magma",
                      begin = 0.2,
                      end = 0.9,
                      direction = -1,
                      drop = TRUE,
                      labels = c("Female, imm.",
                                 "Female, mat.",
                                 "Male, imm.",
                                 "Male, mat.\n(omitted)")) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("DTS", "FRZ"),
                        name = "Type") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +  
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = NULL,
       # x = paste("PC 1 (",
       #           round(FVnMM_HP_PCA_perc_1, 1),
       #           "%)",
       #           sep=""),
       y = paste("PC 2 (",
                 round(FVnMM_HP_PCA_perc_2, 1),
                 "%)",
                 sep=""),
       color = "Sex, Maturity") +
  coord_fixed(
    xlim = c(-0.5, 0.4),
    ylim = c(-0.25, 0.375)
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        legend.box = "vertical",
        axis.text.y.right = element_blank(),
        # axis.text.y.right = element_text(colour="grey50",
        #                                  size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50")) +
  guides(color = guide_legend(override.aes = list(size = 2),
                              nrow = 1),
         shape = guide_legend(override.aes = list(size = 2)))

plot_FVnMM_HP_PCA

# > > > > 3.2.2.6.2 : By Individual ----
plot_FVnMM_HP_PCA2 = ggplot() +
  geom_polygon(aes(x = dim1,
                   y = dim2,
                   color = Tag_no_Crab_ID,
                   group = Tag_no_Crab_ID),
               linewidth = 1,
               fill = NA,
               data = df_FVnMM_HP_RPC_hull2) +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = Tag_no_Crab_ID,
                 shape = Type),
             alpha = 0.5,
             size = 2,
             data = df_FVnMM_HP_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = Tag_no_Crab_ID,
  #               label = Tag_no_Crab_ID),
  #           alpha = 0.5,
  #           size = 3,
  #           data = df_FVnMM_HP_RPC[df_FVnMM_HP_RPC$Type == "Tissue",]) +
  # geom_point(aes(x = ifelse(dim1 > 0,
  #                           dim1 + 0.075 * cos(angle * pi/180),
  #                           dim1 - 0.075 * cos(angle * pi/180)),
  #                y = ifelse(dim1 > 0,
#                           dim2 + 0.075 * sin(angle * pi/180),
#                           dim2 - 0.075 * sin(angle * pi/180))),
#            colour = "grey50",
#            size = 4,
#            data = df_FVnMM_HP_arrows)+
# geom_text(aes(x = ifelse(dim1 > 0,
#                          dim1 + 0.075 * cos(angle * pi/180),
#                          dim1 - 0.075 * cos(angle * pi/180)),
#               y = ifelse(dim1 > 0,
#                          dim2 + 0.075 * sin(angle * pi/180),
#                          dim2 - 0.075 * sin(angle * pi/180)),
#               label = step),
#           colour = "white",
#           size = unit(2, "pt"),
#           data = df_FVnMM_HP_arrows)+
# geom_segment(aes(x = 0,
#                  y = 0,
#                  xend = FVnMM_HP_CRD_scale[ ,1],
#                  yend = FVnMM_HP_CRD_scale[ ,2]),
#              color = "grey50",
#              arrow = arrow(length = unit(0.2, "cm"))) +
scale_color_viridis_d(option = "mako",
                      begin = 0.2,
                      end = 0.9,
                      direction = -1,
                      guide = "none",
                      drop = FALSE) +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type",
                     guide = "none") +  
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FVnMM_HP_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       # y = paste("PC 2 (",
       #           round(FVnMM_HP_PCA_perc_2, 1),
       #           "%)",
       #           sep=""),
       y = NULL,
       color = "Crab ID") +
  coord_fixed(
    xlim = c(-0.5, 0.4),
    ylim = c(-0.25, 0.375)
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        axis.text.y = element_blank(),
        legend.box = "vertical",
        axis.text.y.right = element_blank(),
        # axis.text.y.right = element_text(colour="grey50",
        #                                  size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50")) 

plot_FVnMM_HP_PCA2

# > > > > 3.2.2.6.3 : Baseplot with Vectors ----
plot_FVnMM_HP_PCA_base = 
  ggplot() +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FVnMM_HP_CRD_scale[ ,1],
                   yend = FVnMM_HP_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(aes(x = ifelse(dim1 > 0,
                            dim1 + 0.075 * cos(angle * pi/180),
                            dim1 - 0.075 * cos(angle * pi/180)),
                 y = ifelse(dim1 > 0,
                            dim2 + 0.075 * sin(angle * pi/180),
                            dim2 - 0.075 * sin(angle * pi/180))),
             colour = "grey50",
             size = 3,
             data = df_FVnMM_HP_arrows)+
  geom_text(aes(x = ifelse(dim1 > 0,
                           dim1 + 0.075 * cos(angle * pi/180),
                           dim1 - 0.075 * cos(angle * pi/180)),
                y = ifelse(dim1 > 0,
                           dim2 + 0.075 * sin(angle * pi/180),
                           dim2 - 0.075 * sin(angle * pi/180)),
                label = step),
            colour = "white",
            size = unit(2, "pt"),
            data = df_FVnMM_HP_arrows) +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FVnMM_HP_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       # y = paste("PC 2 (",
       #           round(FVnMM_HP_PCA_perc_2, 1),
       #           "%)",
       #           sep=""),
       y = NULL,
       color = "Sex, Maturity") +
  coord_fixed(
    xlim = c(-0.5, 0.4),
    ylim = c(-0.25, 0.375)
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        axis.text.y = element_blank(),
        legend.box = "vertical",
        axis.text.y.right = element_text(colour="grey50",
                                         size = unit(6, "pt")),
        axis.text.x.top = element_text(colour="grey50",
                                       size = unit(6, "pt")),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50")) +
  guides(color = guide_legend(override.aes = list(size = 2),
                              nrow = 1),
         shape = guide_legend(override.aes = list(size = 2)))
plot_FVnMM_HP_PCA_base

# > > > 3.2.2.7 : Plot acyclic graph ----
plot_FVnMM_HP_LR_graph = ggplot() +
  geom_segment(aes(x = dim1.x,
                   y = dim2.x,
                   xend = dim1.y,
                   yend = dim2.y,
                   color = step),
               linewidth = 1,
               data = FVnMM_HP_FR) +
  geom_label(aes(x = dim1,
                 y = dim2,
                 label = FA_labs),
             size = unit(3, "pt"),
             data = FVnMM_HP_PiR_dim) +
  geom_point(aes(x = (dim1.x + dim1.y)/2,
                 y = (dim2.x + dim2.y)/2,
                 color = step),
             size = 6,
             data = FVnMM_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2,
                y = (dim2.x + dim2.y)/2,
                label = step),
            color = "white",
            size = unit(3, "pt"),
            data = FVnMM_HP_FR) +
  geom_text(aes(x = (dim1.x + dim1.y)/2 + 0.125,
                y = (dim2.x + dim2.y)/2 - 0.125,
                label = round(R2*100, 0),
                color = step),
            size = unit(2, "pt"),
            data = FVnMM_HP_FR) +
  scale_color_gradient(high = "dodgerblue4",
                       low = "dodgerblue1") +
  coord_cartesian(xlim = c(0.5, 4.5),
                  clip = "off") +
  # coord_fixed(xlim = c(0.25, 4.75),
  #             clip = "off") +
  theme_void() +
  theme(legend.position = "none")
plot_FVnMM_HP_LR_graph 

# > > > 3.2.2.8 : Export acyclic graph and PCA ----
# tiff("plot_C_opilio_FVnMM_HP_LR_doc.tiff",
#      height = res * 6.5,
#      width = res * 7.5,
#      res = res)
# plot_FVnMM_HP_LR_graph / (plot_FVnMM_HP_PCA | plot_FVnMM_HP_PCA2) +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(heights = c(1, 2))
# dev.off()
# 
# plot_FV_HP_LR_graph + plot_FV_HP_PCA + plot_FV_HP_PCA2 +
#   plot_FVnMM_HP_LR_graph + plot_FVnMM_HP_PCA + plot_FVnMM_HP_PCA2 +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(ncol = 3,
#               heights = c(1,1),
#               widths = c(1.5, 1, 1))
# 
# tiff("plot_Field_Validation.tiff",
#      width = 12 * res,
#      height = 9 * res,
#      res = res)
# plot_FV_HP_LR_graph + plot_FV_HP_PCA + plot_FV_HP_PCA2 +
#   plot_FVnMM_HP_LR_graph + plot_FVnMM_HP_PCA + plot_FVnMM_HP_PCA2 +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(ncol = 3,
#               heights = c(1,1),
#               widths = c(1.5, 1, 1))
# dev.off()

# > > > 3.2.2.9 : Univariate plot ----
df_FVnMM_HP_LR10w = df_FVnMM_HPw[,1:4]
df_FVnMM_HP_LR10w$SexDev = paste(df_FVnMM_HP_LR10w$Sex,
                                 df_FVnMM_HP_LR10w$Dev,
                                 sep = "_")

df_FVnMM_HP_LR10w$SexDev = factor(df_FVnMM_HP_LR10w$SexDev,
                                  levels = c("F_imm",
                                             "F_mat",
                                             "M_imm"))

df_FVnMM_HP_LR10w = cbind(df_FVnMM_HP_LR10w,
                          FVnMM_HP_LR10)

df_FVnMM_HP_LR10 = df_FVnMM_HP_LR10w |>
  pivot_longer(cols = 6:length(df_FVnMM_HP_LR10w),
               names_to = "FAs",
               values_to = "LR")

df_FVnMM_HP_LR10$FAs = str_replace_all(df_FVnMM_HP_LR10$FAs,
                                       FA_labeller)

df_FVnMM_HP_LR10$FAs = factor(df_FVnMM_HP_LR10$FAs,
                              levels = str_replace_all(FVnMM_HP_FR$Ratio,
                                                       FA_labeller))
df_FVnMM_HP_LR10_sum = df_FVnMM_HP_LR10 |>
  group_by(Type, SexDev, FAs) |>
  summarise(M = mean(LR),
            SD = sd(LR))

# > > > 3.2.2.9.1 : Plot ----
plot_FVnMM_HP_LR_cat = ggplot() +
  geom_hline(aes(yintercept = y),
             linetype = 2,
             color = "grey75",
             data = labels_logratios) +
  geom_errorbar(aes(x = SexDev,
                    color = Type,
                    ymin = M - SD,
                    ymax = M + SD),
                width = 0,
                alpha = 0.5,
                position = position_dodge(0.5),
                data = df_FVnMM_HP_LR10_sum) +
  geom_point(aes(x = SexDev,
                 color = Type,
                 y = M),
             size = 2,
             alpha = 0.5,
             shape = 9,
             position = position_dodge(0.5),
             data = df_FVnMM_HP_LR10_sum) +
  geom_point(aes(x = SexDev,
                 y = LR,
                 color = Type),
             size = 0.5,
             position = position_dodge(0.5),
             data = df_FVnMM_HP_LR10) +
  scale_color_manual(values = c("firebrick",
                                "dodgerblue")) +
  scale_x_discrete(labels = c("Female, imm.",
                              "Female, mat.",
                              "Male, imm.",
                              "Male, mat.")) +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),
                                         breaks = exp(labels_logratios$y),
                                         labels = labels_logratios$label,
                                         name = "Ratio"),
                     limits = c(floor(min(df_FVnMM_HP_LR10$LR) / 0.5) * 0.5,
                                ceiling(max(df_FVnMM_HP_LR10$LR) / 0.5) * 0.5)) +
  labs(x = "Sex, Maturity",
       y = "Logratio",
       color = "Sample Type") +
  theme_classic() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text = element_text(size = unit(4, "pt")),
        axis.title = element_text(size = unit(6, "pt")),
        axis.text.y.right = element_text(color="grey70"),
        axis.ticks.y.right = element_line(color="grey70"),
        axis.line.y.right = element_line(color="grey70"),
        axis.title.y.right = element_text(color = "grey70"),
        legend.background = element_blank(),
        legend.text = element_text(size = unit(4, "pt")),
        legend.title = element_text(size = unit(6, "pt")),
        legend.spacing = unit(0, "pt"),
        strip.text = element_text(size = unit(4, "pt"))) +
  facet_grid(.~ FAs)
plot_FVnMM_HP_LR_cat

# > > > 3.2.2.9.2 : Export ----
# tiff("plot_C_opilio_FVnMM_HP_LR_cat.tiff",
#      height = res * 4,
#      width = res * 7.5,
#      res = res)
# plot_FVnMM_HP_LR_cat
# dev.off()

# > > > 3.2.2.10: PERMANOVA on LRs -----

# Permutation structure
FVnMM_LR10_perm = how(within = Within(type = "free"), # samples are free to shuffle within individual
                      plots = 
                        Plots(strata = df_FVnMM_HP_LR10w$Tag_no_Crab_ID, 
                              type = "free"), # samples cannot shuffle across individual
                      nperm = 9999)

# Check number of possible permutations
numPerms(df_FVnMM_HP_LR10w,
         FVnMM_LR10_perm)

# Run PERMANOVA
# We include these terms to allow proper estimation of pseudo-F statistics,
# which rely on proper denominators relevant to each term. Using the residual
# as done with adonis2() default is note (always) appropriate (Anderson and Ter Braak
# 2003 Journal of Statistical Computation and Simulation). Should be based on
# exchangeable units. The method used below is from Bakker 2024, Applied 
# Multivariate Statistics in R.

# F_Type = MS_Type / MS_Residual
# F_SexDev = MS_SexDev / MS_Tag_no_Crab_ID

# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a 
# nesting factor/random-ish effect

# We run it once to check our general structure. We do not use this output for
# inference because we know the pseudo-Fs are not calculated how we want.
FVnMM_LR10_PERMANOVA_base = adonis2(FVnMM_HP_LR10 ~ 
                                      Type +
                                      SexDev + # SexDev must come first, as its variation would be otherwise subsumed by Crab ID
                                      Tag_no_Crab_ID,
                                    data = df_FVnMM_HP_LR10w,
                                    method = "euclidean",
                                    permutations = FVnMM_LR10_perm,
                                    by = "terms")
FVnMM_LR10_PERMANOVA_base

set.seed(207)

FVnMM_LR10_perms = rbind(1:nrow(df_FVnMM_HP_LR10w),
                         shuffleSet(n = nrow(df_FVnMM_HP_LR10w),
                                    control = FVnMM_LR10_perm,
                                    nset = 9999))

FVnMM_LR10_PERMANOVA = matrix(nrow = nrow(FVnMM_LR10_perms),
                              ncol = 5)

colnames(FVnMM_LR10_PERMANOVA) = c("Type",
                                   "SexDev",
                                   "Tag_no_Crab_ID",
                                   "Residual",
                                   "Total")

for (i in 1:nrow(FVnMM_LR10_perms)) {
  df_temp = df_FVnMM_HP_LR10w[FVnMM_LR10_perms[i, ], ]
  res_temp = adonis2(FVnMM_HP_LR10 ~ 
                       Type +
                       SexDev + # SexDev must come first, as its variation would be otherwise subsumed by Crab ID
                       Tag_no_Crab_ID,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  FVnMM_LR10_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-Fs
FVnMM_LR10_PERMANOVA = 
  FVnMM_LR10_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 14),
         F_SexDev = (SexDev / 2) / (Tag_no_Crab_ID / 12))

# Calculate p-values
FVnMM_LR10_PERMANOVA$F_Type[1] # pseudo-F of data
with(FVnMM_LR10_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

FVnMM_LR10_PERMANOVA$F_SexDev[1] # pseudo-F of data
with(FVnMM_LR10_PERMANOVA, sum(F_SexDev >= F_SexDev[1]) / length(F_SexDev))

set.seed(NULL)

# > > > 3.2.2.11: nMDS and PERMANOVA on raw percentages -----
# Do not include replacements for 0 values.

# > > > > 3.2.2.11.1: nMDS -----
# > > > > > 3.2.2.11.1.1: Fit NMDS ----
FVnMM_nMDS = metaMDS(ma_FVnMM_HPw,
                     k = 2,
                     distance = "euclidean",
                     wascores = F,
                     noshare = F,
                     autotransform = F)
NMDS_stress = FVnMM_nMDS$stress

# View Stressplot
stressplot(FVnMM_nMDS)

# Extract scores from fit
df_FVnMM_nMDS_plot = scores(FVnMM_nMDS, display = c("sites")) # pull scores from MDS fit 
df_FVnMM_nMDS_plot = data.frame(df_FVnMM_HPw[,1:4], df_FVnMM_nMDS_plot)

# get arrows of FAs correlating to first two dimensions
FVnMM_nMDS_envfit = envfit(FVnMM_nMDS, ma_FVnMM_HPw, perm = 99999) 
# fit vector arrows of FAs correlating with dimensions

# pull scores vector arrows of FAs correlating with dimensions
FA_arrows = scores(FVnMM_nMDS_envfit, display = c("vectors")) 

# calculate scaling multiplier based on min range of points in NMDS
scalefactor = 0.35 * min(max(df_FVnMM_nMDS_plot$NMDS1) - min(df_FVnMM_nMDS_plot$NMDS1), 
                         max(df_FVnMM_nMDS_plot$NMDS2) - min(df_FVnMM_nMDS_plot$NMDS2))
# makes scaling multiplier based on min range of points

# scaling of arrows to fit plot area (with vegan's ordiArrowMul function)
#FA_arrows=FA_arrows*ordiArrowMul(FVnMM_nMDS_plot) 

# scaling of arrows to fit plot area based on scalefactor above
FA_arrows = FA_arrows * scalefactor 

# save those as a dataframe
FA_arrows = as.data.frame(FA_arrows) 

# Use onlyFA with over 5% representation for arrows:
FA_arrows = FA_arrows[rownames(FA_arrows) %in% df_FVnMM_HP_list_5,]

# FA arrow text formatting
arrow_angle = 90 - ((atan2(FA_arrows$NMDS1, FA_arrows$NMDS2) * 180) / pi)
arrow_angle = ifelse(arrow_angle > 90 & arrow_angle < 270,
                     arrow_angle + 180,
                     arrow_angle)
arrow_hjust = ifelse(arrow_angle > 90,
                     1,
                     0)
FA_arrows$labs = as.character(rownames(FA_arrows))
FA_arrows$labs = gsub("C", "", FA_arrows$labs)
FA_arrows$labs = gsub("B", "", FA_arrows$labs)
FA_arrows$labs = sub("c", "*c", FA_arrows$labs)
FA_arrows$labs = sub("t", "*t", FA_arrows$labs)
FA_arrows$labs = sub("\\.", ":", FA_arrows$labs)
FA_arrows$labs = sub("n", "*omega*-", FA_arrows$labs)

# Stress annotation
fontfamily = "sans"

stressanot = grobTree(textGrob(paste("stress = ",
                                     round(NMDS_stress,3),
                                     sep = ''),
                               x = 0.03,
                               y = 0.975,
                               hjust = 0,
                               gp = gpar(fontsize = 6,
                                         fontfamily = fontfamily)))

# create SexDev
df_FVnMM_nMDS_plot$SexDev = paste(df_FVnMM_nMDS_plot$Sex,
                                  df_FVnMM_nMDS_plot$Dev,
                                  sep = "_")

df_FVnMM_nMDS_plot$SexDev = factor(df_FVnMM_nMDS_plot$SexDev,
                                   levels = c("F_imm",
                                              "F_mat",
                                              "M_imm",
                                              "M_mat"))

# Aggregate by type and storage time
df_FVnMM_nMDS_MSD = df_FVnMM_nMDS_plot |>
  group_by(Type, SexDev) |>
  summarize(NMDS1M = mean(NMDS1),
            NMDS2M = mean(NMDS2),
            NMDS1SD = sd(NMDS1),
            NMDS2SD = sd(NMDS2))

# Create hulls
df_FVnMM_nMDS_hull = df_FVnMM_nMDS_plot |>
  group_by(Type, SexDev) |>
  slice(chull(NMDS1, NMDS2))

df_FVnMM_nMDS_hull2 = df_FVnMM_nMDS_plot |>
  group_by(Tag_no_Crab_ID) |>
  slice(chull(NMDS1, NMDS2))

# > > > > > 3.2.2.11.1.1: Plot NMDS ----
# > > > > > > 3.2.2.11.1.1.1: By Type, SexDev ----
plot_FVnMM_HP_nMDS = ggplot()+
  # geom_polygon(aes(x = NMDS1,
  #                  y = NMDS2,
  #                  color = SexDev,
  #                  linetype = Type,
  #                  group = interaction(SexDev, Type)),
  #              linewidth = 1,
  #              fill = NA,
  #              data = df_FVnMM_nMDS_hull) +
  geom_point(data = df_FVnMM_nMDS_plot,
             aes(NMDS1,
                 NMDS2,
                 colour = SexDev,
                 shape = Type),
             size = 2,
             alpha = 0.5)+
  # geom_segment(data = FA_arrows,
  #              aes(x = 0,
  #                  y = 0,
  #                  xend = NMDS1,
  #                  yend = NMDS2),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              colour = "gray50",
  #              alpha = 1)+
  # geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
  #           aes(NMDS1, NMDS2, label =  FA_arrows$labs),
  #           color = "grey50",
#           size = unit(2, "pt"),
#           parse = T,
#           angle = arrow_angle,
#           hjust = arrow_hjust) +
geom_errorbar(aes(x = NMDS1M,
                  ymin = NMDS2M - NMDS2SD,
                  ymax = NMDS2M + NMDS2SD,
                  color = SexDev),
              width = 0,
              linewidth = 0.5,
              data = df_FVnMM_nMDS_MSD) +
  geom_errorbarh(aes(y = NMDS2M,
                     xmin = NMDS1M - NMDS1SD,
                     xmax = NMDS1M + NMDS1SD,
                     color = SexDev),
                 height = 0,
                 linewidth = 0.5,
                 data = df_FVnMM_nMDS_MSD) +
  geom_point(aes(x = NMDS1M,
                 y = NMDS2M,
                 color = SexDev,
                 shape = Type),
             size = 3,
             data = df_FVnMM_nMDS_MSD) +
  scale_color_viridis_d(option = "magma",
                        begin = 0.2,
                        end = 0.9,
                        direction = -1,
                        drop = TRUE,
                        labels = c("Female, imm.",
                                   "Female, mat.",
                                   "Male, imm.",
                                   "Male, mat.\n(omitted)")) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("DTS", "FRZ"),
                        name = "Type") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type") +  
  theme_classic() +
  coord_fixed(xlim = c(-5.5, 3.75),
              ylim = c(-3, 4.5)) +
  labs(x = NULL,
       # x = "nMDS dim 1",
       y = "nMDS dim 2",
       color = "Sex, Maturity") +
  annotation_custom(stressanot) +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size = 2),
                              nrow = 1),
         shape = guide_legend(override.aes = list(size = 2)))

plot_FVnMM_HP_nMDS

# > > > > > > 3.2.2.11.1.1.2: By Individual ----
plot_FVnMM_HP_nMDS2 = ggplot()+
  geom_polygon(aes(x = NMDS1,
                   y = NMDS2,
                   color = Tag_no_Crab_ID),
               linewidth = 1,
               fill = NA,
               data = df_FVnMM_nMDS_hull2) +
  geom_point(data = df_FVnMM_nMDS_plot,
             aes(NMDS1,
                 NMDS2,
                 colour = Tag_no_Crab_ID,
                 shape = Type),
             size = 2,
             alpha = 0.5)+
  # geom_segment(data = FA_arrows,
  #              aes(x = 0,
  #                  y = 0,
  #                  xend = NMDS1,
  #                  yend = NMDS2),
  #              arrow = arrow(length = unit(0.2, "cm")),
  #              colour = "gray50",
  #              alpha = 1)+
  # geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
  #           aes(NMDS1, NMDS2, label =  FA_arrows$labs),
  #           color = "grey50",
#           size = unit(2, "pt"),
#           parse = T,
#           angle = arrow_angle,
#           hjust = arrow_hjust) +
scale_color_viridis_d(option = "mako",
                      begin = 0.2,
                      end = 0.9,
                      direction = -1,
                      guide = "none",
                      drop = FALSE) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("DTS", "FRZ"),
                        name = "Type") +
  scale_shape_manual(values = c(17, 16),
                     labels = c("DTS", "FRZ"),
                     name = "Type",
                     guide = "none") +  
  theme_classic() +
  coord_fixed(xlim = c(-5.5, 3.75),
              ylim = c(-3, 4.5)) +
  labs(x = "nMDS dim 1",
       y = NULL,
       # y = "nMDS dim 2",
       color = "Crab ID") +
  annotation_custom(stressanot)+
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical")

plot_FVnMM_HP_nMDS2

# > > > > > > 3.2.2.11.1.1.3: Baseplot with Vectors ----
plot_FVnMM_HP_nMDS_base = 
  ggplot() +
  geom_segment(data = FA_arrows,
               aes(x = 0,
                   y = 0,
                   xend = NMDS1,
                   yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "gray50",
               alpha = 1)+
  geom_text(data = as.data.frame(FA_arrows[,c("NMDS1","NMDS2")] * 1.1),
            aes(NMDS1, NMDS2, label =  FA_arrows$labs),
            color = "grey50",
            size = unit(2, "pt"),
            parse = T,
            angle = arrow_angle,
            hjust = arrow_hjust) +
  theme_classic() +
  coord_fixed(xlim = c(-5.5, 3.75),
              ylim = c(-3, 4.5)) +
  labs(x = NULL,
       y = NULL,
       # x = "nMDS dim 1",
       # y = "nMDS dim 2",
       color = "Sex, Maturity") +
  annotation_custom(stressanot) +
  theme(legend.position = "bottom",
        legend.spacing = unit(0, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = unit(6, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size = 2),
                              nrow = 1),
         shape = guide_legend(override.aes = list(size = 2)))

plot_FVnMM_HP_nMDS_base

# > > > > > > 3.2.2.11.1.1.4: Combine and Export ----
plot_FVnMM_HP_PCA + 
  plot_FVnMM_HP_PCA2 +
  plot_FVnMM_HP_PCA_base +
  plot_FVnMM_HP_nMDS +
  plot_FVnMM_HP_nMDS2 +
  plot_FVnMM_HP_nMDS_base +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(20, "pt"))

ggsave("CrabDTS_Fig04.tiff",
       device = "tiff",
       height = 4.75,
       width = 7.5,
       dpi = res)

# > > > 3.1.11.2: PERMANOVA -----

# We use the same contraints on permutations to account for individual as in the
# PERMANOVA using the selected logratios
# create SexDev
df_FVnMM_HPw$SexDev = paste(df_FVnMM_HPw$Sex,
                            df_FVnMM_HPw$Dev,
                            sep = "_")

df_FVnMM_HPw$SexDev = factor(df_FVnMM_HPw$SexDev,
                             levels = c("F_imm",
                                        "F_mat",
                                        "M_imm",
                                        "M_mat"))

# Permutation structure for only male crab
FVnMM_perc_perm = how(within = Within(type = "free"), # samples are free to shuffle within individual
                      plots = 
                        Plots(strata = df_FVnMM_HPw$Tag_no_Crab_ID, 
                              type = "free"), # samples cannot shuffle across individual
                      nperm = 9999)

# Check number of possible permutations
numPerms(df_FVnMM_HPw,
         FVnMM_perc_perm)

# Run PERMANOVA
# We include these terms to allow proper estimation of pseudo-F statistics,
# which rely on proper denominators relevant to each term. Using the residual
# as done with adonis2() default is note (always) appropriate (Anderson and Ter Braak
# 2003 Journal of Statistical Computation and Simulation). Should be based on
# exchangeable units. The method used below is from Bakker 2024, Applied 
# Multivariate Statistics in R.

# F_Type = MS_Type / MS_Residual
# F_SexDev = MS_SexDev / MS_Tag_no_Crab_ID
# The pseudo-F of Tag_no_Crab_ID is not evaluated as it is essentially a 
# nesting factor/random-ish effect

# We run it once to check our general structure. We do not use this output for
# inference because we know the pseudo-Fs are not calculated how we want.
FVnMM_perc_PERMANOVA_base = adonis2(ma_FVnMM_HPw ~ 
                                      Type +
                                      SexDev + # SexDev must come first, as its variation would be otherwise subsumed by Crab ID
                                      Tag_no_Crab_ID,
                                    data = df_FVnMM_HPw,
                                    method = "euclidean",
                                    permutations = FVnMM_perc_perm,
                                    by = "terms")
FVnMM_perc_PERMANOVA_base

set.seed(207)

FVnMM_perc_perms = rbind(1:nrow(df_FVnMM_HPw),
                         shuffleSet(n = nrow(df_FVnMM_HPw),
                                    control = FVnMM_perc_perm,
                                    nset = 9999))

FVnMM_perc_PERMANOVA = matrix(nrow = nrow(FVnMM_perc_perms),
                              ncol = 5)

colnames(FVnMM_perc_PERMANOVA) = c("Type",
                                   "SexDev",
                                   "Tag_no_Crab_ID",
                                   "Residual",
                                   "Total")

for (i in 1:nrow(FVnMM_perc_perms)) {
  df_temp = df_FVnMM_HPw[FVnMM_perc_perms[i, ], ]
  res_temp = adonis2(ma_FVnMM_HPw ~ 
                       Type +
                       SexDev + # SexDev must come first, as its variation would be otherwise subsumed by Crab ID
                       Tag_no_Crab_ID,
                     data = df_temp,
                     by = "terms",
                     method = "euclidean",
                     permutations = 0)
  FVnMM_perc_PERMANOVA[i, ] = t(res_temp$SumOfSqs)
}

# Calculate pseudo-Fs
FVnMM_perc_PERMANOVA = 
  FVnMM_perc_PERMANOVA |>
  data.frame() |>
  mutate(F_Type = (Type / 1) / (Residual / 14),
         F_SexDev = (SexDev / 2) / (Tag_no_Crab_ID / 12))

# Calculate p-values
FVnMM_perc_PERMANOVA$F_Type[1] # pseudo-F of data
with(FVnMM_perc_PERMANOVA, sum(F_Type >= F_Type[1]) / length(F_Type))

FVnMM_perc_PERMANOVA$F_SexDev[1] # pseudo-F of data
with(FVnMM_perc_PERMANOVA, sum(F_SexDev >= F_SexDev[1]) / length(F_SexDev))

set.seed(NULL)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####

# 4: Copeman et al. 2021 Polar Biology -----
df_Copeman =
  df_Copeman |>
  mutate(LR_16.0_16.1n7 = log(C16.0 / C16.1n7))

ggplot() +
  geom_hline(aes(yintercept = y),
             linetype = 2,
             color = "grey75",
             data = labels_logratios) +
  geom_hline(aes(yintercept = LR_16.0_16.1n7,
                 linetype = factor(Year)),
             color = "goldenrod2",
             data = df_Copeman) +
  geom_label_repel(aes(x = -0.5,
                       y = LR_16.0_16.1n7,
                       label = Region),
                   color = "goldenrod2",
                   data = df_Copeman) +
  geom_errorbar(aes(x = SexDev,
                    color = Type,
                    ymin = M - SD,
                    ymax = M + SD),
                width = 0,
                alpha = 1,
                position = position_dodge(0.5),
                data = df_FVnMM_HP_LR10_sum[df_FVnMM_HP_LR10_sum$FAs == "Palmitic/Palmitoleic",]) +
  geom_point(aes(x = SexDev,
                 color = Type,
                 y = M),
             size = 3,
             alpha = 1,
             shape = 9,
             position = position_dodge(0.5),
             data = df_FVnMM_HP_LR10_sum[df_FVnMM_HP_LR10_sum$FAs == "Palmitic/Palmitoleic",]) +
  geom_point(aes(x = SexDev,
                 y = LR,
                 color = Type),
             size = 1,
             position = position_dodge(0.5),
             data = df_FVnMM_HP_LR10[df_FVnMM_HP_LR10$FAs == "Palmitic/Palmitoleic",]) +
  scale_color_manual(values = c("firebrick",
                                "dodgerblue")) +
  scale_x_discrete(labels = c("Female, imm.",
                              "Female, mat.",
                              "Male, imm.",
                              "Male, mat."))+
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),
                                         breaks = exp(labels_logratios$y),
                                         labels = labels_logratios$label,
                                         name = "Ratio"),
                     limits = c(floor(min(df_Copeman$LR_16.0_16.1n7) / 0.5) * 0.5,
                                ceiling(max(df_Copeman$LR_16.0_16.1n7) / 0.5) * 0.5)) +
  labs(x = "Sex, Maturity",
       y = "log(16:0 / 16:1ω7)",
       color = "Sample Type",
       linetype = "Year") +
  theme_classic() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text = element_text(size = unit(8, "pt")),
        axis.title = element_text(size = unit(10, "pt")),
        axis.text.y.right = element_text(color="grey70"),
        axis.ticks.y.right = element_line(color="grey70"),
        axis.line.y.right = element_line(color="grey70"),
        axis.title.y.right = element_text(color = "grey70"),
        legend.background = element_blank(),
        legend.text = element_text(size = unit(8, "pt")),
        legend.title = element_text(size = unit(10, "pt")),
        legend.spacing = unit(0, "pt"),
        strip.text = element_text(size = unit(8, "pt")))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####
# END OF SCRIPT #####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####