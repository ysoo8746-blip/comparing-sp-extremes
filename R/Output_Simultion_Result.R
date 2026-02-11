######################################################
## Plot Multiple test results
######################################################
source('/Functions/Draw_ggplot.R')
library(ggpubr)
library(ggplot2)

m1= c('BH', 'JS', 'Mirror', 'BiCLfdr-1');m2= c('BH', 'JS', 'Mirror', 'BiCLfdr-2')
m4= c('BH', 'JS', 'Mirror', 'BiCLfdr-4');m12= c('BH', 'JS', 'Mirror', 'BiCLfdr-12')
m= c('BH', 'JS', 'Mirror', 'BiCLfdr')
m_n0 = c("BH", "JS", "Mirror","BiCLfdr-4", "BiCLfdr-12")
m_n = c(m1,  m4[4], m12[4])
s = c("Weak", "Moderate", "Strong")

###### Simulation Boxplot results for Location parameter ########
Loc_Multi_Strong2 = MyMulti_outputNeiLoc("Loc","Strong",Method=m2, Nei=2)
Loc_Multi_Mod2 = MyMulti_outputNeiLoc("Loc","Moderate",Method=m2, Nei=2)
Loc_Multi_Weak2 = MyMulti_outputNeiLoc("Loc","Weak",Method=m2, Nei=2)

Loc_Multi_Strong4 = MyMulti_outputNeiLoc("Loc","Strong",Method=m4, Nei=4)
Loc_Multi_Mod4 = MyMulti_outputNeiLoc("Loc","Moderate",Method=m4, Nei=4)
Loc_Multi_Weak4 = MyMulti_outputNeiLoc("Loc","Weak",Method=m4, Nei=4)

Loc_Multi_Strong12 = MyMulti_outputNeiLoc("Loc","Strong",Method=m12, Nei=12)
Loc_Multi_Mod12 = MyMulti_outputNeiLoc("Loc","Moderate",Method=m12, Nei=12)
Loc_Multi_Weak12 = MyMulti_outputNeiLoc("Loc","Weak",Method=m12, Nei=12)

Loc_ResOpt_2 = list(Weak = Loc_Multi_Weak2, Moderate= Loc_Multi_Mod2, Strong=Loc_Multi_Strong2)
Loc_ResOpt_4 = list(Weak = Loc_Multi_Weak4, Moderate= Loc_Multi_Mod4, Strong=Loc_Multi_Strong4)
Loc_ResOpt_12 = list(Weak = Loc_Multi_Weak12, Moderate= Loc_Multi_Mod12, Strong=Loc_Multi_Strong12)

### ggplot output ### 
G.Loc_N2_4 = My_GG_MultiTestLoc(Loc_ResOpt_2, Loc_ResOpt_4, "Loc", Method=m_n0, Sig=s)
G.Loc_N4_12 = My_GG_MultiTestLoc(Loc_ResOpt_4, Loc_ResOpt_12, "Loc", Method=m_n0, Sig=s)

###### Simulation Boxplot results for Scale parameter ########
Scale_Multi_Strong2 = MyMulti_outputNeiN(Param="Scale",signal="Strong",Method=m2, Nei=2)
Scale_Multi_Mod2 = MyMulti_outputNeiN(Param="Scale",signal="Moderate",Method=m2, Nei=2)
Scale_Multi_Weak2 = MyMulti_outputNeiN(Param="Scale",signal="Weak",Method=m2, Nei=2)

Scale_Multi_Strong4 = MyMulti_outputNeiN(Param="Scale",signal="Strong",Method=m4, Nei=4)
Scale_Multi_Mod4 = MyMulti_outputNeiN(Param="Scale",signal="Moderate",Method=m4, Nei=4)
Scale_Multi_Weak4 = MyMulti_outputNeiN(Param="Scale",signal="Weak",Method=m4, Nei=4)

Scale_Multi_Strong12 = MyMulti_outputNeiN(Param="Scale",signal="Strong",Method=m12, Nei=12)
Scale_Multi_Mod12 = MyMulti_outputNeiN(Param="Scale",signal="Moderate",Method=m12, Nei=12)
Scale_Multi_Weak12 = MyMulti_outputNeiN(Param="Scale",signal="Weak",Method=m12, Nei=12)

Scale_ResOpt_2 = list(Weak = Scale_Multi_Weak2, Moderate= Scale_Multi_Mod2, Strong=Scale_Multi_Strong2)
Scale_ResOpt_4 = list(Weak = Scale_Multi_Weak4, Moderate= Scale_Multi_Mod4, Strong=Scale_Multi_Strong4)
Scale_ResOpt_12 = list(Weak = Scale_Multi_Weak12, Moderate= Scale_Multi_Mod12, Strong=Scale_Multi_Strong12)

### ggplot output ### 
G.Scale_N2_4 = My_GG_MultiTest(Scale_ResOpt_2,  Scale_ResOpt_4, "Scale", Method=m_n0, Sig=s)
G.Scale_N4_12 = My_GG_MultiTest(Scale_ResOpt_4,  Scale_ResOpt_12, "Scale", Method=m_n0, Sig=s)

###### Simulation Boxplot results for Shape parameter ########
Shape_Multi_Strong2 = MyMulti_outputNeiN(Param="Shape",signal="Strong",Method=m2, Nei=2)
Shape_Multi_Mod2 = MyMulti_outputNeiN(Param="Shape",signal="Moderate",Method=m2, Nei=2)
Shape_Multi_Weak2 = MyMulti_outputNeiN(Param="Shape",signal="Weak",Method=m2, Nei=2)

Shape_Multi_Strong4 = MyMulti_outputNeiN(Param="Shape",signal="Strong",Method=m4, Nei=4)
Shape_Multi_Mod4 = MyMulti_outputNeiN(Param="Shape",signal="Moderate",Method=m4, Nei=4)
Shape_Multi_Weak4 = MyMulti_outputNeiN(Param="Shape",signal="Weak",Method=m4, Nei=4)

Shape_Multi_Strong12 = MyMulti_outputNeiN(Param="Shape",signal="Strong",Method=m12, Nei=12)
Shape_Multi_Mod12 = MyMulti_outputNeiN(Param="Shape",signal="Moderate",Method=m12, Nei=12)
Shape_Multi_Weak12 = MyMulti_outputNeiN(Param="Shape",signal="Weak",Method=m12, Nei=12)

Shape_ResOpt_2 = list(Weak = Shape_Multi_Weak2, Moderate= Shape_Multi_Mod2, Strong=Shape_Multi_Strong2)
Shape_ResOpt_4 = list(Weak = Shape_Multi_Weak4, Moderate= Shape_Multi_Mod4, Strong=Shape_Multi_Strong4)
Shape_ResOpt_12 = list(Weak = Shape_Multi_Weak12, Moderate= Shape_Multi_Mod12, Strong=Shape_Multi_Strong12)

### ggplot output ### 
G.Shape_N2_4 = My_GG_MultiTest(Shape_ResOpt_2,  Shape_ResOpt_4, "Shape", Method=m_n0, Sig=s)
G.Shape_N4_12 = My_GG_MultiTest(Shape_ResOpt_4,  Shape_ResOpt_12, "Shape", Method=m_n0, Sig=s)

############ Parameter and return level estimation ############ 
##### ---- Figure 3 : Scale Difference ----- ####
# Weak Signal
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig2')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Weak = rt.diff.GEVP[,2];R50_Scale_diff_Weak = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Weak = param.diff.GEVP[,2] # j = 8

# Moderate Signal 
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig2.5')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Moderate = rt.diff.GEVP[,2];R50_Scale_diff_Moderate = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Moderate = param.diff.GEVP[,2] # j = 8

# Strong Signal 
setwd('Result/Scale_Rt_Size2Diff_P33_R3S0.2B9_Sig3')
datafile_names0 = list.files();load(datafile_names0[7])
R10_Scale_diff_Strong = rt.diff.GEVP[,2];R50_Scale_diff_Strong = rt.diff.GEVP[,4] ## j = 1/ j = 7
Scale_diff_Strong = param.diff.GEVP[,2] # j = 8
# Apply the function to each case
R10_df_weak     <- reshape_loc_diff(R10_Scale_diff_Weak, "Weak", delta.tr)
R10_df_moderate <- reshape_loc_diff(R10_Scale_diff_Moderate, "Moderate", delta.tr)
R10_df_strong   <- reshape_loc_diff(R10_Scale_diff_Strong, "Strong", delta.tr)

R50_df_weak     <- reshape_loc_diff(R50_Scale_diff_Weak, "Weak", delta.tr)
R50_df_moderate <- reshape_loc_diff(R50_Scale_diff_Moderate, "Moderate", delta.tr)
R50_df_strong   <- reshape_loc_diff(R50_Scale_diff_Strong, "Strong", delta.tr)

# Build one combined data frame
plot_df <- bind_rows(
  make_df(R10_Scale_diff_Weak,     R50_Scale_diff_Weak,     "Weak",     delta.tr),
  make_df(R10_Scale_diff_Moderate, R50_Scale_diff_Moderate, "Moderate", delta.tr),
  make_df(R10_Scale_diff_Strong,   R50_Scale_diff_Strong,   "Strong",   delta.tr)
)

# Plot: 2 (rows: 10/50-year) × 3 (cols: Weak/Moderate/Strong)
# Note: `scales = "free_y"` lets the 10-year and 50-year rows use different y-axis ranges.
plot_df$delta <- factor(plot_df$delta, levels = c("Null", "Alternative"))
RT_scalediff = ggplot(plot_df, aes(x = loc_id, y = diff)) +
  
  # Null locations (background)
  geom_point(
    data = subset(plot_df, delta == "Null"),
    aes(shape = delta),
    color = "grey60",
    size  = 1.0,
    alpha = 0.7
  ) +
  
  # Alternative locations (foreground)
  geom_point(
    data = subset(plot_df, delta == "Alternative"),
    aes(shape = delta),
    color = "black",
    size  = 1.6,
    alpha = 0.9
  ) +
  
  scale_shape_manual(
    values = c("Null" = 1,        # open circle
               "Alternative" = 16 # filled circle
    ),
    name = "Location Type"
  ) +
  
  facet_grid(horizon ~ strength, scales = "free_y") +
  
  labs(
    title = NULL,
    x = "Location Index",
    y = "Return-Level Difference"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))

##### ---- Figure 5 : Shape Difference ----- ####
# Weak signal
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig0.8')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Weak = rt.diff.GEVP[,2];R50_Shape_diff_Weak = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Weak = param.diff.GEVP[,3] # j = 3

# Moderate Signal
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig1')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Moderate = rt.diff.GEVP[,2];R50_Shape_diff_Moderate = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Moderate = param.diff.GEVP[,3] # j = 3

# Strong Signal 
setwd('Result/Shape_Rt_Size2Diff_P33_R3S0.2B9_Sig1.2')
datafile_names0 = list.files();load(datafile_names0[4])
R10_Shape_diff_Strong = rt.diff.GEVP[,2];R50_Shape_diff_Strong = rt.diff.GEVP[,4] # j =3/j= 4
Shape_diff_Strong = param.diff.GEVP[,3] # j = 3

# Apply the function to each case
R10_df_weak     <- reshape_loc_diff(R10_Shape_diff_Weak, "Weak", delta.tr)
R10_df_moderate <- reshape_loc_diff(R10_Shape_diff_Moderate, "Moderate", delta.tr)
R10_df_strong   <- reshape_loc_diff(R10_Shape_diff_Strong, "Strong", delta.tr)

R50_df_weak     <- reshape_loc_diff(R50_Shape_diff_Weak, "Weak", delta.tr)
R50_df_moderate <- reshape_loc_diff(R50_Shape_diff_Moderate, "Moderate", delta.tr)
R50_df_strong   <- reshape_loc_diff(R50_Shape_diff_Strong, "Strong", delta.tr)

# Build one combined data frame
plot_df <- bind_rows(
  make_df(R10_Shape_diff_Weak,     R50_Shape_diff_Weak,     "Weak",     delta.tr),
  make_df(R10_Shape_diff_Moderate, R50_Shape_diff_Moderate, "Moderate", delta.tr),
  make_df(R10_Shape_diff_Strong,   R50_Shape_diff_Strong,   "Strong",   delta.tr)
)

# Plot: 2 (rows: 10/50-year) × 3 (cols: Weak/Moderate/Strong)
# Note: `scales = "free_y"` lets the 10-year and 50-year rows use different y-axis ranges.
plot_df$delta <- factor(plot_df$delta, levels = c("Null", "Alternative"))
RT_shapediff = ggplot(plot_df, aes(x = loc_id, y = diff)) +
  
  # Null locations (background)
  geom_point(
    data = subset(plot_df, delta == "Null"),
    aes(shape = delta),
    color = "grey60",
    size  = 1.0,
    alpha = 0.7
  ) +
  
  # Alternative locations (foreground)
  geom_point(
    data = subset(plot_df, delta == "Alternative"),
    aes(shape = delta),
    color = "black",
    size  = 1.6,
    alpha = 0.9
  ) +
  
  scale_shape_manual(
    values = c("Null" = 1,        # open circle
               "Alternative" = 16 # filled circle
    ),
    name = "Location Type"
  ) +
  
  facet_grid(horizon ~ strength, scales = "free_y") +
  
  labs(
    title = NULL,
    x = "Location Index",
    y = "Return-Level Difference"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))



