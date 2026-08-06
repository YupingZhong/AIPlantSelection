
############################################################

##----3.5 plot Relationship of the richness of pollinators of subsampled networks ---------------
#         with the richness of the full network for each subsampling strategy


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
head(result_all)
result_all<-read.csv("result_all_published_phylo.csv",header=TRUE)


plot_rich_site<-result_all[,c("Study_Network_id","Abun_Top10_Rich","Abun_Top3_Rich","Abun_Top5_Rich",
                              "FlwShape_Top5_Rich","FlwShape_Top3_Rich","Total_Rich")]
head(plot_rich_site)

###################################################################

# ---------- 基础主题和颜色 ----------
my_colors <- c(
  "Top 10" = "#66C2A5", "Top 5" = "#4E79A7", "Top 3" = "#F28E2B",
  "flw 5" = "#4E79A7", "flw 3" = "#F28E2B",
  "Random 10" = "#66C2A5", "Random 5" = "#4E79A7", "Random 3" = "#F28E2B"
)

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.position = "none"
  )

# ---------- 函数：计算每个分组回归显著性 ----------
get_linetype <- function(df, x_col, y_col, group_col){
  
  df %>%
    group_by(across(all_of(group_col))) %>%
    summarise(
      lm_res = list(
        lm(
          as.formula(paste(y_col, "~", x_col)),
          data = cur_data()
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p = sapply(lm_res, \(x) summary(x)$coefficients[2,4]),
      linetype = ifelse(p < 0.05, "solid", "dashed")
    ) %>%
    dplyr::select(all_of(group_col), linetype)
}




# ---------- reshape Abundance, Flower-shape, Random ----------
rich_abun_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, Abun_Top10_Rich, Abun_Top5_Rich, Abun_Top3_Rich) %>%
  pivot_longer(cols = c(Abun_Top10_Rich, Abun_Top5_Rich, Abun_Top3_Rich),
               names_to = "Top_Plant_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(Top_Plant_Level, 
                                  levels = c("Abun_Top10_Rich","Abun_Top5_Rich","Abun_Top3_Rich"),
                                  labels = c("Top 10","Top 5","Top 3")))

rich_shape_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, FlwShape_Top5_Rich, FlwShape_Top3_Rich) %>%
  pivot_longer(cols = c(FlwShape_Top5_Rich, FlwShape_Top3_Rich),
               names_to = "flw_Plant_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(flw_Plant_Level,
                                  levels = c("FlwShape_Top5_Rich","FlwShape_Top3_Rich"),
                                  labels = c("flw 5","flw 3")))

rich_random_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, Random_10_Rich, Random_5_Rich, Random_3_Rich) %>%
  pivot_longer(cols = c(Random_10_Rich, Random_5_Rich, Random_3_Rich),
               names_to = "Random_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(Random_Level,
                                  levels = c("Random_10_Rich","Random_5_Rich","Random_3_Rich"),
                                  labels = c("Random 10","Random 5","Random 3")))

# ---------- 计算每组 linetype ----------
linetype_abun <- get_linetype(
  rich_abun_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)

linetype_shape <- get_linetype(
  rich_shape_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)

linetype_random <- get_linetype(
  rich_random_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)



# ---------- 绘图函数，增加 stat_cor ----------
plot_rich <- function(df, linetype_df, ylab,
                      cor_label_x = "left",
                      cor_label_y = NULL) {
  
  n_group <- length(unique(df$Top_Plant_Level))
  
  if (is.null(cor_label_y)) {
    cor_label_y <- seq(0.98, 0.98 - 0.02*(n_group-1), length.out = n_group)
  }
  
  ggplot(df, aes(x = Total_Rich, y = Richness_Value, color = Top_Plant_Level)) +
    geom_point(shape = 1, size = 2,
               position = position_jitter(width = 0.03, height = 0.5)) +
    geom_smooth(
      method = "lm", se = FALSE,
      aes(linetype = Top_Plant_Level),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    scale_linetype_manual(values = setNames(linetype_df$linetype, linetype_df[[1]])) +
    stat_cor(
      aes(group = Top_Plant_Level, color = Top_Plant_Level),
      method = "pearson", size = 4,
      label.x.npc = cor_label_x,
      label.y.npc = cor_label_y
    ) +
    scale_color_manual(values = my_colors) +
    labs(x = "Pollinator richness in full network", y = ylab) +
    base_theme+
    theme(
      plot.margin = margin(t = 0, r = 10, b = 0, l = 0)
    )
}


# ---------- 绘制三张图 ----------
library(cowplot)
library(ggplot2)

p_abun <- plot_rich(rich_abun_long, linetype_abun, ylab = "") +
  ggtitle("A. Abundance-based") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))


p_shape <- plot_rich(rich_shape_long, linetype_shape, ylab = "") +
  ggtitle("B. Abundance + flower shape") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))

p_random <- plot_rich(rich_random_long, linetype_random, ylab = "", cor_label_x = 0.33) +
  ggtitle("C. Random") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))

# ---------- 单独生成 legend ----------
legend_plot <- p_abun +
  labs(colour = "Number of plant species subsampled") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

legend <- cowplot::get_legend(legend_plot)


# ---------- 拼接三张图 + legend ----------
plots_row <- plot_grid(
  p_abun, p_shape, p_random,
  ncol = 3,
  align = "h",    # 水平对齐
  axis = "tb",    # 上下轴对齐
  rel_widths = c(1, 1, 1) # 三张图宽度相等
)

final_plot <- plot_grid(
  plots_row,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

# ---------- 添加大 Y 轴标题 ----------
final_plot <- ggdraw(final_plot) +
  draw_label(
    "Pollinator richness in subsampled network",
    x = 0, y = 0.5,          # 左侧居中纵向
    angle = 90,              # 竖直
    fontface = "bold",
    size = 12,
    hjust = 0.41,             # 完全水平居中
    vjust = 1.5              # 完全纵向居中
  )

final_plot


ggsave("/Chap1_TargetPlant_to_monitor/result_260526/rich_3.png", final_plot, width = 9.8, height = 5, units = "in", dpi = 300)








############################################################
## Visualization of nestedness results
############################################################

# Noted: “R” and “p” indicate the Pearson correlation coefficient and its corresponding p-value. 
#          Lines represent linear model fits for visualization purposes. 

# Define color palette for different top-N plant levels
############

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)


# read data
all_nestedness_df <- readRDS(
  "data/processed/all_nestedness_df.rds"
)

all_connectance_df<- readRDS(
  "data/processed/all_connectance_df.rds"
)

############################################################
# colors
############################################################

my_colors <- c(
  "Top 10" = "#66C2A5",
  "Top 3"  = "#E9C46A",
  "Top 5"  = "#B07AA1",
  "flw 3"  = "#E9C46A",
  "flw 5"  = "#B07AA1",
)



############################################################
# Figure 1
# Abundance + Phylogeny top N vs full network
############################################################

long_df <- all_nestedness_df %>%
  pivot_longer(
    cols = c(
      Nestedness_Abun10,
      Nestedness_Abun5,
      Nestedness_Abun3
    ),
    names_to = "Top_Plant_Level",
    values_to = "Nestedness_Value"
  ) %>%
  filter(!is.na(Nestedness_Value)) %>%
  mutate(
    Top_Plant_Level = factor(
      Top_Plant_Level,
      levels = c(
        "Nestedness_Abun10",
        "Nestedness_Abun5",
        "Nestedness_Abun3"
      ),
      labels = c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )


# regression significance

reg_significance <- long_df %>%
  group_by(Top_Plant_Level) %>%
  summarise(
    lm_res = list(
      lm(
        Nestedness_Value ~ Nestedness_orig,
        data = cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value =
      summary(lm_res)$coefficients[2,4],
    signif_line =
      ifelse(p_value < 0.05,
             "solid",
             "dashed")
  ) %>%
  select(
    Top_Plant_Level,
    signif_line
  )



set.seed(2025)

p1 <- ggplot(
  long_df,
  aes(
    x = Nestedness_orig,
    y = Nestedness_Value,
    color = Top_Plant_Level
  )
)+
  geom_point(
    shape = 1,
    size = 2,
    position = position_jitter(
      width = 0.03,
      height = 0.5
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(
      linetype = Top_Plant_Level
    )
  )+
  scale_linetype_manual(
    values =
      setNames(
        reg_significance$signif_line,
        reg_significance$Top_Plant_Level
      )
  )+
  scale_y_continuous(
    limits=c(0,60)
  )+
  stat_cor(
    aes(color=Top_Plant_Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  labs(
    y="NODF in Subnetwork",
    x="NODF of Full Network",
    color="Subsampling method",
    linetype="Regression significance"
  )+
  theme_classic(base_size=12)



############################################################
# Figure 2
# Flower shape comparison
############################################################
long_df2 <- all_nestedness_df %>%
  pivot_longer(
    cols = c(
      Nestedness_FlwShape5,
      Nestedness_FlwShape3
    ),
    names_to = "flw_Plant_Level",
    values_to = "Nestedness_Value"
  ) %>%
  filter(!is.na(Nestedness_Value)) %>%
  mutate(
    Top_Plant_Level = factor(
      flw_Plant_Level,
      levels = c(
        "Nestedness_FlwShape5",
        "Nestedness_FlwShape3"
      ),
      labels = c(
        "flw 5",
        "flw 3"
      )
    )
  )


reg_significance2 <- long_df2 %>%
  group_by(Top_Plant_Level) %>%
  summarise(
    lm_res=list(
      lm(
        Nestedness_Value~Nestedness_orig,
        data=cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise()%>%
  mutate(
    p_value=
      summary(lm_res)$coefficients[2,4],
    signif_line=
      ifelse(
        p_value<0.05,
        "solid",
        "dashed"
      )
  )%>%
  select(
    Top_Plant_Level,
    signif_line
  )



p2 <- ggplot(
  long_df2,
  aes(
    x=Nestedness_orig,
    y=Nestedness_Value,
    color=Top_Plant_Level
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.03,
      height=0.5
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(
      linetype=Top_Plant_Level
    )
  )+
  scale_linetype_manual(
    values=
      setNames(
        reg_significance2$signif_line,
        reg_significance2$Top_Plant_Level
      )
  )+
  scale_y_continuous(
    limits=c(0,60)
  )+
  stat_cor(
    aes(color=Top_Plant_Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  labs(
    y="NODF in Flower-based Subnetwork",
    x="NODF of Full Network",
    color="Subsampling method",
    linetype="Regression significance"
  )+
  theme_classic(base_size=12)

####
#Figure 3 phylo
###############3

long_df3 <- all_nestedness_df %>%
  pivot_longer(
    cols = c(
      Nestedness_Pylo10,
      Nestedness_Pylo5,
      Nestedness_Pylo3
    ),
    names_to = "Pylo_Level",
    values_to = "Nestedness_Value"
  ) %>%
  filter(!is.na(Nestedness_Value)) %>%
  mutate(
    Pylo_Level = factor(
      Pylo_Level,
      levels = c(
        "Nestedness_Pylo10",
        "Nestedness_Pylo5",
        "Nestedness_Pylo3"
      ),
      labels = c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )

reg_significance3 <- long_df3 %>%
  group_by(Pylo_Level) %>%
  summarise(
    lm_res = list(
      lm(
        Nestedness_Value ~ Nestedness_orig,
        data = cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value =
      summary(lm_res)$coefficients[2,4],
    signif_line =
      ifelse(
        p_value < 0.05,
        "solid",
        "dashed"
      )
  ) %>%
  select(
    Pylo_Level,
    signif_line
  )

p3 <- ggplot(
  long_df3,
  aes(
    x = Nestedness_orig,
    y = Nestedness_Value,
    color = Pylo_Level
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.03,
      height=0.5
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(
      linetype=Pylo_Level
    )
  )+
  scale_linetype_manual(
    values=setNames(
      reg_significance3$signif_line,
      reg_significance3$Pylo_Level
    )
  )+
  scale_y_continuous(
    limits=c(0,60)
  )+
  stat_cor(
    aes(color=Pylo_Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  labs(
    y="NODF in Phylogenetic Subnetwork",
    x="NODF of Full Network",
    color="Subsampling method",
    linetype="Regression significance"
  )+
  theme_classic(base_size=12)

############################################################
# combine
############################################################

plot_grid(
  p1,
  p2,
  p3,
  ncol=3
)



##----------- Visualization -------------

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)


############################################################
# colors (same as Nestedness)
############################################################

my_colors <- c(
  "Top 10" = "#66C2A5",
  "Top 5"  = "#B07AA1",
  "Top 3"  = "#E9C46A",
  "flw 5"  = "#B07AA1",
  "flw 3"  = "#E9C46A",
  "Phylo 10" = "#66C2A5",
  "Phylo 5"  = "#B07AA1",
  "Phylo 3"  = "#E9C46A"
)



############################################################
# Figure 4
# Abundance-based connectance
############################################################

long_connect_abun <- all_connectance_df %>%
  pivot_longer(
    cols = c(
      Connectance_Abun10,
      Connectance_Abun5,
      Connectance_Abun3
    ),
    names_to="Level",
    values_to="Connectance_Value"
  ) %>%
  filter(!is.na(Connectance_Value)) %>%
  mutate(
    Level=factor(
      Level,
      levels=c(
        "Connectance_Abun10",
        "Connectance_Abun5",
        "Connectance_Abun3"
      ),
      labels=c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )


reg_abun <- long_connect_abun %>%
  group_by(Level) %>%
  summarise(
    lm_res=list(
      lm(
        Connectance_Value ~ Connectance_orig,
        data=cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise()%>%
  mutate(
    p_value=summary(lm_res)$coefficients[2,4],
    signif_line=ifelse(
      p_value<0.05,
      "solid",
      "dashed"
    )
  )


p_connect_abun <- ggplot(
  long_connect_abun,
  aes(
    x=Connectance_orig,
    y=Connectance_Value,
    color=Level
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.005,
      height=0.005
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(linetype=Level)
  )+
  scale_linetype_manual(
    values=setNames(
      reg_abun$signif_line,
      reg_abun$Level
    )
  )+
  stat_cor(
    aes(color=Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  scale_y_continuous(
    limits=c(0,1)
  )+
  labs(
    x="Connectance of Full Network",
    y="Connectance in Abundance-based Subnetwork",
    color="Subsampling method",
    linetype="Significance"
  )+
  theme_classic(base_size=12)



############################################################
# Figure 5
# Flower shape connectance
############################################################

long_connect_flw <- all_connectance_df %>%
  pivot_longer(
    cols=c(
      Connectance_FlwShape5,
      Connectance_FlwShape3
    ),
    names_to="Level",
    values_to="Connectance_Value"
  ) %>%
  filter(!is.na(Connectance_Value)) %>%
  mutate(
    Level=factor(
      Level,
      levels=c(
        "Connectance_FlwShape5",
        "Connectance_FlwShape3"
      ),
      labels=c(
        "flw 5",
        "flw 3"
      )
    )
  )


reg_flw <- long_connect_flw %>%
  group_by(Level)%>%
  summarise(
    lm_res=list(
      lm(
        Connectance_Value~Connectance_orig,
        data=cur_data()
      )
    ),
    .groups="drop"
  )%>%
  rowwise()%>%
  mutate(
    p_value=summary(lm_res)$coefficients[2,4],
    signif_line=ifelse(
      p_value<0.05,
      "solid",
      "dashed"
    )
  )


p_connect_flw <- ggplot(
  long_connect_flw,
  aes(
    x=Connectance_orig,
    y=Connectance_Value,
    color=Level
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.005,
      height=0.005
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(linetype=Level)
  )+
  scale_linetype_manual(
    values=setNames(
      reg_flw$signif_line,
      reg_flw$Level
    )
  )+
  stat_cor(
    aes(color=Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  scale_y_continuous(
    limits=c(0,1)
  )+
  labs(
    x="Connectance of Full Network",
    y="Connectance in Flower-based Subnetwork",
    color="Subsampling method",
    linetype="Significance"
  )+
  theme_classic(base_size=12)



############################################################
# Figure 6
# Phylogenetic connectance
############################################################

long_connect_phylo <- all_connectance_df %>%
  pivot_longer(
    cols=c(
      Connectance_Pylo10,
      Connectance_Pylo5,
      Connectance_Pylo3
    ),
    names_to="Level",
    values_to="Connectance_Value"
  ) %>%
  filter(!is.na(Connectance_Value)) %>%
  mutate(
    Level=factor(
      Level,
      levels=c(
        "Connectance_Pylo10",
        "Connectance_Pylo5",
        "Connectance_Pylo3"
      ),
      labels=c(
        "Phylo 10",
        "Phylo 5",
        "Phylo 3"
      )
    )
  )


p_connect_phylo <- ggplot(
  long_connect_phylo,
  aes(
    x=Connectance_orig,
    y=Connectance_Value,
    color=Level
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.005,
      height=0.005
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(linetype=Level)
  )+
  stat_cor(
    aes(color=Level),
    method="pearson",
    size=4
  )+
  scale_color_manual(
    values=my_colors
  )+
  scale_y_continuous(
    limits=c(0,1)
  )+
  labs(
    x="Connectance of Full Network",
    y="Connectance in Phylogenetic Subnetwork",
    color="Subsampling method",
    linetype="Significance"
  )+
  theme_classic(base_size=12)



############################################################
# check plots
############################################################

p_connect_abun
p_connect_flw
p_connect_phylo

############################################################
# Random-N plant subsampling visualization
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


############################################################
# Read random results
############################################################

random_all <- readRDS(
  "data/processed/random_all_metrics.rds"
) %>%
  mutate(
    Method = trimws(Method),
    Method = case_when(
      Method %in% c("random10", "Random10", "random 10") ~ "Top 10",
      Method %in% c("random5", "Random5", "random 5") ~ "Top 5",
      Method %in% c("random3", "Random3", "random 3") ~ "Top 3",
      TRUE ~ Method
    )
  )


random_all$Method <- factor(
  random_all$Method,
  levels = c(
    "Top 10",
    "Top 5",
    "Top 3"
  )
)



############################################################
# Prepare data
############################################################


random_nodf_df <- random_all %>%
  select(
    Study_Network_id,
    Method,
    Random_NODF
  )


random_conn_df <- random_all %>%
  select(
    Study_Network_id,
    Method,
    Random_Connectance
  )



# NODF

nodf_random_df <- all_nestedness_df %>%
  select(
    Study_Network_id,
    Nestedness_orig
  ) %>%
  left_join(
    pivot_wider(
      random_nodf_df,
      names_from = Method,
      values_from = Random_NODF,
      values_fn = mean
    ),
    by = "Study_Network_id"
  )



# Connectance

conn_random_df <- all_connectance_df %>%
  select(
    Study_Network_id,
    Connectance_orig
  ) %>%
  left_join(
    pivot_wider(
      random_conn_df,
      names_from = Method,
      values_from = Random_Connectance,
      values_fn = mean
    ),
    by = "Study_Network_id"
  )




############################################################
# Long format
############################################################


random_long_nodf <- nodf_random_df %>%
  pivot_longer(
    cols = c(
      `Top 10`,
      `Top 5`,
      `Top 3`
    ),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )



random_long_conn <- conn_random_df %>%
  pivot_longer(
    cols = c(
      `Top 10`,
      `Top 5`,
      `Top 3`
    ),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )



############################################################
# regression significance
############################################################


random_nodf_sig <- random_long_nodf %>%
  group_by(Method) %>%
  summarise(
    lm_res = list(
      lm(
        Value ~ Nestedness_orig,
        data = cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value =
      summary(lm_res)$coefficients[2,4],
    signif_line =
      ifelse(
        p_value < 0.05,
        "solid",
        "dashed"
      )
  ) %>%
  select(
    Method,
    signif_line
  )



random_conn_sig <- random_long_conn %>%
  group_by(Method) %>%
  summarise(
    lm_res = list(
      lm(
        Value ~ Connectance_orig,
        data = cur_data()
      )
    ),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value =
      summary(lm_res)$coefficients[2,4],
    signif_line =
      ifelse(
        p_value < 0.05,
        "solid",
        "dashed"
      )
  ) %>%
  select(
    Method,
    signif_line
  )



############################################################
# Colors (same as all figures)
############################################################

random_colors <- c(
  "Top 10" = "#66C2A5",
  "Top 5"  = "#B07AA1",
  "Top 3"  = "#E9C46A"
)



############################################################
# Random NODF plot
############################################################


p_random_nodf <- ggplot(
  random_long_nodf,
  aes(
    x = Nestedness_orig,
    y = Value,
    color = Method
  )
)+
  geom_point(
    shape = 1,
    size = 2,
    position = position_jitter(
      width = 0.02,
      height = 0.2
    )
  )+
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    aes(
      linetype = Method
    )
  )+
  scale_linetype_manual(
    values=setNames(
      random_nodf_sig$signif_line,
      random_nodf_sig$Method
    )
  )+
  stat_cor(
    aes(color=Method),
    method="pearson",
    size=4
  )+
  scale_y_continuous(
    limits=c(0,60)
  )+
  scale_color_manual(
    values=random_colors
  )+
  labs(
    x="NODF of Full Network",
    y="NODF in Random-N Subnetwork",
    color="Method",
    linetype="Significance"
  )+
  theme_classic(base_size=12)



############################################################
# Random Connectance plot
############################################################


p_random_conn <- ggplot(
  random_long_conn,
  aes(
    x = Connectance_orig,
    y = Value,
    color = Method
  )
)+
  geom_point(
    shape=1,
    size=2,
    position=position_jitter(
      width=0.005,
      height=0.005
    )
  )+
  geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8,
    aes(
      linetype=Method
    )
  )+
  scale_linetype_manual(
    values=setNames(
      random_conn_sig$signif_line,
      random_conn_sig$Method
    )
  )+
  stat_cor(
    aes(color=Method),
    method="pearson",
    size=4
  )+
  scale_y_continuous(
    limits=c(0,1)
  )+
  scale_color_manual(
    values=random_colors
  )+
  labs(
    x="Connectance of Full Network",
    y="Connectance in Random-N Subnetwork",
    color="Method",
    linetype="Significance"
  )+
  theme_classic(base_size=12)



# display

p_random_nodf
p_random_conn


############################################################
# Combine all network metric plots
############################################################

library(cowplot)


############################################################
# Remove legends
############################################################

# Connectance
p_connect_abun.1 <- p_connect_abun +
  theme(legend.position = "none")

p_connect_flw.1 <- p_connect_flw +
  theme(legend.position = "none")

p_connect_phylo.1 <- p_connect_phylo +
  theme(legend.position = "none")

p_random_conn.1 <- p_random_conn +
  theme(legend.position = "none")


# Nestedness

p_abun_nodf.1 <- p1 +
  theme(legend.position = "none")

p_flw_nodf.1 <- p2 +
  theme(legend.position = "none")

p_phylo_nodf.1 <- p3 +
  theme(legend.position = "none")

p_random_nodf.1 <- p_random_nodf +
  theme(legend.position = "none")



############################################################
# Arrange Connectance row
############################################################


connectance_row <- plot_grid(
  
  p_connect_abun.1,
  p_connect_flw.1,
  p_connect_phylo.1,
  p_random_conn.1,
  
  ncol = 4,
  align = "hv",
  axis = "tblr"
)



############################################################
# Arrange Nestedness row
############################################################


nestedness_row <- plot_grid(
  
  p_abun_nodf.1,
  p_flw_nodf.1,
  p_phylo_nodf.1,
  p_random_nodf.1,
  
  ncol = 4,
  align = "hv",
  axis = "tblr"
)



############################################################
# Unified legend
############################################################


legend_plot <- p1 +
  theme(
    legend.position="bottom"
  ) +
  guides(
    color = guide_legend(
      nrow = 1
    ),
    linetype = guide_legend(
      nrow = 1
    )
  )


legend <- cowplot::get_legend(
  legend_plot
)



############################################################
# Final figure
############################################################


combined_plot <- plot_grid(
  
  connectance_row,
  nestedness_row,
  legend,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1,
    0.12
  )
)



combined_plot


ggsave(
  "./result_260526/network.matrix.png",
  combined_plot,
  width = 18,
  height = 9.5,
  units = "in",
  dpi = 300
)