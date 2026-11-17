# Replacement for 0.2.plotting.network.locations.R
# Source this file in RStudio. Project inputs are never changed.
# Default: the SAME processed dataset as the original map (267 annual networks).
# To map the full database instead, set dataset_scope <- 'raw' below.
project_dir <- 'D:/Chap1_TargetPlant_to_monitor/Project_PlantSelection'
output_dir <- 'D:/Chap1_TargetPlant_to_monitor/result_revision'
dataset_scope <- 'processed'  # 'processed' or 'raw'
required <- c('dplyr','ggplot2','sf','giscoR','patchwork')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf); library(patchwork)})
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
stopifnot(dataset_scope %in% c('processed','raw'))
input <- if(dataset_scope=='processed') 'data/processed/data_interact_published.rds' else 'data/raw/Interaction_data_published.rds'
d <- readRDS(file.path(project_dir,input))
if(dataset_scope=='processed') {
  d$Network_key <- as.character(d$Study_Network_id)
} else {
  # Same study/network unit as the continental classifier analysis; years pooled.
  stopifnot(!anyNA(d$Study_id),!anyNA(d$Network_id))
  d$Network_key <- paste(d$Study_id,d$Network_id,sep='_')
}
stopifnot(!anyNA(d$Network_key),all(nzchar(d$Network_key)))
d <- d %>% mutate(Bioregion=trimws(as.character(Bioregion)),
  Bioregion=if_else(is.na(Bioregion)|Bioregion=='','Unknown',Bioregion),
  Latitude=suppressWarnings(as.numeric(as.character(Latitude))),
  Longitude=suppressWarnings(as.numeric(as.character(Longitude))))
# Do not silently count a network in multiple bioregions.
regions <- d %>% distinct(Network_key,Bioregion)
conflicts <- regions %>% count(Network_key) %>% filter(n>1)
if(nrow(conflicts)) {
 write.csv(semi_join(regions,conflicts,by='Network_key'),file.path(output_dir,'CHECK_bioregion_conflicts.csv'),row.names=FALSE)
 stop('Networks have conflicting bioregions; inspect CHECK_bioregion_conflicts.csv.')
}
counts <- regions %>% count(Bioregion,name='Networks') %>% arrange(desc(Networks)) %>%
 mutate(Percent=100*Networks/sum(Networks))
stopifnot(sum(counts$Networks)==n_distinct(d$Network_key))
# Retain all distinct reported locations, rather than arbitrarily choosing one per study.
locations <- d %>% distinct(Network_key,Bioregion,Latitude,Longitude) %>%
 filter(is.finite(Latitude),is.finite(Longitude),abs(Latitude)<=90,abs(Longitude)<=180)
missing_coords <- anti_join(regions,distinct(locations,Network_key),by='Network_key')
points <- locations %>% distinct(Bioregion,Latitude,Longitude)
palette <- c(Alpine='#9864B5',Atlantic='#468AC5',Boreal='#309C9C',Continental='#50934B',
 Mediterranean='#E5A522',Pannonian='#CB7968',Steppic='#BCB363',Unknown='#777777',
 Arctic='#8DAAC5',Anatolian='#9B6945','Black Sea'='#D57DA4',Macaronesian='#7C727E')
unmapped <- setdiff(counts$Bioregion,names(palette))
if(length(unmapped)) palette <- c(palette,setNames(grDevices::hcl.colors(length(unmapped),'Dark 3'),unmapped))
# Cache the basemap: first run may require internet; subsequent runs are offline.
cache <- file.path(output_dir,'europe_countries_basemap.rds')
countries <- if(file.exists(cache)) readRDS(cache) else {
 x <- tryCatch(giscoR::gisco_get_countries(resolution='20',year='2020',epsg='4326'),error=function(e) NULL)
 if(is.null(x)) {
   if(!requireNamespace('maps',quietly=TRUE)) stop('GISCO unavailable; install maps for an offline basemap.')
   message('GISCO unavailable: using offline maps::world country boundaries.')
   x <- sf::st_as_sf(maps::map('world',plot=FALSE,fill=TRUE))
   x <- sf::st_transform(x,4326)
 }
 saveRDS(x,cache); x
}
map <- ggplot(countries)+
 geom_sf(fill='floralwhite',colour='grey45',linewidth=0.18)+
 geom_point(data=points,aes(Longitude,Latitude,fill=Bioregion),shape=21,size=2,stroke=0.25,colour='grey20',alpha=0.9)+
 scale_fill_manual(values=palette,breaks=counts$Bioregion,drop=FALSE)+
 coord_sf(crs=st_crs(4326),xlim=c(-13,39),ylim=c(34,71),expand=FALSE)+
 scale_x_continuous(breaks=c(-10,0,10,20,30))+scale_y_continuous(breaks=c(40,50,60,70))+
 labs(x='Longitude',y='Latitude')+theme_bw(base_size=15)+
 theme(panel.background=element_rect(fill='aliceblue'),
       panel.grid=element_blank(),legend.position='none',
       plot.tag = element_text(face = "bold", size = 22),
       axis.text = element_text(size = 13, colour = "black"),
       axis.title = element_text(size = 16),
       plot.margin = margin(3, 0, 2, 0))

bar <- ggplot(counts,aes(Networks,reorder(Bioregion,Networks),fill=Bioregion))+
 geom_col(width=0.67)+geom_text(aes(label=Networks),hjust=-0.2,size=5)+
 scale_fill_manual(values=palette)+scale_x_continuous(
   limits = c(0, max(counts$Networks) * 1.35),
   breaks = c(0, 40, 80, 120),
   expand = expansion(mult = c(0, 0))
 ) +
 labs(x='Number of networks',y=NULL)+theme_bw(base_size=15)+
 theme(legend.position='none',
       panel.grid.major.y=element_blank(),
       panel.grid.minor=element_blank(),
       plot.tag = element_text(face = "bold", size = 22),
       axis.text = element_text(size = 13, colour = "black"),
       axis.title = element_text(size = 15),
       plot.margin = margin(3, 16, 2, 0))

figure <- cowplot::ggdraw() +
  cowplot::draw_plot(
    map, x = 0, y = 0, width = 0.64, height = 0.96
  ) +
  cowplot::draw_plot(
    bar, x = 0.64, y = 0, width = 0.36, height = 0.96
  ) +
  cowplot::draw_plot_label(
    c("a", "b"),
    x = c(0.006, 0.645),
    y = c(0.995, 0.995),
    size = 24, fontface = "bold",
    hjust = 0, vjust = 1
  )

prefix <- paste0('Figure_S11_bioregions_',dataset_scope)
ggsave(file.path(output_dir,paste0(prefix,'.png')),figure,width=10,height=6.3,dpi=400,bg='white')
ggsave(file.path(output_dir,paste0(prefix,'.pdf')),figure,width=10,height=6.3,device=grDevices::cairo_pdf,bg='white')
