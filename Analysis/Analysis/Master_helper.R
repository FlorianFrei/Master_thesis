Align_behv <- function(Data,State){
  # a function to center Data on specific Bhevaioural components
  # finds the min time_bin in the selected State per Trial and calculates a rel_time to this min time
  filter2 <- Data %>% group_by(Trial2) %>%  filter(Trialtype==State) %>% slice_min(time_bin) %>% mutate(minT = time_bin) %>% select(minT,id,Trial) %>% distinct()
  Data <- Data  %>% left_join(filter2) %>% mutate(rel_time = round(time_bin - minT,2) )
  Data
}


load_data <- function(datapath,metapath){
  # loads individual data (from raw_ks_to_df) and meta data (from phy output) and combines them into a df
  # takes names from filename. need to be M[Mouse_id]_[recording session] e.g. M1_2
  list_of_data <- list.files(path = datapath,
                             pattern = "\\.csv$",
                             full.names = TRUE)
  Data <- read_csv(list_of_data, id = 'group')
  Data <- Data %>% mutate(id = str_extract(Data$group, "M\\d+_\\d+"))
  
  list_of_meta <- list.files(path = metapath,
                             pattern = "\\.tsv$",
                             full.names = TRUE)
  meta <- read_tsv(list_of_meta, id = 'id2') %>%  as_tibble()%>% select(id2,cluster_id,n_spikes,depth)
  meta <- meta %>% mutate(id = str_extract(meta$id2, "M\\d+_\\d+"))
  
  Data <- Data %>% select(cluster_id:trial,id) %>% left_join(meta %>% select(!id2),by = join_by(id,cluster_id))
  # adds more data such as the behaviour,depth or Hit/miss based on which Trialtypes where present 
  Data <- Data %>% mutate(Trial = sprintf('%03d',trial))
  Data <- Data %>% unite('Trial2',c(id,Trial),remove=F) %>% group_by(Trial2) %>%  mutate(Behv = case_when(
    c('Audio','W2T_Audio') %in% Trialtype ~ "W2T",
    c('Airpuff','A2L_Audio') %in% Trialtype ~ "A2L",
    c('PC_Audio','Airpuff2') %in% Trialtype ~"PC", 
    .default = "WTF")) %>% mutate(succ = ifelse('HIT' %in% Trialtype,'Hit','Miss')) %>% ungroup() %>%  unite('cluster_id2',c(cluster_id,id),remove=F) %>% mutate(depthcat = ifelse(depth < ((1000)*(2/3)),"L5",'L3'))
  Data
}

Align_all <- function(Data){

  contAir <-Align_behv(Data = Data %>% filter(Behv =='A2L'),State = 'Airpuff') %>% mutate(State = 'Airpuff')
  contWT_Hits <- Align_behv(Data = Data %>% filter(Behv =='W2T' & succ =='Hit'),State = 'HIT') %>% mutate(State = 'Hit')
  contPC <- Align_behv(Data = Data %>% filter(Behv =='PC'),State = 'Airpuff2') %>% mutate(State = 'PC')
 
  combcont <- rbind(contAir,contWT_Hits,contPC)
  combcont <- combcont %>% left_join(combcont %>% group_by(id,Behv) %>% distinct(Trial) %>% summarise(n_trials = n()))
  combcont
}


whisk_vs_hit_df <- function(Data, whiskpath){
  whiskstart <- read_csv(whiskpath)
  contWT_Hits <- combcont %>% filter(Behv='W2T')
  contWhisk <- contWT_Hits %>% filter(Trial2 %in% whiskstart$filt) %>% left_join(whiskstart, by=join_by(Trial2 == filt)) %>% mutate(TTT = (500-mframe)/200, rel_time = rel_time - TTT) %>% mutate(State = 'whisk', Behv = 'W2T_W')
  wvh <- rbind(contWT_Hits,contWhisk)
  wvh <- wvh %>% left_join(wvh %>% group_by(id,Behv) %>% distinct(Trial) %>% summarise(n_trials = n()))
  wvh
}

avg_trials <- function(Data,bin_size=0.01){
  Trialless <- combcont %>% group_by(Behv,rel_time, cluster_id2)  %>% 
  mutate(rate = (sum(event_count)/n_trials), rate = rate/bin_size) %>% ungroup() %>% 
  distinct(depthcat,Behv,rel_time,cluster_id2,rate,depth,n_spikes,id) %>% group_by(Behv,cluster_id2) %>% 
  mutate(Zscore = (rate - mean(rate))/sd(rate))
}


Plot_overview <- function(Trialless,orderby_Behv = 'A2L',by_time = F,theme_smo = theme_minimal()) {
  if(by_time ==F) {
    ord <- Trialless  %>% filter(between(rel_time,0,0.1)) %>% group_by(cluster_id2)  %>%
      filter(Behv == orderby_Behv) %>% group_by(cluster_id2) %>%
      summarise(ord = mean(Zscore)) %>%   arrange( ord)  %>% mutate(cluster_id2 = as.factor(cluster_id2)) %>% pull(cluster_id2)
  }
  if(by_time ==T){
    ord <-  Trialless %>% group_by(cluster_id2)  %>%
      filter(between(rel_time,-1,1)) %>%  filter(Behv == orderby_Behv ) %>% group_by(cluster_id2) %>% 
      slice_max(Zscore,n=1,with_ties=F) %>% summarise(ord = rel_time) %>%   arrange( ord)  %>% distinct() %>% 
      mutate(cluster_id2 = as.factor(cluster_id2)) %>% pull(cluster_id2)
  }
  
  
  
  
  
  # z score and filter
  g1 <- Trialless  %>% filter(between(rel_time,-1,1)) %>% select(Behv,cluster_id2,Zscore,rel_time)
  

  
  #apply order 
  g1$cluster_id2 <- factor(g1$cluster_id2, levels = ord)
  
  
  #heatmap
  map <-g1 %>% ggplot(aes(rel_time,cluster_id2,fill=Zscore))+
    geom_tile()+
    geom_vline(xintercept = 0, col ='grey',linetype = 'dashed')+
    facet_wrap(~Behv,scales = 'free')+
    scale_fill_gradient2(low = 'darkblue', mid = 'beige', high = 'red', midpoint = 0, limits = c(round(min(g1$Zscore)), round(max(g1$Zscore)*0.5)), breaks = c(round(min(g1$Zscore)),  round(max(g1$Zscore)*0.5)) ,na.value = "darkred")+
    theme_void()+
    theme(plot.margin = unit(c(0,0.2,0,1), 'lines'),
          strip.text.x = element_text(size = 30))+
    labs(fill='Z- scored FR') 
  
  
  #activity across all units
  smo <- g1 %>% group_by(Behv,rel_time) %>% summarise( rate3 = mean(Zscore), sd = sd(Zscore), se = sd/sqrt(n())) %>% ggplot(aes(rel_time,rate3))+
    geom_line()+
    geom_ribbon(aes(ymin = rate3-se, ymax = rate3+se),alpha = 0.4,fill = 'black')+
    geom_vline(xintercept = 0, col ='darkgrey',linetype = 'dashed')+
    facet_wrap(~Behv)+
    theme_smo+
    theme(plot.margin = unit(c(0.1,0.2,0,1), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=20,face="bold"))+
    xlab("time from stimulus")+
    ylab("Z-scored FR")
  
  # arrange plots
  plot <- ggpubr::ggarrange(map, smo, heights = c(2, 0.8),vjust= c(0.5,0.5),
                    ncol = 1, nrow = 2, align = "v",common.legend=T,legend = 'right')
  plot
  
}


plot_overview_perID <- function(Trialless,orderby_Behv = 'A2L',by_time = F,theme_smo = theme_minimal(),save = NULL){
  plots <- unique(Trialless$id) %>% 
    map(~ {
      Datatemp <- Trialless %>% filter(str_detect(cluster_id2, as.character(.x)))
      plot <- Plot_overview(Datatemp, orderby_Behv = orderby_Behv, by_time = by_time, theme_smo = theme_smo)
      ggpubr::annotate_figure(plot, top = ggpubr::text_grob(paste("ID:", .x), 
                                                            color = "red", face = "bold", size = 14)) 
    })
  
  if (!is.null(save)) {
    if (!dir.exists(save)) dir.create(save, recursive = TRUE)
    walk2(plots, unique(Trialless$id), ~ {
      ggsave(filename = file.path(save, paste0("plot_ID_", .y, ".png")), plot = .x, width = 10, height = 8,bg='white')
    })
  }
  
  # Display each plot
  walk(plots, print)
}

add_layers <- function(Data, depthpath){
  list_of_depth <- list.files(path = depthpath,
                              pattern = "\\.csv$",
                              full.names = TRUE)
  depth <- read_csv(list_of_depth, id = 'id',col_names = F)
  depth <- depth %>% mutate(id = str_extract(depth$id, "M\\d+_\\d+")) %>% rename('surface'=X1)
  
  
  Layerd <- Data %>% separate_wider_delim(cluster_id2, delim ='_', names = c('del', 'Mouse','Day'), cols_remove = F) %>% unite('id', c(Mouse,Day)) %>%
    select(!del) %>% 
    left_join(depth) %>% replace_na(list(surface = 950)) %>% 
    mutate(dist_surf = surface - depth) %>% 
    mutate(Layer_fine = case_when(
      between(dist_surf,-500,-20) ~ 'out',
      between(dist_surf,-20,50) ~ 'L1',
      between(dist_surf,50,500) ~ 'L2/3',
      between(dist_surf,500,750) ~ 'L5A',
      between(dist_surf,750,1050) ~ 'L5B',
      between(dist_surf,1050,1210) ~ 'L6',
      .default = "WTF"
    )) %>% 
    mutate(Layer= case_when(
    between(dist_surf,-500,-20) ~ 'out',
    between(dist_surf,-20,500) ~ 'L2/3',
    between(dist_surf,500,1210) ~ 'L5',
    .default = "WTF"
  ))
  Layerd
}

plot_single_unit <- function(Combcont,cluster_id,smoothness = 0.2){
  
  Datatemp <- Combcont %>% filter(cluster_id2 == !!cluster_id) %>%  filter(event_count > 0, between(rel_time, -1, 1)) %>%  group_by(Behv) %>%
    mutate(t2 = dense_rank(Trial))
  
  if (nrow(Datatemp) == 0) {
    stop("No data found for the specified cluster_id and filtering criteria.")
  }
  
  
  # Create the PSTH
  map <- Datatemp %>% 
    ggplot() +
    geom_jitter(aes(rel_time,t2), height = 0.2,size=0.5) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0.2,0,1), 'lines'),
          strip.text.x = element_text(size = 30))+
    scale_x_discrete(labels = NULL, breaks = NULL)+
    labs(x = NULL)
  
  #density plot
  smo <- Datatemp  %>%
    ggplot() +
    geom_density(aes(rel_time), bw = smoothness) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
   theme(plot.margin = unit(c(0.1,0.2,0,1), 'lines'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))+
    xlab("time from stimulus")+
    ylab("Z-scored FR")
  
  plot <- ggpubr::ggarrange(map,smo,  heights = c(2, 0.8), vjust = c(0.5, 0.5),
                            ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend = 'right')
  
  plot
}

plot_single_unit_top <- function(Combcont,cluster_id,smoothness = 0.03){
  
  Datatemp <- Combcont %>% filter(cluster_id2 == !!cluster_id) %>%  filter(event_count > 0, between(rel_time, -1, 1)) %>%  group_by(Behv) %>%
    mutate(t2 = dense_rank(Trial))
  # Create the PSTH
  map <- Datatemp %>% 
    ggplot() +
    geom_jitter(aes(rel_time,t2), height = 0.2,size=0.5) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0.2,0,1), 'lines'),
          strip.text.x = element_text(size = 20),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank())+
    geom_vline(xintercept = 0)
  
  #density plot
  smo <- Datatemp  %>%
    ggplot() +
    geom_density(aes(rel_time), bw = smoothness) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0.2,0,1), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank())+
    xlab("time from stimulus")+
    ylab("Point Density")+
    geom_vline(xintercept = 0)
  
  plot <- ggpubr::ggarrange(map,smo,  heights = c(2, 1),
                            ncol = 1, nrow = 2, align = "v")
  plot
}

plot_single_unit_mid <- function(Combcont,cluster_id,smoothness = 0.03){
  
  Datatemp <- Combcont %>% filter(cluster_id2 == !!cluster_id) %>%  filter(event_count > 0, between(rel_time, -1, 1)) %>%  group_by(Behv) %>%
    mutate(t2 = dense_rank(Trial))
  
  
  # Create the PSTH
  map <- Datatemp %>% 
    ggplot() +
    geom_jitter(aes(rel_time,t2), height = 0.2,size=0.5) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0.2,0.2,0,0.2), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank())
  
  
  #density plot
  smo <- Datatemp  %>%
    ggplot() +
    geom_density(aes(rel_time), bw = smoothness) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0.2,0.2,0.2), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank())
  
  
  plot <- ggpubr::ggarrange(map,smo,  heights = c(2, 1),
                            ncol = 1, nrow = 2, align = "v")
  
  
  plot
}

plot_single_unit_bot<- function(Combcont,cluster_id,smoothness = 0.03){
  
  Datatemp <- Combcont %>% filter(cluster_id2 == !!cluster_id) %>%  filter(event_count > 0, between(rel_time, -1, 1)) %>%  group_by(Behv) %>%
    mutate(t2 = dense_rank(Trial))
  
  
  # Create the PSTH
  map <- Datatemp %>% 
    ggplot() +
    geom_jitter(aes(rel_time,t2), height = 0.2,size=0.5) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0.2,0.2,0,1), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank())+
    geom_vline(xintercept = 0)
  
  #density plot
  smo <- Datatemp  %>%
    ggplot() +
    geom_density(aes(rel_time), bw = smoothness) +
    facet_wrap(~Behv, scales = 'free_y')+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text=element_text(size=12),
          axis.text.y=element_blank(),
          axis.title=element_blank())+
    xlab("time from stimulus")+
    ylab("Point Density")+
    geom_vline(xintercept = 0)
  
  plot <- ggpubr::ggarrange(map,smo,  heights = c(2, 1),
                            ncol = 1, nrow = 2, align = "v")
  
  
  plot
}

Plot_all_units <- function(Data, output_dir ){
  plots <- unique(combcont$cluster_id2) %>%
  walk(~ {
    # Filter data for the current cluster_id2
    Datatemp <- combcont %>% filter(cluster_id2 == .x)

    # Create the density plot
    map <- Datatemp %>%
      filter(event_count > 0, between(rel_time, -1, 1)) %>%
      group_by(Behv) %>%
      mutate(t2 = dense_rank(Trial)) %>%
      ggplot() +
      geom_density(aes(rel_time), bw = 0.03) +
      facet_wrap(~Behv, scales = "free_y") +
      theme_minimal()

    # Create the histogram plot
    smo <- Datatemp %>%
      filter(event_count > 0, between(rel_time, -1, 1)) %>%
      group_by(Behv) %>%
      mutate(t2 = dense_rank(Trial)) %>%
      ggplot() +
      geom_jitter(aes(rel_time, t2), height = 0.2, size = 0.5) +
      facet_wrap(~Behv, scales = "free_y") +
      theme_minimal()

    # Combine the two plots vertically with ggarrange
    plot <- ggpubr::ggarrange(smo, map,
      heights = c(2, 0.8), vjust = c(0.5, 0.5),
      ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend = "right"
    )

    # Add a title with the cluster_id2
    plot <- ggpubr::annotate_figure(plot,
      top = ggpubr::text_grob(paste("ID:", .x),
        color = "red", face = "bold", size = 14
      )
    )

    # Save the plot to the output directory with a filename based on cluster_id2 without displaying
    suppressMessages(
      ggsave(
        filename = paste0(output_dir, "/", .x, "_plot.png"),
        bg = "white"
      )
    )
  })
}
