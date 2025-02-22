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
  meta <- read_tsv(list_of_meta, id = 'id2') %>%  as_tibble()%>% select(id2,cluster_id,n_spikes,depth,group,fr)
  meta <- meta %>% mutate(id = str_extract(meta$id2, "M\\d+_\\d+"))
  
  Data <- Data %>% select(cluster_id:trial,id) %>% left_join(meta %>% select(!id2),by = join_by(id,cluster_id))
  # adds more data such as the behaviour,depth or Hit/miss based on which Trialtypes where present 
  Data <- Data %>% mutate(Trial = sprintf('%03d',trial)) %>% filter(group!='noise')
  Data <- Data %>% unite('Trial2',c(id,Trial),remove=F) %>% group_by(Trial2) %>%  mutate(Behv = case_when(
    any(Trialtype %in% c('Audio', 'W2T_Audio')) ~ "W2T",
    any(Trialtype %in% c('Airpuff', 'A2L_Audio')) ~ "A2L",
    any(Trialtype %in% c('PC_Audio', 'Airpuff2')) ~ "PC", 
    .default = "WTF")) %>% mutate(succ = ifelse('HIT' %in% Trialtype,'Hit','Miss')) %>% ungroup() %>%  unite('cluster_id2',c(cluster_id,id),remove=F) 
  Data
}

load_opto <- function(datapath,metapath){
  
  # loads individual data (from raw_ks_to_df) and meta data (from phy output) and combines them into a df
  # takes names from filename. need to be M[Mouse_id]_[recording session] e.g. M1_2
  list_of_data <- list.files(path = datapath,
                             pattern = "\\.csv$",
                             full.names = TRUE)
  Data <- read_csv(list_of_data, id = 'group') 
  Data <- Data %>% mutate(id = str_extract(Data$group, "O\\d+_[A-Za-z0-9]"))
  
  list_of_meta <- list.files(path = metapath,
                             pattern = "\\.tsv$",
                             full.names = TRUE)
  meta <- read_tsv(list_of_meta, id = 'id2') %>%  as_tibble()%>% select(id2,cluster_id,n_spikes,depth,group,fr)
  meta <- meta %>% mutate(id = str_extract(meta$id2, "O\\d+_[A-Za-z0-9]"))
  
  Data
  
  Data <- Data %>% select(cluster_id:id) %>% left_join(meta %>% select(!id2),by = join_by(id,cluster_id))
  # adds more data such as the behaviour,depth or Hit/miss based on which Trialtypes where present 
  Data <- Data %>% mutate(Trial = sprintf('%03d',trial)) %>% filter(group!='noise')
  Data <- Data %>% unite('Trial2',c(id,Trial),remove=F) %>% group_by(Trial2) %>% mutate(succ = ifelse('HIT' %in% Trialtype,'Hit','Miss')) %>% ungroup() %>%  unite('cluster_id2',c(cluster_id,id),remove=F) 
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




avg_trials <- function(Data,bin_size=0.01){
  Trialless <- combcont %>% group_by(Behv,rel_time, cluster_id2)  %>% 
    mutate(rate = (sum(event_count)/n_trials), rate = rate/bin_size) %>% ungroup() %>% 
    distinct(Behv,rel_time, cluster_id2,id,rate,n_spikes,depth,group,fr,n_trials) %>% group_by(Behv,cluster_id2) %>% 
    mutate(Zscore = (rate - mean(rate))/sd(rate))
  
}


load_data_memoryLimited <- function(datapath, metapath, output_folder, chunk_size = 1) {
  # Ensure the output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  # Helper function to process and save a chunk of CSV files
  process_and_save_chunk <- function(file_chunk, meta, output_folder, chunk_index) {
    # Load and process all files in the chunk
    chunk_data <- file_chunk %>%
      map_dfr(~ read_csv(.x, id = 'group')) %>%
      mutate(id = str_extract(group, "M\\d+_\\d+"))
    
    # Add meta data
    chunk_data <- chunk_data %>%
      select(cluster_id:trial, id) %>%
      left_join(meta %>% select(-id2), by = join_by(id, cluster_id)) %>%
      mutate(
        Trial = sprintf('%03d', trial)) %>%  unite('Trial2', c(id, Trial), remove = FALSE) %>% group_by(Trial2) %>% 
      mutate(
        Behv = case_when(
          any(Trialtype %in% c('Audio', 'W2T_Audio')) ~ "W2T",
          any(Trialtype %in% c('Airpuff', 'A2L_Audio')) ~ "A2L",
          any(Trialtype %in% c('PC_Audio', 'Airpuff2')) ~ "PC",
          .default = "WTF"
        ),
        succ = ifelse('HIT' %in% Trialtype, 'Hit', 'Miss')
      ) %>% ungroup() %>%
      unite('cluster_id2', c(cluster_id, id), remove = FALSE)
    
    # Save chunk to output folder
    chunk_name <- paste0(output_folder, "/chunk_", chunk_index, ".csv")
    write_csv(chunk_data, chunk_name)
    message(paste("Processed and saved chunk:", chunk_name))
    
    rm(chunk_data)
    gc() # Trig
  }
  
  # Load and prepare meta data
  list_of_meta <- list.files(path = metapath, pattern = "\\.tsv$", full.names = TRUE)
  meta <- map_dfr(list_of_meta, read_tsv, id = 'id2') %>%
    as_tibble() %>%
    select(id2, cluster_id, n_spikes, depth, group, fr) %>%
    mutate(id = str_extract(id2, "M\\d+_\\d+"))
  
  # Get list of CSV files
  list_of_data <- list.files(path = datapath, pattern = "\\.csv$", full.names = TRUE)
  
  # Process files in chunks
  chunk_index <- 1
  for (i in seq(1, length(list_of_data), by = chunk_size)) {
    file_chunk <- list_of_data[i:min(i + chunk_size - 1, length(list_of_data))]
    process_and_save_chunk(file_chunk, meta, output_folder, chunk_index)
    chunk_index <- chunk_index + 1
  }
  
  # Combine all processed chunks into a single data frame
  processed_files <- list.files(output_folder, pattern = "\\.csv$", full.names = TRUE)
  final_data <- map_dfr(processed_files, read_csv)
  
  message("All files processed and combined.")
  final_data
}


  