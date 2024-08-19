
mat_to_tibble <- function(matpath){
  list_files <- list.files(path = matpath,full.names = TRUE)
  datalist = lapply(list_files, readMat)
  
  
  flist = list()
  times <- list()
  
  for (x in 1:length(datalist))
  {
    i=1
    data <- datalist[[x]]
    time <- hms::as_hms(data[["SessionData"]][[2]][[6]])
    dat <- data[[1]][[4]][[1]]
    df=tibble()
    for (i in 1:dim(dat)[2])
    {
      zz <- dat %>% flatten()  %>% pluck(i,1) %>%  enframe() %>%   unnest_wider(value,names_sep="_",names_repair='unique')
      df <- rbind(df,tibble(name = zz$name, 
                            start=zz$value_1[,1], 
                            end=zz$value_1[,2],Trial=i,
                            id= str_extract( list_files[x],"(?<=/)[^/]+(?=\\.mat)"),
                            times = time))
      flist[[x]] = df   
    }
    times[[x]]= time
  }
  
  
  
  df <- bind_rows(flist, .id = "column_label") %>% na.omit() 
  Data2 <- df %>% separate_wider_delim(id,delim ='_',names = c('mouse','day')) %>%  unite('Trial2',c(mouse,day,Trial),remove=F) %>% group_by(Trial2)%>%
    mutate(Behv = case_when(
      'Audio' %in% name ~ "W2T",
      'Airpuff' %in% name ~ "A2L",
      'Audio2Air' %in% name ~"PC",
      .default = "WTF")) %>% mutate(succ = ifelse('HIT' %in% name,'Hit','Miss')) %>% ungroup() %>% unite('id',mouse,day)
  Data2
}

load_creas_times <- function(timespath){
  list2 <- list.files(path = timespath,full.names = TRUE)
  creas <- read_csv(list2,id = "id")
  
  creas_times <- creas %>% mutate(id = str_extract(id, "(?<=/)[^/]+(?=\\.csv)"), video = as.numeric(str_extract(FilePath, "\\d+(?=\\.mp4v$)")))  %>% select(id,CreationTime,video)
  creas_times
}

match_within_group <- function(crea_group, hit_group) {
  matched_hits <- integer(0)
  
  find_closest_unmatched <- function(crea_time, hit_times, matched_hits) {
    unmatched_hits <- setdiff(seq_along(hit_times), matched_hits)
    closest_idx <- unmatched_hits[which.min(abs(hit_times[unmatched_hits] - crea_time))]
    return(closest_idx)
  }
  
  closest_hit_indices <- integer(nrow(crea_group))
  
  for (i in seq_len(nrow(crea_group))) {
    crea_time <- crea_group$CreationTime[i]
    closest_idx <- find_closest_unmatched(crea_time, hit_group$timediff, matched_hits)
    closest_hit_indices[i] <- closest_idx
    matched_hits <- append(matched_hits, closest_idx)
  }
  
  crea_group <- crea_group %>%
    mutate(closest_hit_idx = closest_hit_indices,
           closest_hit_time = hit_group$timediff[closest_hit_idx])
  
  return(crea_group)
}

add_times <- function(data){
  data2 <- data %>% group_by(id) %>% mutate(len = end - start) %>% mutate(cumsec = round(cumsum(len),2), cum2 = hms::as_hms(cumsec), timediff = hms::as_hms(cum2 + times)) %>% ungroup()
  data2
}