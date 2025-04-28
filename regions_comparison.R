

suppressPackageStartupMessages(library(dplyr))

# ia toate fisierel cu rezultate blast filtrate
blast_results_filtered <- list.files(pattern = "_filtr\\.txt$")

out_of_centr_files <- list.files(pattern = "out_of_centr.*_filtr\\.txt$")
in_centr_files <- list.files(pattern = "in_centr.*_filtr\\.txt$")


# functie care extrage chr din numele fisierelor
get_chr_name <- function(file_name) {
  chr <- sub("^(.*)_blast.*_filtr\\.txt$", "\\1", file_name)  
}

if (length(out_of_centr_files) == 0) {
  # grupeaza fisierele dupa chr
  blast_results_grouped <- split(blast_results_filtered, 
                               sapply(blast_results_filtered, get_chr_name))



  # loop prin grupurile de fișiere (per cromozom)
  for (chr_name in names(blast_results_grouped)) {
  
    #tabel gol pt rezultate
    final_table <- data.frame()
  
    # extrage fișierele pentru acest cromozom
    chr_files <- blast_results_grouped[[chr_name]]
  
    for (file in chr_files){
    
      # ia din nume tipul aliniamentului (ref, cont, random)
      type <- sub(".*blast_(.*)_filtr\\.txt$", "\\1", file)
    
      # citeste fisierul
      transp_table <- read.table(file, header = FALSE)
    
      # numara de cate ori apare fiecare tranpozon
      transp_table <- transp_table %>% count(V1)
      colnames(transp_table) <- c("transposon", type)
    
      # adauga in tabel
      if (nrow(final_table) == 0) {
        final_table <- transp_table
      } else {
        final_table <- full_join(final_table, transp_table, by = "transposon")
      }
    }

  
    final_table <- final_table %>% arrange(transposon)
  
    final_table[is.na(final_table)] <- 0
  
    final_table <- final_table %>%
      select(transposon, ref, cont, random_set1, random_set2, random_set3)
  
    output_file <- paste0(chr_name, "_transposons_comparison.csv")
  
    write.csv(final_table, file = output_file, quote = FALSE, row.names = FALSE)
  
  }
}

if (length(out_of_centr_files) > 0) {
  # grupeaza fisierele dupa chr
  blast_results_grouped <- split(out_of_centr_files, 
                               sapply(out_of_centr_files, get_chr_name))



  # loop prin grupurile de fișiere (per cromozom)
  for (chr_name in names(blast_results_grouped)) {
  
    #tabel gol pt rezultate
    final_table <- data.frame()
  
    # extrage fișierele pentru acest cromozom
    chr_files <- blast_results_grouped[[chr_name]]
  
    for (file in chr_files){
    
      # ia din nume tipul aliniamentului (ref, cont, random)
      type <- sub(".*blast_(.*)_filtr\\.txt$", "\\1", file)
    
      # citeste fisierul
      transp_table <- read.table(file, header = FALSE)
    
      # numara de cate ori apare fiecare tranpozon
      transp_table <- transp_table %>% count(V1)
      colnames(transp_table) <- c("transposon", type)
    
      # adauga in tabel
      if (nrow(final_table) == 0) {
        final_table <- transp_table
      } else {
        final_table <- full_join(final_table, transp_table, by = "transposon")
      }
    }

  
    final_table <- final_table %>% arrange(transposon)
  
    final_table[is.na(final_table)] <- 0
  
    final_table <- final_table %>%
      select(transposon, out_of_centr_ref, out_of_centr_cont, out_of_centr_random_set1, 
      out_of_centr_random_set2, out_of_centr_random_set3) %>%
    rename(
      reference = out_of_centr_ref,
      contigs = out_of_centr_cont,
      random_set1 = out_of_centr_random_set1,
      random_set2 = out_of_centr_random_set2,
      random_set3 = out_of_centr_random_set3
    )
  
    output_file <- paste0(chr_name, "_out_of_centr_transposons_comparison.csv")
  
    write.csv(final_table, file = output_file, quote = FALSE, row.names = FALSE)
  
  }
}

if (length(in_centr_files) > 0) {
  # grupeaza fisierele dupa chr
  blast_results_grouped <- split(in_centr_files, 
                               sapply(in_centr_files, get_chr_name))

  for (chr_name in names(blast_results_grouped)) {
  
    #tabel gol pt rezultate
    final_table <- data.frame()
  
    # extrage fișierele pentru acest cromozom
    chr_files <- blast_results_grouped[[chr_name]]
  
    for (file in chr_files){
    
      # ia din nume tipul aliniamentului (ref, cont, random)
      type <- sub(".*blast_(.*)_filtr\\.txt$", "\\1", file)
    
      # citeste fisierul
      transp_table <- read.table(file, header = FALSE)
    
      # numara de cate ori apare fiecare tranpozon
      transp_table <- transp_table %>% count(V1)
      colnames(transp_table) <- c("transposon", type)
    
      # adauga in tabel
      if (nrow(final_table) == 0) {
        final_table <- transp_table
      } else {
        final_table <- full_join(final_table, transp_table, by = "transposon")
      }
    }

  
    final_table <- final_table %>% arrange(transposon)
  
    final_table[is.na(final_table)] <- 0
  
    final_table <- final_table %>%
      select(transposon, in_centr_ref, in_centr_cont) %>%
    rename(
      reference = in_centr_ref,
      contigs = in_centr_cont
    )
  
    output_file <- paste0(chr_name, "_in_centr_transposons_comparison.csv")
  
    write.csv(final_table, file = output_file, quote = FALSE, row.names = FALSE)
  
  }
}

