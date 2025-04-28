

packages <- c("dplyr", "tidyr")

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(tidyr)
suppressPackageStartupMessages(library(dplyr))


# ia fisierele paf cu alinierile 
centr_files <- list.files(pattern = "_alignment_(in_centr|out_of_centr)_seq_to_cont\\.paf$")
normal_files <- list.files(pattern = "_alignment_seq_to_cont\\.paf$")
ref_alignment_files <- list.files(pattern = "_alignment_cont_to_ref\\.paf$")


# Function to extract chr name
get_chr_name <- function(file_name) {
  chr <- sub("^(.*)_alignment_?.*?_?seq_to_cont\\.paf$", "\\1", file_name)
  return(chr)
}

get_chr_from_ref_file <- function(file_name) {
  chr <- sub("^(.*)_alignment_cont_to_ref\\.paf$", "\\1", file_name)
  return(chr)
}


# Define mode and matching function
if (length(centr_files) > 0) {
  suffixes <- c("in_centr", "out_of_centr")
} else {
  suffixes <- c("")
}


# ia pe rand fisierele paf cu alinierea secv la cont

for (suffix in suffixes) {
  if (suffix == "") {
    reg_alignment_files <- list.files(pattern = "_alignment_seq_to_cont\\.paf$")
  } else {
    reg_alignment_files <- list.files(pattern = paste0("_alignment_", suffix, "_seq_to_cont\\.paf$"))
  }
  
  for (reg_alignment in reg_alignment_files) {
    reg_chr <- get_chr_name(reg_alignment)  # ia chr din numele fisierului
    
    # ia pe rand fisierele paf cu alinierea cont la ref
    for (ref_alignment in ref_alignment_files) {
      ref_chr <- get_chr_from_ref_file(ref_alignment)  # ia chr din numele fisierului
      
      # continua daca numele chr se potrivesc 
      if (reg_chr == ref_chr) {
        
        chr <- reg_chr  # salveaza numele chr pt fiserele output
        
        # Dynamic output filenames
        output_tag <- if (suffix == "") "" else paste0(suffix, "_")
        output_file_csv <- paste0(chr, "_", output_tag, "contig_coordinates_table.csv")
        coordinates_file <- paste0(chr, "_", output_tag, "contig_coordinates.txt")
        
        ##############################
        
        # ia doar coloanele necesare din fisierele paf
        
        reg_df <- read.table(reg_alignment)
        
        reg_df <- reg_df %>% select(V1, V3, V4, V6, V8, V9)
        
        reg_df <- separate(reg_df, V1, into = 
                             c("ref_name", "ref_start", "ref_end"), sep = "[:\\-]", convert = TRUE)
        
        reg_df <- reg_df %>% rename(contig = V6, contig_start = V8, 
                                    contig_end = V9, algn_start = V3, algn_end = V4)
        
        
        
        ref_df <- read.table(ref_alignment)
        
        ref_df <- ref_df %>% select(V1, V8, V9)
        
        ref_df <- ref_df %>%
          rename(contig = V1, ref_start = V8, ref_end = V9)
        
        
        ##############################
        
        # tabel gol pt rezultate
        results_table <- data.frame(contig = c(), left_coord = c(), right_coord = c(),
                                    contig_length = c(), ref_start = c(),
                                    ref_end = c(), ref_length = c())
        
        # ia pe rand fiecare rand din cele doua tabele
        for (i in 1:nrow(reg_df)){
          for (j in 1:nrow(ref_df)){
            
            reg <- reg_df[i,]
            ref <- ref_df[j,]
            
            # daca alinierile sunt pt acelasi contig
            if (ref$contig == reg$contig){
              
              # verifica daca secv din ref care s-a aliniat la contig se afla in
              # intervalul in care s-a aliniat contigul la crz
              if (reg$ref_start >= ref$ref_start && reg$ref_start <= ref$ref_end ||
                  reg$ref_end >= ref$ref_start && reg$ref_end <= ref$ref_end){
                
                # lungimea secventelor din contiguri si din ref
                contig_length <- reg$contig_end - reg$contig_start
                ref_length <- reg$ref_end - reg$ref_start
                
                # adauga in tabel coordonatele, lungimea secventelor etc.
                results_table <- dplyr::bind_rows(results_table, data.frame(
                  contig = reg$contig,
                  left_coord = reg$contig_start,
                  right_coord = reg$contig_end,
                  contig_length = contig_length,
                  ref_start = reg$ref_start, 
                  ref_end = reg$ref_end,
                  ref_length = ref_length))
              }
            }
          }
        }
        
        # pt fiecare rezultat pt un anumit contig si o anume secventa din ref
        # pastreaza alinierea cea mai buna 
        results_table <- results_table %>%
          group_by(contig, ref_start) %>%
          slice(1) %>% # Keep only the first row within each group
          ungroup()
        
        # ordoneaza tabelul
        results_table <- results_table %>% arrange(ref_start)
        
        # scrie fisierele de output
        write.csv(results_table, file = output_file_csv, row.names = FALSE)
        
        coordinates <- results_table %>% select(contig, left_coord, right_coord)
        
        write.table(coordinates, file = coordinates_file, quote = FALSE, 
                    col.names = FALSE, row.names = FALSE)
        
      }
    }
  }
}


