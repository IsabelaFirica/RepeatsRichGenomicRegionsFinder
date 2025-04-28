

# functie care genereaza coord random
generate_random_intervals <- function(chrom_size, seq_lengths) {
  n <- length(seq_lengths)
  intervals <- data.frame(start = integer(n), end = integer(n))
  
  # set de pozitii ocupate
  occupied_positions <- integer(0) 
  
  for (i in seq_along(seq_lengths)) {
    repeat {
      # genereaza coordonate random
      start_pos <- sample(1:(chrom_size - seq_lengths[i] + 1), 1) 
      end_pos <- start_pos + seq_lengths[i] - 1
      
      # verifica suprapunerea
      if (!any(occupied_positions %in% start_pos:end_pos)) {
        intervals[i, ] <- c(start_pos, end_pos)
        occupied_positions <- c(occupied_positions, start_pos:end_pos)
        break
      }
    }
  }
  
  return(intervals)
}

# ia fisierele
intervals_files <- list.files(pattern = "chr_limits\\.txt$")
centr_intervals_files <- list.files(pattern = "out_of_centr_limits\\.txt$")

if (length(intervals_files) > 0) {

  chrom_sizes <- read.table("chrom_sizes.txt", stringsAsFactors = FALSE)
  # ia fiecare fisier pe rand
  for (file in intervals_files) {
  
    # extrage numele cromozomului
    chr_limits <- sub("^(.*)_limits.*\\.txt$", "\\1", file)  
  
    # ia fiecare rand din fiserul cu lungimile crz
    for (i in 1:nrow(chrom_sizes)) {
    
      # ia numele crz de pe rand
      chr_fasta <- as.character(chrom_sizes[i, 1])  
    
      # verifica daca crz se potrivesc
      if (chr_limits == chr_fasta) {
      
        # citeste tabelul 
        interval_data <- read.table(file, header = FALSE)
      
        # calculeaza lungimile secventelor
        seq_lengths <- interval_data[,3] - interval_data[,2]
      
        # le ordoneaza descrescator
        seq_lengths <- sort(seq_lengths, decreasing = TRUE)
      
        # ia lungimea cromozomului actual
        chr_size <- chrom_sizes[i, 2]
      
        # generează seturile de coordonate random
        set1 <- generate_random_intervals(chr_size, seq_lengths)
        set2 <- generate_random_intervals(chr_size, seq_lengths)
        set3 <- generate_random_intervals(chr_size, seq_lengths)
      
        # adauga numele crz pe prima coloana
        set1 <- cbind(chr = chr_limits, set1)
        set2 <- cbind(chr = chr_limits, set2)
        set3 <- cbind(chr = chr_limits, set3)
      
        # creeaza numele fisierelor de output
        file1 <- paste0(chr_limits, "_random_coords_set1.txt")
        file2 <- paste0(chr_limits, "_random_coords_set2.txt")
        file3 <- paste0(chr_limits, "_random_coords_set3.txt")
      
        # salveaza coordonatele random
        write.table(set1, file = file1, quote = FALSE, col.names = FALSE, row.names = FALSE)
        write.table(set2, file = file2, quote = FALSE, col.names = FALSE, row.names = FALSE)
        write.table(set3, file = file3, quote = FALSE, col.names = FALSE, row.names = FALSE)
      }
    }
  }
}

if (length(centr_intervals_files) > 0) {

  centr_start <- read.table("centromere_coordinates.txt", stringsAsFactors = FALSE) 
  # ia fiecare fisier pe rand
  for (file in centr_intervals_files) {
  
    # extrage numele cromozomului
    chr_limits <- sub("^(.*)_out_of_centr_limits.*\\.txt$", "\\1", file)  
  
    #ia fiecare rand din fiserul cu lungimile crz
    for (i in 1:nrow(centr_start)) {
    
      #ia numele crz de pe rand
      chr_fasta <- as.character(centr_start[i, 1])  
    
      # verifica daca crz se potrivesc
      if (chr_limits == chr_fasta) {
      
        # citeste tabelul 
        interval_data <- read.table(file, header = FALSE)
      
        # calculeaza lungimile secventelor
        seq_lengths <- interval_data[,3] - interval_data[,2]
      
        # le ordoneaza descrescator
        seq_lengths <- sort(seq_lengths, decreasing = TRUE)
      
        # ia lungimea cromozomului actual
        centromere_start <- centr_start[i, 2]
      
        # generează seturile de coordonate random
        set1 <- generate_random_intervals(centromere_start, seq_lengths) 
        set2 <- generate_random_intervals(centromere_start, seq_lengths)
        set3 <- generate_random_intervals(centromere_start, seq_lengths)
      
        # adauga numele crz pe prima coloana
        set1 <- cbind(chr = chr_limits, set1)
        set2 <- cbind(chr = chr_limits, set2)
        set3 <- cbind(chr = chr_limits, set3)
      
        # creeaza numele fisierelor de output
        file1 <- paste0(chr_limits, "_random_coords_out_of_centr_set1.txt")
        file2 <- paste0(chr_limits, "_random_coords_out_of_centr_set2.txt")
        file3 <- paste0(chr_limits, "_random_coords_out_of_centr_set3.txt")
      
        # salveaza coordonatele random
        write.table(set1, file = file1, quote = FALSE, col.names = FALSE, row.names = FALSE)
        write.table(set2, file = file2, quote = FALSE, col.names = FALSE, row.names = FALSE)
        write.table(set3, file = file3, quote = FALSE, col.names = FALSE, row.names = FALSE)
      }
    }
  }
}

