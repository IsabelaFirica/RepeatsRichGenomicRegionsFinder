

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)


# ia toate fisierele mapped_contigs.csv
files <- list.files(pattern = "^mapped_contigs.*\\.csv$")

data_list <- lapply(files, read.csv)


# nr minim de hit-uri de pe o coloana
min_count <- 10 

# distanta maxima dintre 2 col luate in acelasi grup
group_gap <- 5 


# ia pe rand fiecare tabel
for (file in files) {
  data <- read.csv(file)

  # salveaza numele chr din numele fisierului
  chr <- sub("^mapped_contigs_(.*)\\.csv$", "\\1", file)
  
  # fisierele de output
  output_file <- paste0(chr, "_columns_selected.txt")
  limits_file <- paste0(chr, "_limits.txt")
  output_excel <- paste0(chr, "_colored.xlsx")


  colnames_data <- colnames(data) 
  num_cols <- ncol(data) 


  # liste goale
  results <- list()
  limits <- list()
  columns <- list() 
  boundaries <- list() 

  i <- 2 
  while (i <= num_cols) { 
    count <- sum(data[[i]] == "______", na.rm = TRUE) 
    
    if (count >= min_count) { # daca nr de "______" e mai mare ca min ales
      current_group <- c(i)  # adauga i la grup
      j <- i + 1  # verifica urmatoarea coloana
      gap_count <- 0  # distanta dintre col e 0
    
      while (j <= num_cols) { 
        next_count <- sum(data[[j]] == "______", na.rm = TRUE) 
        
        if (next_count >= min_count) {  # daca nr e mai mare ca min
          current_group <- c(current_group, j)  # adauga j la grup
          gap_count <- 0  # distanta ramane 0
          
        } else {
          gap_count <- gap_count + 1 # daca nu, gap_count creste cu 1
          
          if (gap_count > group_gap) {  # daca dist e mai mare ca cea min, se opreste
            break 
          }
        }
        j <- j + 1 
      }
    
      # lista cu toate col. care indeplinesc conditia (indicii lor)
      columns <- append(columns, current_group)
    
      # daca grupul are mai mult de o col, ia limitele la 2 col la stanga/dreapta
      if (length(current_group) > 1) {
        left_boundary <- max(2, min(current_group) - 2)
        right_boundary <- min(num_cols, max(current_group) + 2)
        
      } else {
        # daca e o singura col, limitele sunt prima col la stanga/dreapta
        left_boundary <- max(2, current_group[1] - 1)
        right_boundary <- min(num_cols, current_group[1] + 1)
      }
    
      # lista cu toate col. limita (indici)
      boundaries <- append(boundaries, left_boundary)
      boundaries <- append(boundaries, right_boundary)
    
      # adauga grupurile si limitele in lista (numele coloanelor)
      results <- append(results, list(
        list(
          "Columns Meeting Condition" = colnames_data[current_group],
          "Left Boundary Column" = colnames_data[left_boundary],
          "Right Boundary Column" = colnames_data[right_boundary]
        )
      ))
    
      # lista doar cu coordonatele limita
      left_limit <- gsub(".*\\.(\\d+)\\..*", "\\1", colnames_data[left_boundary])
      right_limit <- gsub(".*\\.(\\d+)\\.(\\d+)", "\\2", colnames_data[right_boundary])
    
      limits <- c(limits, paste(chr, left_limit, right_limit))
    
      i <- max(current_group) + 1
    } else {
      i <- i + 1
    }
  }


  output <- capture.output({
    for (group in results) {
      cat("Group of columns:\n")
      cat("  Columns meeting condition: ", paste(group[["Columns Meeting Condition"]], 
          collapse = ", "), "\n")
      cat("  Left boundary column: ", group[["Left Boundary Column"]], "\n")
      cat("  Right boundary column: ", group[["Right Boundary Column"]], "\n")
    }
  })
  
  #################################
  
  if (file.exists("centromere_coordinates.txt")) {
    centromere_coord <- read.table("centromere_coordinates.txt")
    colnames(centromere_coord) <- c("chr", "coord")
   
    centr_coord <- as.numeric(centromere_coord[centromere_coord$chr == chr, "coord"])
    
    before_centromere <- c()
    in_centromere <- c()
    
    for (limit in limits) {
      fields <- strsplit(limit, " ")[[1]]
      this_chr <- fields[1]
      start <- as.numeric(fields[2])
      end <- as.numeric(fields[3])
      
      if (this_chr == chr) {
        if (start < centr_coord) {
          before_centromere <- c(before_centromere, limit)
        } else {
          in_centromere <- c(in_centromere, limit)
        } 
      }
    }
    
    out_limits <- paste0(chr, "_out_of_centr_limits.txt")
    in_limits <- paste0(chr, "_in_centr_limits.txt")
    
    writeLines(unlist(before_centromere), con = out_limits)
    writeLines(unlist(in_centromere), con = in_limits)
    
    }
  else {
    writeLines(unlist(limits), con = limits_file)
    }

  ###############################
  
  writeLines(output, con = output_file)


  # face tabelul cu coloanele colorate

  wb <- createWorkbook()

  addWorksheet(wb, "Sheet1")

  writeData(wb, "Sheet1", data)

  color1 <- createStyle(fgFill = "lightblue1")
  color2 <- createStyle(fgFill = "lightblue3")

  num_rows <- nrow(data)

  # scoate coloana 2 (prima col.) din boundaries
  boundaries <- lapply(boundaries, function(x) setdiff(x, 2))
  boundaries <- boundaries[lengths(boundaries) > 0]

  addStyle(wb, "Sheet1", color1, rows = 2:num_rows, cols = columns, gridExpand = TRUE)
  addStyle(wb, "Sheet1", color2, rows = 2:num_rows, cols = boundaries, gridExpand = TRUE)

  saveWorkbook(wb, output_excel, overwrite = TRUE)

}


