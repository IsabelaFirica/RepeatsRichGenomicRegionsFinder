#!/bin/bash

required_packages=(seqtk minimap2 r-base ncbi-blast+)

for pkg in "${required_packages[@]}"; do
  if ! dpkg -s "$pkg" &> /dev/null; then
    sudo apt update > /dev/null 2>&1
    sudo apt install -y "$pkg" > /dev/null 2>&1 && echo "Installed: $pkg"
  fi
done

# verifica daca au fost date 4 argumente
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <dscaff_folder> <chr_fasta_folder> <transposons_fasta> <transposons_lengths"
  exit 1
fi

# salveaza argumentele date
dscaff_folder="$1"
chr_fasta_folder="$2"
transp_fasta="$3"
transp_lengths="$4"

# impartit in intervale out of/inside centromere region
echo " "
echo "Analyze the the intervals in/out of centromere region separately? (yes/no)"

read -r answer

if [[ "$answer" == "yes" ]]; then
  echo " "
  echo "Introduce the start coordinate of the centromere region for each chromosome:"
  output_file="centromere_coordinates.txt"
  > "$output_file"  

  exec 3< "$MAIN_DIR/chr_list.txt"  

  while read -r chr <&3; do
    while true; do
      echo " "
      read -e -r -p "Enter coordinate for $chr: " coordinate

      if [[ "$coordinate" =~ ^[0-9]+$ ]]; then
        echo "$chr $coordinate" >> "$output_file"
        echo " "
        echo "Saved: $chr $coordinate"
        break
      else
        echo " "
        echo "Invalid input. Please enter a numeric value for $chr."
      fi
    done
  done

elif [[ "$answer" == "no" ]]; then
  echo " "
  echo "The analysis will be done for the whole chromosome."
else
  echo " "
  echo "Invalid response. Please answer yes or no."
fi

echo " "
echo "getting coordinates from reference..."
echo " "

# folderul principal
MAIN_DIR="$PWD"

# folder pt rezultatele pt referinta
ref_coord="$MAIN_DIR/reference_coordinates"
mkdir -p $ref_coord

#intra in folderul de dscaff
cd "$dscaff_folder"

# gaseste fisierele mapped_contigs
find . -type f -path "*/chromosome/mapped_contigs.csv" | while read file; do
    # ia numele crz din numele folderului in care se afla mapped_contigs
    chr=$(basename "$(dirname "$(dirname "$file")")")
    
    # salveaza numele crz intr-un fisier text
    echo "$chr" >> "$MAIN_DIR/chr_list.txt"
    
    # aduga crz in numele fisierelor mapped_contigs
    new_name="mapped_contigs_${chr}.csv"
    
    # copiaza fisierele mapped_contigs in folderul nou
    cp "$file" "$ref_coord/$new_name"
done

sed -i 's/\r$//' "$MAIN_DIR/chr_list.txt"


# duce scriptul R in folderul nou
cd $MAIN_DIR
cp ref_coord.R $ref_coord
cd $ref_coord


# ruleaza scriptul R
Rscript ref_coord.R 


########################################################
#######################################################

echo " "
echo "extracting sequences from reference..."
echo " "

# ia fisierele fasta ale crz
cd $MAIN_DIR/$chr_fasta_folder
cp *.fasta $ref_coord

cd $ref_coord 

# aduga "_chr" la numele fisierelor (ca sa fie toate numite la fel)
for file in *.fasta; do
  base="${file%.fasta}"             
  mv "$file" "${base}_chr.fasta"   
done

# schimba numele headerelor sa contina doar numele crz 
for file in *.fasta; do
  # salveaza numele fisierului
  base_name=$(basename "$file" .fasta)

  # fisier temporar in care sa salveze rezultatele
  temp_file=$(mktemp)

  # schimba headerul
  awk -v header=">$base_name" 'BEGIN {header_set=0} 
       /^>/ { 
           if (header_set == 0) { 
               print header; 
               header_set=1; 
           }
       } 
       !/^>/ {print}' "$file" > "$temp_file"
  
  mv "$temp_file" "$file"
done

#######################################################

# extrage secventele din ref pe baza coordonatelor 

if [[ "$answer" == "no" ]]; then

  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}.fasta"            
    coords="${chr}_limits.txt"      
    out="${chr}_extracted_sequences.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do

      # verifica daca headerul se potriveste cu crz actual
      if [[ "$header" == "$chr" ]]; then
      
        seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
      fi
    done < "$coords"
  done < "$MAIN_DIR/chr_list.txt"

fi


if [[ "$answer" == "yes" ]]; then

  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}.fasta"            
    coords="${chr}_in_centr_limits.txt"      
    out="${chr}_in_centr_extr_seq.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do

      # verifica daca headerul se potriveste cu crz actual
      if [[ "$header" == "$chr" ]]; then
      
        seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
      fi
    done < "$coords"
  done < "$MAIN_DIR/chr_list.txt"

  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}.fasta"            
    coords="${chr}_out_of_centr_limits.txt"      
    out="${chr}_out_of_centr_extr_seq.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do

      # verifica daca headerul se potriveste cu crz actual
      if [[ "$header" == "$chr" ]]; then
      
        seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
      fi
    done < "$coords"
  done < "$MAIN_DIR/chr_list.txt"

fi

rm *.csv
rm ref_coord.R


#######################################################
#######################################################

echo " "
echo "getting coordinates from contigs..."
echo " "

cd $MAIN_DIR

# folder nou pentru secventele din contiguri
mkdir contigs_coordinates
cont_coord="$MAIN_DIR/contigs_coordinates"

# ia fisiere fasta cu contigurile din folderul dscaff
find "$dscaff_folder" -type f -name "*assembly.fasta" -exec cp {} "$cont_coord" \;

# ia fisierele necesare din ref_coord
cd $ref_coord
cp *sequences.fasta $cont_coord 2>/dev/null
cp *out_of_centr_extr_seq.fasta $cont_coord 2>/dev/null
cp *in_centr_extr_seq.fasta $cont_coord 2>/dev/null
cp *chr.fasta $cont_coord

#rm *_chr.fasta

cd $cont_coord

#######################################################

# schimba headerele fisierelor fasta cu contigurile, sa aiba doar numele contigului
for fasta_file in "$cont_coord"/*assembly.fasta; do
  
  # pastreaza al doilea cuvant din header
  awk '/^>/ {print ">" $2; next} {print}' "$fasta_file" > temp.fasta && mv temp.fasta "$fasta_file"
  
done

#######################################################

# aliniere cu minimap2 intre contiguri si crz
for chr_fasta in "$cont_coord"/*chr.fasta; do
  # ia numele crz fiecarui fisier fasta
  chr_name=$(basename "$chr_fasta" .fasta)

  # gaseste fisierul corespunzator cu contiguri
   assembly_fasta=$(ls "$cont_coord" | grep "^$chr_name.*assembly.fasta$")

  # daca exista, le aliniaza
  if [ -f "$assembly_fasta" ]; then
    
    output_paf="$cont_coord/${chr_name}_alignment_cont_to_ref.paf"

    minimap2 -x asm5 --secondary=no "$chr_fasta" "$assembly_fasta" > "$output_paf" 2> /dev/null
  fi
done

#######################################################

# aliniere cu minimap2 intre secventele extrase deja din ref si contiguri

if [[ "$answer" == "no" ]]; then
  for seq_fasta in "$cont_coord"/*sequences.fasta; do
    # ia numele crz fiecarui fisier fasta cu secvente
    chr_name=$(basename "$seq_fasta" _extracted_sequences.fasta)

    # gaseste fisierul corespunzator cu contiguri
     assembly_fasta=$(ls "$cont_coord" | grep "^$chr_name.*assembly.fasta$")

    # daca exista, le aliniaza
    if [ -f "$assembly_fasta" ]; then
    
      output_paf="$cont_coord/${chr_name}_alignment_seq_to_cont.paf"

      minimap2 -x asm5 --secondary=no "$assembly_fasta" "$seq_fasta" > "$output_paf" 2> /dev/null
    fi
  done
fi

if [[ "$answer" == "yes" ]]; then
  for seq_fasta in "$cont_coord"/*in_centr_extr_seq.fasta; do
    # ia numele crz fiecarui fisier fasta cu secvente
    chr_name=$(basename "$seq_fasta" _in_centr_extr_seq.fasta)

    # gaseste fisierul corespunzator cu contiguri
    assembly_fasta=$(ls "$cont_coord" | grep "^$chr_name.*assembly.fasta$")

    # daca exista, le aliniaza
    if [ -f "$assembly_fasta" ]; then
    
      output_paf="$cont_coord/${chr_name}_alignment_in_centr_seq_to_cont.paf"

      minimap2 -x asm5 --secondary=no "$assembly_fasta" "$seq_fasta" > "$output_paf" 2> /dev/null
    fi
  done

  for seq_fasta in "$cont_coord"/*out_of_centr_extr_seq.fasta; do
    # ia numele crz fiecarui fisier fasta cu secvente
    chr_name=$(basename "$seq_fasta" _out_of_centr_extr_seq.fasta)

    # gaseste fisierul corespunzator cu contiguri
    assembly_fasta=$(ls "$cont_coord" | grep "^$chr_name.*assembly.fasta$")

    # daca exista, le aliniaza
    if [ -f "$assembly_fasta" ]; then
    
      output_paf="$cont_coord/${chr_name}_alignment_out_of_centr_seq_to_cont.paf"

      minimap2 -x asm5 --secondary=no "$assembly_fasta" "$seq_fasta" > "$output_paf" 2> /dev/null
    fi
  done
fi


# copiaza scriptul R in folderul actual
cd $MAIN_DIR
cp contig_coord_new.R $cont_coord
cd $cont_coord

# ruleaza scriptul R 
Rscript contig_coord_new.R


echo " "
echo "extracting sequences from contigs..."
echo " "

# extrage secventele din contiguri pe baza coordonatelor 
if [[ "$answer" == "no" ]]; then
  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}_chromosome_dScaff_assembly.fasta"            
    coords="${chr}_contig_coordinates.txt"      
    out="${chr}_contigs_sequences.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do
      
      seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
    
    done < "$coords"

  done < "$MAIN_DIR/chr_list.txt"
fi

if [[ "$answer" == "yes" ]]; then
  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}_chromosome_dScaff_assembly.fasta"            
    coords="${chr}_in_centr_contig_coordinates.txt"      
    out="${chr}_in_centr_contigs_sequences.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do
      
      seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
    
    done < "$coords"

  done < "$MAIN_DIR/chr_list.txt"

  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}_chromosome_dScaff_assembly.fasta"            
    coords="${chr}_out_of_centr_contig_coordinates.txt"      
    out="${chr}_out_of_centr_contigs_sequences.fasta"

    > "$out"

    # extrage secventele 
    while IFS=' ' read -r header start end; do
      
      seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
    
    done < "$coords"

  done < "$MAIN_DIR/chr_list.txt"
fi


#######################################################
#######################################################

echo " "
echo "getting random coordinates from ref..."
echo " "

transp_comp="$MAIN_DIR/transposons_comparison"
mkdir -p $transp_comp

cd $MAIN_DIR

cp $transp_fasta $transp_lengths random_coords.R transposons_comparison.R $transp_comp

cd $ref_coord
cp *limits.txt $transp_comp
cp *chr.fasta $transp_comp
cp centromere_coordinates.txt $transp_comp 2>/dev/null

cd $transp_comp

# Numărăm nucleotidele din fiecare fișier FASTA

for fasta in *.fasta; do
    base_name=$(basename "$fasta" .fasta)
    count=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)
    echo "$base_name $count" >> chrom_sizes.txt
done


Rscript random_coords.R 

#######################################################

echo " "
echo "extracting random sequences from ref..."
echo " "

if [[ "$answer" == "no" ]]; then
  # extragere secvente 
  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}.fasta"
  
    for set in 1 2 3; do
      coords="${chr}_random_coords_set${set}.txt"
      out="${chr}_extr_seq_random_set${set}.fasta"

      > "$out"

      # extrage secventele 
      while IFS=' ' read -r header start end; do

        # verifica daca headerul se potriveste cu crz actual
        if [[ "$header" == "$chr" ]]; then
      
          seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
        fi
      done < "$coords"
  
    done

  done < "$MAIN_DIR/chr_list.txt"
fi

if [[ "$answer" == "yes" ]]; then
  # extragere secvente 
  while read -r chr; do
    # ia pe rand fisierele pt fiecare crz
    fasta="${chr}.fasta"
  
    for set in 1 2 3; do
      coords="${chr}_random_coords_out_of_centr_set${set}.txt"
      out="${chr}_extr_seq_out_of_centr_random_set${set}.fasta"

      > "$out"

      # extrage secventele 
      while IFS=' ' read -r header start end; do

        # verifica daca headerul se potriveste cu crz actual
        if [[ "$header" == "$chr" ]]; then
      
          seqtk subseq "$fasta" <(echo -e "$header $start $end") >> "$out"
        fi
      done < "$coords"
  
    done

  done < "$MAIN_DIR/chr_list.txt"
fi

cd $ref_coord
cp *extracted_sequences.fasta $transp_comp 2>/dev/null
cp *extr_seq.fasta $transp_comp 2>/dev/null

cd $cont_coord
cp *contigs_sequences.fasta $transp_comp

cd $transp_comp

#######################################################

echo " "
echo "finding transposons..."
echo " "

if [[ "$answer" == "no" ]]; then
  while read -r chr; do

    # ref
    makeblastdb -in "${chr}_extracted_sequences.fasta" -dbtype nucl -out "${chr}_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_extr_seq_db" -out "${chr}_blast_ref.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_ref.txt" > "${chr}_blast_ref_filtr.txt"

    # contigs
    makeblastdb -in "${chr}_contigs_sequences.fasta" -dbtype nucl -out "${chr}_cont_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_cont_extr_seq_db" -out "${chr}_blast_cont.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_cont.txt" > "${chr}_blast_cont_filtr.txt"

    # random
    for set in 1 2 3; do
   
      makeblastdb -in "${chr}_extr_seq_random_set${set}.fasta" -dbtype nucl -out "${chr}_extr_seq_random_set${set}_db" > /dev/null 2>&1
   
      blastn -query $transp_fasta -db "${chr}_extr_seq_random_set${set}_db" -out "${chr}_blast_random_set${set}.txt" -outfmt 6
     
      awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_random_set${set}.txt" > "${chr}_blast_random_set${set}_filtr.txt"

    done

  done < "$MAIN_DIR/chr_list.txt"
fi

if [[ "$answer" == "yes" ]]; then
  while read -r chr; do
    # ref - in centromere
    makeblastdb -in "${chr}_in_centr_extr_seq.fasta" -dbtype nucl -out "${chr}_in_centr_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_in_centr_extr_seq_db" -out "${chr}_blast_in_centr_ref.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_in_centr_ref.txt" > "${chr}_blast_in_centr_ref_filtr.txt"

    # ref - out of centromere
    makeblastdb -in "${chr}_out_of_centr_extr_seq.fasta" -dbtype nucl -out "${chr}_out_of_centr_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_out_of_centr_extr_seq_db" -out "${chr}_blast_out_of_centr_ref.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_out_of_centr_ref.txt" > "${chr}_blast_out_of_centr_ref_filtr.txt"

    # contigs - in centromere
    makeblastdb -in "${chr}_in_centr_contigs_sequences.fasta" -dbtype nucl -out "${chr}_in_centr_cont_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_in_centr_cont_extr_seq_db" -out "${chr}_blast_in_centr_cont.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_in_centr_cont.txt" > "${chr}_blast_in_centr_cont_filtr.txt"

    # contigs - out of centromere
    makeblastdb -in "${chr}_out_of_centr_contigs_sequences.fasta" -dbtype nucl -out "${chr}_out_of_centr_cont_extr_seq_db" > /dev/null 2>&1

    blastn -query $transp_fasta -db "${chr}_out_of_centr_cont_extr_seq_db" -out "${chr}_blast_out_of_centr_cont.txt" -outfmt 6

    awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_out_of_centr_cont.txt" > "${chr}_blast_out_of_centr_cont_filtr.txt"

    # random 
    for set in 1 2 3; do
   
      makeblastdb -in "${chr}_extr_seq_out_of_centr_random_set${set}.fasta" -dbtype nucl -out "${chr}_extr_seq_out_of_centr_random_set${set}_db" > /dev/null 2>&1
   
      blastn -query $transp_fasta -db "${chr}_extr_seq_out_of_centr_random_set${set}_db" -out "${chr}_blast_out_of_centr_random_set${set}.txt" -outfmt 6
     
      awk 'NR==FNR {len[$1]=$2; next} $1 in len {if ($4 >= 0.15 * len[$1]) print}' $transp_lengths "${chr}_blast_out_of_centr_random_set${set}.txt" > "${chr}_blast_out_of_centr_random_set${set}_filtr.txt"

    done

  done < "$MAIN_DIR/chr_list.txt"
fi

rm *_db.*



# ultimul script R care face tabelele csv

Rscript transposons_comparison.R


#######################################################

rm *chr.fasta *limits.txt chrom_sizes.txt random_coords.R *comparison.R $transp_fasta $transp_lengths *sequences.fasta
rm centromere_coordinates.txt 2>/dev/null 

# pune fisierle in foldere separate pt fiecare crz
while read chr; do
  mkdir -p "$transp_comp/$chr"  
  mv "$transp_comp"/"$chr"* "$transp_comp/$chr/" 2>/dev/null  
done < "$MAIN_DIR/chr_list.txt"

cd $ref_coord

rm *chr.fasta 
rm centromere_coordinates.txt 2>/dev/null

while read chr; do
  mkdir -p "$ref_coord/$chr"  
  mv "$ref_coord"/"$chr"* "$ref_coord/$chr/" 2>/dev/null 
done < "$MAIN_DIR/chr_list.txt"

cd $cont_coord

rm *chr.fasta *assembly.fasta contig_coord_new.R
rm *extracted_sequences.fasta *extr_seq.fasta 2>/dev/null

while read chr; do
  mkdir -p "$cont_coord/$chr"  
  mv "$cont_coord"/"$chr"* "$cont_coord/$chr/" 2>/dev/null 
done < "$MAIN_DIR/chr_list.txt"

cd $MAIN_DIR
rm chr_list.txt

echo " "
echo "done"
echo " "


