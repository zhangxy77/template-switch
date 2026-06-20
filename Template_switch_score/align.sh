input_dir="./muscle_data"
output_dir="./muscle_align"


mkdir -p "$output_dir"


for fasta_file in "$input_dir"/*.fasta; do
    
    base_name=$(basename "$fasta_file" .fasta)
    
    
    output_file="$output_dir/${base_name}_aligned.fasta"
    
    
    muscle -in "$fasta_file" -out "$output_file"
    
    
    echo "Processed: $fasta_file -> $output_file"
done

echo "All files processed!"