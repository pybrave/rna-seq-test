
echo "{{index}} {{fastq1}} {{fastq2}}" > {{output_dir}}/{{meta.file_name}}.bam


cat > {{output_dir}}/outputs.json <<EOF
{
  "bam": " {{output_dir}}/{{meta.file_name}}.bam"
}
EOF