
echo {{ bam | join('   ') }} > {{output_dir}}/vcf.bam
cat > {{output_dir}}/outputs.json <<EOF
{
  "vcf": " {{output_dir}}/{{output_dir}}/vcf.bam"
}
EOF

