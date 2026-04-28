
echo "{{reference_genome}}" > output/mouse.1.bt2
echo "{{reference_genome}}" > output/mouse.2.bt2
echo "{{reference_genome}}" > output/mouse.3.bt2
echo "{{reference_genome}}" > output/mouse.4.bt2
echo "{{reference_genome}}" > output/mouse.rev.1.bt2
echo "{{reference_genome}}" > output/mouse.rev.2.bt2


cat > {{output_dir}}/outputs.json <<EOF
{
  "index": "{{output_dir}}/mouse.1.bt2"
}
EOF