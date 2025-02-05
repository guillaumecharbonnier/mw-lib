#!/bin/bash
#Thomas script:
#https://mail.google.com/mail/u/0/#search/modkit/FMfcgzQVxbhHLCGbPqJkQpMFKVGbhdWR
# Only for reference, not used in my workfow

# Start timer
start_time=$(date +%s.%N)

# Input arguments
name=$1
time=$2

case_dir="/data/longreads_basic_pipe/${name}/${time}"
reference="/data/ref/hs1/chm13v2.0.fa"
cpg_bed="/data/ref/hs1/hs1_CpG.bed"

bam="${case_dir}/${name}_${time}_hs1.bam"
output_dir="${case_dir}/5mC_5hmC"

if [ ! -d "$output_dir" ]; then
  mkdir "$output_dir"
fi

out_bed="${output_dir}/${name}_${time}_5mC_5hmC.bed"
if [ ! -f "$out_bed.gz" ]; then
  modkit pileup -t 150 "$bam" "$out_bed" --log-filepath "${output_dir}/pileup.log"
  crabz -f bgzf "$out_bed" -o "$out_bed.gz"
  tabix "$out_bed.gz"
  rm "$out_bed"
fi

out_dmr="/data/longreads_basic_pipe/${name}/diag/5mC_5hmC/${name}_diag_DMR_5mC_5hmC.bed"
if [ "$time" == "diag" ]; then
  mrd_bed="/data/longreads_basic_pipe/${name}/mrd/5mC_5hmC/${name}_mrd_5mC_5hmC.bed.gz"
  modkit dmr pair -m C \
    -a "$mrd_bed" -b "${out_bed}.gz" \
    --ref "$reference" \
    -o "$out_dmr" \
    -t 160 \
    -r "$cpg_bed"
fi

# Stop timer and calculate execution time
end_time=$(date +%s.%N)
execution_time=$(echo "$end_time - $start_time" | bc)
echo "Modkit for ${name} ${time} execution time: $execution_time seconds" >&2
"""

