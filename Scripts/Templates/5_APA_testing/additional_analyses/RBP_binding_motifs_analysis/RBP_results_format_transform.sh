#!/bin/zsh

dir='/Users/ninhle/Desktop/Research/scPASU_pipeline_runs/Ureter10_scPASU_run/outputs/differentiation_stage_cellranger_peakcount/rbp_binding_motifs/U10APAshortening3UTRNULL_Strict/'

fprefix=${dir}U10APAlengthening3UTR_Strict_All_Predictions
fprefix=${dir}U10APAlengthening3UTRNULL_Strict_All_Predictions
fprefix=${dir}U10APAshortening3UTR_Strict_All_Predictions
fprefix=${dir}U10APAshortening3UTRNULL_Strict_All_Predictions

input_file=${fprefix}.csv  # must be csv
output_file=${fprefix}.txt  # only txt

if [[ ! -f "$input_file" ]]; then
  echo "Error: Input file '$input_file' not found."
  exit 1
fi

# Process the file using awk
awk '
BEGIN {
  FS = ",";
  OFS = ",";
  current_seq = "";
  line_count = 0;
  header_printed = 0;
  old_header = "";
}
# Match any line that does not contain a comma which should be seq name
/^[^,]+$/ {
  current_seq = $0;
  # table_line_count = 0; # Reset line count for a new table
  next;
}
# Skip blank lines
/^$/ {
  next;
}
{
  if (NF == 7 && header_printed == 0) {
    print "SeqName," $0;
    old_header = $0;
    header_printed = 1;
    next;
  }
  
  if (NF == 7 && $0 != old_header) {
    print current_seq, $0;
  }
}' "$input_file" > "$output_file"

echo "Transformation complete. The new file is '$output_file'."
