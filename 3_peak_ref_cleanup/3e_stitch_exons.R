library(Biostrings)
library(argparser)

parser<-arg_parser(name="3e_stitch_exons.R",description="Stitch exons into transcripts")

parser<-add_argument(
  parser,
  arg='--input_file',
  short = '-i',
  type="character",
  help="Enter input file path. Format: path/to/input/file/")

parser<-add_argument(
  parser,
  arg='--output_file',
  short = '-o',
  type="character",
  help="Enter output file path. Format: path/to/output/file/")

args <- parse_args(parser)

input_file <- args$input_file
output_file <- args$output_file

sequences <- readDNAStringSet(input_file)
names(sequences) <- sapply(strsplit(names(sequences), ':'), '[[', 1)

grouped_sequences <- lapply(split(sequences, names(sequences)), function(group) {
  unlist(group)
})

grouped_sequences <- DNAStringSet(unlist(grouped_sequences))

writeXStringSet(grouped_sequences, output_file, format = "fasta")
