suppressMessages(library(Biostrings))
suppressMessages(library(argparse))


parser <- ArgumentParser()
parser$add_argument("-f", "--fasta", required=TRUE, help="Multi-sequence FASTA file")
parser$add_argument("-c", "--contiglist", required=TRUE, help="Sorted contiglist for merging")
parser$add_argument("-o", "--overlap", default=500, type="integer", help="Overlap length")
parser$add_argument("-n", "--contigname", default='seq', type="character", help="Basename of the merged contigs")
parser$add_argument("-out", "--output", required=TRUE, help="Output file for merged sequences")

args <- parser$parse_args()

args = list()
args$fasta = "/Users/Z364220/projects/15q13/trio2/child_t2_line.fa"
args$contiglist = "/Users/Z364220/projects/15q13/trio2/t2_child_utigs.txt"
args$overlap = 500
args$contigname = 'seq'
args$output = ''#/Users/Z364220/projects/mini-programs/contig-merger/data/trio4/mother/mf/mf.fa'

merge_sequences <- function(seqs, overlap_length, max.mismatch = 3) {
  idx_baseseq <- 1
  
  while (1) {  # Run until there are no sequences left

    found <- FALSE
    
    current_base_seq <- seqs[idx_baseseq]
    last_part <- subseq(current_base_seq, start = width(current_base_seq) - overlap_length + 1)

    for (i in seq_along(seqs)) {
      if (i == idx_baseseq) next # Do not compare to self
      candidate_seq <- seqs[i]
      candidate_parts <- c(candidate_seq, reverseComplement(candidate_seq))

        for (n_candidate in seq_along(candidate_parts)) {

          candidate = candidate_parts[n_candidate]
          match_position <- matchPattern(last_part[[1]], candidate[[1]])
          if ((length(match_position) > 1)){
            stop("Multiple matches found within one seq. Consider increasing window size, or check your fasta.")
          }
          
          if ((length(match_position) == 0)){
            next # Next candidate
          } 
            
          if (end(match_position) == width(candidate)){
            next #Next candidate
          }
          

          end_base_seq = subseq(seqs[idx_baseseq], start = max(1,width(current_base_seq) - (start(match_position) + overlap_length) + 2))
          start_candidate = subseq(candidate, start = 1, end = end(match_position))

          # Define the maximum length allowed for pattern matching
          max_length <- 20000
          
          # Trim the sequences to keep only the last max_length bases
          end_base_seq_trimmed <- if (width(end_base_seq) > max_length) {
            subseq(end_base_seq, start = width(end_base_seq) - max_length + 1)
          } else {
            end_base_seq
          }
          
          start_candidate_trimmed <- if (width(start_candidate) > max_length) {
            subseq(start_candidate, start = width(start_candidate) - max_length + 1)
          } else {
            start_candidate
          }
          
          match_certified = matchPattern(end_base_seq_trimmed[[1]], start_candidate_trimmed[[1]], max.mismatch = max.mismatch, with.indels = T)
          
          if ((length(match_certified) == 0)){
            print(paste0('Rejecting a merge between ', names(seqs)[idx_baseseq], ' and ', names(seqs)[i]))
            #browser()
            next # Next candidate
          }              
                       
          start_pos <- start(match_position)[1] # Take the first match's start position
          new_seq <- paste0(as.character(current_base_seq), as.character(subseq(candidate, start = start_pos + overlap_length)))
          current_name <- names(seqs)[idx_baseseq]
          candidate_name <- names(seqs)[i]
          
          if(n_candidate == 2){
            candidate_name = paste0(candidate_name, '_r')
          }
          new_name <- paste(current_name, candidate_name, sep="->")

          seqs[idx_baseseq] <- DNAStringSet(DNAString(new_seq))
          
          if (width(seqs[idx_baseseq]) < width(seqs[i])){
            # Print which sequences are affected
            print(names(seqs[i]))
            print(names(seqs[idx_baseseq]))

            print("Seq has been shortened by merge! Likely error")
            browser()
          }
          names(seqs)[idx_baseseq] <- new_name
          
          seqs <- seqs[-i]
          
          # Patch to fix fases where the seqs[-i] kicks out a previous sequence, 
          # changing all index values of the later sequences.
          if (i < idx_baseseq){
            idx_baseseq = idx_baseseq - 1
          }
          found <- TRUE
          
          print(paste0("Merged ",
                      current_name,
                       ' and ',
                      candidate_name,
                       ' to form ',
                       new_name
                ))

          # print(paste0("Merged ",
          #              as.character(current_base_seq),
          #              ' and ',
          #              substring(as.character(candidate), start_pos + overlap_length),
          #              ' to form ',
          #              as.character(new_seq)
          #       ))
          # 
          break # Stop checking the other candidate
        }
      
      if (found) break # Break and start again with that idx_base
    }
    
    # If we were successful, we want to give the idx another go. 
    # If not, we move to the next.
    if (!found){
      idx_baseseq <- idx_baseseq + 1
    }
    if (length(seqs) == 1 || idx_baseseq > length(seqs)) break
  }

  return(seqs)
}

seqs <- readDNAStringSet(filepath = args$fasta)
names_contigs_to_merge <- read.table(args$contiglist, header = FALSE, stringsAsFactors = FALSE)$V1



seqs_contigs_to_merge = seqs[names_contigs_to_merge,]
merged_seqs <- 
    merge_sequences(seqs_contigs_to_merge, args$overlap, max.mismatch=10)



dim(merged_seqs)
#merged_seqs <- merge_sequences(merged_seqs, args$overlap, max.mismatch = 13)

#seqs_contigs_to_merge[1] = reverseComplement(seqs_contigs_to_merge[1])
#merged_seqs <- merge_sequences(seqs_contigs_to_merge, args$overlap)

names(merged_seqs)



names(merged_seqs) = paste0(args$contigname, seq(1, length(merged_seqs)))

# Save as fasta file
Biostrings::writeXStringSet(merged_seqs, args$output)

# Inform the user what was done and where the saved file is
print(paste0("Merged ", length(seqs_contigs_to_merge), " sequences into ", length(merged_seqs), " sequences."))
print(paste0("Saved merged sequences to ", args$output))
