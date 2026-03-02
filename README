# sequence_merger

This repository contains a lightweight script to merge ordered
unitigs or contigs based on sequence overlaps.

The tool was developed to reproduce the *final path selection and
simplification stages* of long-read assemblers such as
[hifiasm](https://github.com/chhylp123/hifiasm) and conceptually similar
graph-based merging tools like Rukki. It is intended as a small, transparent post-assembly helper for manual
curation or the 'duo-phasing' approach where we use a child to phase parental assemblies.

------------------------------------------------------------------------

## Concept

Assemblers produce unitigs embedded in a graph. In complex regions —
especially during duo/trio phasing — you may want to:

- Inspect the graph manually  
- Decide on the biologically correct traversal  
- Explicitly merge unitigs along that path  

`sequence_merger` performs only the final step:

- Reads a multi-FASTA file  
- Takes an ordered list of contig names  
- Merges them sequentially using a defined overlap length  
- Outputs a curated merged sequence  

No graph inference is performed internally — the biological logic is
provided by the user.

------------------------------------------------------------------------

## Suggested Duo-Phasing Workflow

### 1. Generate a parental assembly

Assemble the parent using **hifiasm**.

### 2. Identify inherited unitigs

Run `yak trioeval` using the child as a pseudo-parent to obtain
parental unitigs likely inherited by the child.

### 3. Inspect the graph

Visualize the assembly graph in **Bandage** and determine the correct
unitig order.

### 4. Create an ordered list

Write the selected unitig names into a simple `.txt` file (see example_contiglist.txt)


### 5. Merge

Run the script to obtain a single curated merged assembly.

------------------------------------------------------------------------

## Usage

```bash
Rscript sequence_merger.R \
  -f sequences.fasta \
  -c contiglist.txt \
  -o 500 \   # minimal overlap length to allow contig merging
  -n merged_seq \
  -out output.fasta
```

------------------------------------------------------------------------

## Scope

This repository is mostly intended for documentation and reproducibility of curated assembly workflows in Höps et al. 2026, rather than as an
automated tool for general use. If you still want to use it please go ahead, and let me know if you need help. Rukki might be a more viable alternative for most usecases, though.

-----------------------------------------------------------------------

## Concept

MIT License. See LICENSE for details.

------------------------------------------------------------------------

## Contact

wolfram.hops@radboudumc.nl

