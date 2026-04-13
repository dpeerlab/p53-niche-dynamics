#!/usr/bin/bash
m=$(pwd)
f=$m/$1
o=$f/MotifData/homer_tile_based
#srun --time=12:00 --mem=64 --partition=cpu "findMotifsGenome.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 homer_motif_scan -mask -find $m/custom_motifs/homerKnownMotifs.motif > $o/TileData_homerKnownMotif_masked.position"
#srun --time=12:00 --mem=64 --partition=cpu "findMotifsGenome.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 homer_motif_scan -find $m/custom_motifs/homerKnownMotifs.motif > $o/TileData_homerKnownMotif_unmasked.position"
#srun --time=12:00 --mem=64 --partition=cpu "annotatePeaks.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 -noann -nogene -mask -mbed $o/TileData_homerKnownMotif_masked.bed -m $m/custom_motifs/homerKnownMotifs.motif > $o/TileData_homerKnownMotif_masked.ann"
#srun --time=12:00 --mem=64 --partition=cpu "annotatePeaks.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 -noann -nogene -mbed $o/TileData_homerKnownMotif_unmasked.bed -m $m/custom_motifs/homerKnownMotifs.motif > $o/TileData_homerKnownMotif_unmasked.ann"
srun --time=12:00 --mem=64 --partition=cpu "findMotifsGenome.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 homer_motif_scan -mask -find $m/custom_motifs/p53_custom.motif > $o/TileData_homerp53_masked.position"
srun --time=12:00 --mem=64 --partition=cpu "findMotifsGenome.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 homer_motif_scan -find $m/custom_motifs/p53_custom.motif > $o/TileData_homerp53_unmasked.position"
srun --time=12:00 --mem=64 --partition=cpu "annotatePeaks.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 -noann -nogene -mask -mbed $o/TileData_homerp53_masked.bed -m $m/custom_motifs/p53_custom.motif > $o/TileData_homerp53_masked.ann"
srun --time=12:00 --mem=64 --partition=cpu "annotatePeaks.pl $f/TileData/Tile_seacell_tile_100_mode_insertion_extended_20.bed mm10 -noann -nogene -mbed $o/TileData_homerp53_unmasked.bed -m $m/custom_motifs/p53_custom.motif > $o/TileData_homerp53_unmasked.ann"
