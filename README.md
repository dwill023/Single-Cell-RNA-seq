# Single-Cell-RNA-seq

One of the most popular single cell isolation techniques is with the Drop-seq approach, where a single cell is encapsulated into a micordroplet containing a bead with unique barcodes, primers and enzymes where cDNA synthesis and library generation is performed.

![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/drop-seq.png) <sub>[Macosko et al., 2015](https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415005498%3Fshowall%3Dtrue)</sub>

One of the most widely used technologies for the Drop-seq approach is the 10x Genomics Chromium platform. There has been a study by [Svensson et al,2017](https://www.nature.com/articles/nmeth.4220) comparing the different scRNA methods. If transcript level quantificaion and detection is what your experiment is aiming for, then the drop-seq (10x chromium) is fine. If transcript isoforms are part of the study then the full-length method (Smart-seq2) is required. 

Since my lab focuses on comparing transcript levels I will only outline the experimental steps for the 10x chromium protocol. Below is a schematic of the chromium chip where the cells that have been enzymatically broken up and mixed with reagent into the 2nd row of the chip are added. The microfluidics will mix these cells with the beads and oil to form the gel emulsion (GEMs) at the top row. The GEMs from this row are taken off the chip and put on a thermocycler for cDNA synthesis which will also incorporate a unique barcode. Then the GEMs are broken open, amplified, and cleaned up for sequencing. The remainder of the schematic outlines the data processing and visualization steps. 

![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/workflow.png)

When setting up the experiment it is **important** that samples should be balanced, evenly distributed, across all stages of the experiment. This will reduce sources of technical variation in the experiment. 

For example, you have samples on Day 0 sequenced on one flowcell and Day 7 samples run on another flowcell. The variation you observe can not be determined to be from the stage, or from technical variation in the sequencing run. Therefore, samples from day 0 and day 7 should be run on both flowcells. 


