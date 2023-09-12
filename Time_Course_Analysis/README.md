Time course loops analysis
=============================

The code is for  performing a time course analysis given a set of interactions

In the parameter section, specify:

k27ac_up: File with coordinates of upregulated H3K27ac loops in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format. These files contained the observed/expected (OE) counts in day 0, day 3 and day 18 after IKAROS induction. OE counts can be extracted from the .hic file using the [straw](https://github.com/aidenlab/straw) library.

k27ac_down: File with coordinates of downregulated H3K27ac loops in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format. These files contained the observed/expected (OE) counts in day 0, day 3 and day 18 after IKAROS induction.

ctcf_up: File with coordinates of upregulated CTCF loops in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format. These files contained the observed/expected (OE) counts in day 0, day 3 and day 18 after IKAROS induction.

ctcf_down: File with coordinates of downregulated CTCF loops in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format. These files contained the observed/expected (OE) counts in day 0, day 3 and day 18 after IKAROS induction.


