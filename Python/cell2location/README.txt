The seurat visium HD interface is a huge pain in the ass to deal with so here is what I've done

2026_06_18 renaming barcodes
- Take the barcodes output from loupe and change them to be unique to each image, in the visium object they just have random suffixes so I changed it to be from each image
- Rename the barcodes in the visium object to be of the same format
- Subset the visium object to ONLY include barcodes in the dissection & ONLY include the sketch assay so I can actually run it on my computer
-Output each part of the seurat object so we can reconstruct the adata object in model_fit.py

THIS is then input into the model_fit py which was fit on the multiome object with the final_clusters level


I did this in both windows and mac