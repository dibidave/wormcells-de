
# C. elegans single cell gene expression

The `wormcells-de` app allows you to perform differential expression on data from C. elegans single cell RNA sequencing (scRNAseq). 21 experiments from 3 different studies were integrated and can be compared. That means you can select cells from two different experiments and perform differential expression on them!


Just select cell types and experiments to compare, some genes to highlight in your volcano plot, and leave
    your email. You will receive an [interactive volcano plot](https://scvi-differential-expression.s3.us-east-2.amazonaws.com/plots/eduardo%40wormbase.org%4020200227-233946-results.html), a
    [csv file with results](https://scvi-differential-expression.s3.us-east-2.amazonaws.com/csv/eduardo%40wormbase.org%4020200227-233946-results.csv) and [a csv file with the selected groups and genes](https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/eduardo%40wormbase.org%2520200227-233946.csv").
        The example files linked show the result of the example submission below

 Results should arrive in less than 15 minutes. If they take more than an hour, something broke, so let me know by writing to eduardo@wormbase.org. Also feel free to write me if you have any feedback.


The single cell gene count matrixes were processed using a machine learning framework called [Single-cell Variational Inference (scVI)](https://github.com/YosefLab/scVI). The scVI framework enables integrating data from different sources (different experiments, batches and technologies), clustering and label transfer, and performing differential expression between clusters.

To learn how the plots were made, check out [this Jupyter notebook](https://colab.research.google.com/drive/1hF7KSujhhHcyxzWkzAHy9lazXLexainr) which performs differential expression as is here here. It runs on Google Colab so you can start playing right away.

The `wormcells-de` app is still in development and on this tool will inform how WormBase may incorporate and display single cell data in the future. 
