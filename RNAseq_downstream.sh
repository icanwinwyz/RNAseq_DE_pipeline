#########################################################################
# File Name: RNAseq_analysis.sh
# Author: Di Wu
# mail: di.wu@cshs.org
# Created Time: Feb., 2019
#########################################################################
#!/bin/bash
#!/bin/bash
display_usage() {
        echo -e "NAME:\n  RNAseq_analysis."
        echo -e "\nDESCRIPTION:\n   This pipeline will use COUNT file to do doenstream differentially expressed genes (DEGs) analysis."
        echo -e "\nUsage:\n   bash RNAseq_analysis.sh Count_file.csv Sample_info.csv conparisons.csv "
    echo "Input options:"
    echo "   -h|--help    show this help"
    echo "Input files:"
    echo "   The first input is Count_file.csv, check demo_COUNTS.csv  as an example format"
    echo "   The second input is Sample_info.csv, check demo_sample_info.csv as an example format"
    echo "   The third input is Conparisons.csv, check demo_comparisons.csv as an example format"
    echo "   The fourth input is project_ID, for example: demo"
        }
# check whether user had supplied -h or --help . If yes display usage
if [[ ($1 == "--help") ||  ($1 == "-h") ]]
then
        display_usage
        exit 0
fi

# get PCA plots for all samples and DEGs table, interactive report for each comparison
#/hpc/apps/R/3.4.1/bin/Rscript /common/genomics-core/data/Temp/Di_RNA_seq_test/downstream_test/RNAseq_tier2.R $1 $2 $3 $4
/usr/bin/Rscript /home/genomics/genomics/data/Temp/Di_RNA_seq_test/downstream_test/RNAseq_tier2.R $1 $2 $3 $4
