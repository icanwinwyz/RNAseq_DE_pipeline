#########################################################################
# File Name: RNAseq_analysis.sh
# Author: Yizhou Wang
#########################################################################
#!/bin/bash
#!/bin/bash
display_usage() {
	echo -e "NAME:\n  RNAseq_analysis."
	echo -e "\nDESCRIPTION:\n   This pipeline will use COUNT file to do doenstream differentially expressed genes (DEGs) analysis."
	echo -e "\nUsage:\n   bash RNAseq_analysis.sh Count_file.csv Sample_info.csv comparisons.csv "
    
    echo "Input options:"
    echo "   -h|--help    show this help"
    
    echo "Input files:"
    echo "   The first input is Count_file.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_COUNTS.csv  as an example format"
    echo "   The second input is Sample_info.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_sample_info.csv as an example format"
    echo "   The third input is Conparisons.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_comparisons.csv as an example format"
    echo "   The fourth input is project_ID, for example: AP-5782–11–08–2018"
	}
# check whether user had supplied -h or --help . If yes display usage
if [[ ($1 == "--help") ||  ($1 == "-h") ]]
then
	display_usage
	exit 0
fi

# get PCA plots for all samples and DEGs table, interactive report for each comparison
#/hpc/apps/R/3.4.1/bin/Rscript /common/genomics-core/data/Temp/Di_RNA_seq_test/downstream_test/RNAseq_tier2.R $1 $2 $3 $4
/usr/bin/Rscript /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_tier2.R $1 $2 $3 $4
