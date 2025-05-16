#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 16 
#SBATCH --mem 32G			
#SBATCH --time=0-01:00:00				
#SBATCH --output output_%x		
#SBATCH --mail-type=END,FAIL			
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	
#SBATCH -p ei-medium

source ~/.bashrc
mamba activate /hpc-home/mahony/miniforge3/envs/miniconda_dna

#echo "Generating promoter and transcript fasta files"
#sbatch extract_transcripts_and_promoters.sh

#echo "Generating k-mers..."
#python scripts/generate_kmers.py 6

#echo "Testing different common models"
#python scripts/train_model.py --models logistic_regression --cv 5

echo "Selecting top features..."
python scripts/feature_select.py 

#echo "Training models..."
#python scripts/train_model.py --models all --cv 5

#echo "Evaluating model performance..."
#python scripts/evaluate_model.py

#echo "Running SHAP explanation..."
#python scripts/shap_analysis.py --model random_forest

#echo "Pipeline complete."


