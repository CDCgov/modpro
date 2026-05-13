#!/bin/bash -l
##--- Grid Engine Directives ---##
#$ -pe openmpi-fillup 2
#$ -cwd
#$ -o results.txt
#$ -e results.txt
#$ -q gpu.q

#echo "Installing PyTorch and related packages..."
#pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124

#echo "Installing OpenFold3 with CUDA support..."
#pip install "openfold3[cuda]"

#echo "Loading openfold3 environment"
#conda activate openfold3_dockerless

#echo "Loading CUDA module..."
#ml cuda/12.9.1

#echo "Checking NVIDIA GPU status..."
#nvidia-smi

#run_openfold predict --query_json=test_HA.json
time snakemake -s ~/openfold3/openfold-3/scripts/snakemake_msa/MSA_Snakefile --configfile test_msa_gen.json --cores 64 --jobs 3 --nolock --keep-going --latency-wait 120
#echo "Script completed."
