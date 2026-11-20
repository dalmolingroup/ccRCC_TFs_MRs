#!/bin/bash
#SBATCH --job-name=trainnet
#SBATCH --partition=amd-512
#SBATCH --nodes=1
#SBATCH --cpus-per-task=33          # 32 workers do snow + 1 processo principal
#SBATCH --mem=200G
#SBATCH --time=05-00:00:00          # ajuste ao limite da partição (sinfo -p amd-512)
#SBATCH --output=trainnet_%j.log

# Comentário 2 (R2): rede + MRA só no treino, dentro do container rocker 4.5.1
# (o mesmo onde os pacotes de ~/R/rocker-rstudio foram instalados)

ulimit -n 8192
module load singularity

workdir=~/scratch/rstudio
mkdir -p -m 700 ${workdir}/tmp

export CONTAINER_FILE=/opt/npad/shared/containers/rocker-4.5.1-extra.sif
export SINGULARITY_BIND="${workdir}/tmp:/tmp,/scratch/global/$USER"

# --cleanenv apaga as variáveis do host; as que o R precisa vão com o prefixo SINGULARITYENV_
export SINGULARITYENV_R_LIBS_USER=${HOME}/R/rocker-rstudio
export SINGULARITYENV_OMP_NUM_THREADS=1      # o paralelismo vem dos 32 workers do snow

echo "Início: $(date) | nó: $(hostname) | CPUs: ${SLURM_CPUS_PER_TASK}"

singularity exec --cleanenv $CONTAINER_FILE \
  Rscript ${HOME}/03b_rtn_network_train_only.R   # ajuste o caminho do script no HPC

echo "Fim: $(date)"
