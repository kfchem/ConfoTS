#!/bin/bash
#PBS -N confots_job
#PBS -q accel
#PBS -l select=1:ncpus=0:mem=0gb
#PBS -l walltime=168:00:00

cd $PBS_O_WORKDIR

# Activate Python environment (modify as needed)
# source /app/venv/accel/bin/activate

# Optional: copy default templates (only if needed)
# mkdir -p template
# cp -n /app/scripts/confots/template/* template/

# Run ConfoTS workflow steps
python confots_step0.py
python confots_step1.py
python confots_step2.py
python confots_step3.py
python confots_step4.py
python confots_step5.py