#!/bin/bash
#SBATCH -J NN_RCAL
#SBATCH -e NN_RCAL.err
#SBATCH -o NN_RCAL.log
#SBATCH --partition=batch
#SBATCH --ntasks=22
#SBATCH --cpus-per-task=1
#SBATCH --time=99-99:99:99
#SBATCH --nodes=1
export OMP_NUM_THREADS=1

echo
echo `date`
echo "     Directory: `pwd`"
echo "          Host: `hostname`"
echo "----------------" >> NN_RCAL.qlog
echo "started running at `date`" >> NN_RCAL.qlog
python3 generate_nn_database.py   ../../../Build/intel_linux_64/radcal_nn   --input RADCAL_NN_TEST.IN   --output-prefix radcal_nn_test --processes 22
echo "finished running at `date`" >> NN_RCAL.qlog
