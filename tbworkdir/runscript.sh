#! /bin/bash
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1:nsf

# This job's working directory
CPUS=`cat $PBS_NODEFILE | wc -l`
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
NODES=`cat $PBS_NODEFILE`
echo $NODES
echo This job has allocated $CPUS nodes

/sw/local/MATLAB/R2011b/bin/matlab -nodisplay -nodesktop -r "ptags=$tagid;run_tag;exit"
