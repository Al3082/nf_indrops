The objective of this project is to to process data from two SRA datset using a nextflow process. Data are downloaded priorly. They will be processed in three parts v3 stages 8 to 14, v3 stages 16 to 22 and then v2. This is due to disk space constraints on the hpc

These data were generated using indrops V2 and V3 and the Kotov represents a resequencing of libraries from briggs
where as briggs data are already demultiplex, this is not the case of kotov data

The process goes as follows :
Demultiplex the Kotov data using cutadapt and allowing one error. description of sequences to use are in inputs

The libraries demultiplexed from Kotov must be matched with the corresponding briggs ones

Then they should be processed using starsolo with options: --soloType CB_UMI_Complex, --soloCBmatchWLtype ED2

This pipeline should be tested on a small subset of data.

Several prelimary tasks should be attempted. Check the good understanding of the correspondance of the libraries from kotov and briggs
Checking how the tools would be run on an hpc. The pipeline should be hardware independant as much as possible

For the precise HPC that the pipeline will be used on first, it uses module load (ex module load singularity), any info can be asked to user on this matter

paths.txt describes organisation on this cluster

inputs and nf_pipelines come from a prior work that processed the same data and should be taken as such


