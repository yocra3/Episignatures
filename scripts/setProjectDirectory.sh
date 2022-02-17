#'#################################################################################
#'#################################################################################
# Set project folders' structure and link files
#'#################################################################################
#'#################################################################################

## Make folders
mkdir data
mkdir data/tcga_hdf5
mkdir results
mkdir reports
mkdir workflows
mkdir workflows/bin
mkdir -p docker_sources/r_session
mkdir -p docker_sources/python_session

rsync -azvh --progress workflows/ cruizg@marvin.s.upf.edu:/gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures/workflows/
rsync -azvh --progress conf/ cruizg@marvin.s.upf.edu:/gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures/conf/      
rsync -azvh --progress modules/ cruizg@marvin.s.upf.edu:/gpfs42/robbyfs/scratch/lab_laperez/cruizg/Episignatures/modules/
