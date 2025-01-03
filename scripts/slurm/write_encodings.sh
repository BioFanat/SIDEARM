#!/bin/csh -f
#SBATCH --mem=12G
#SBATCH --output=script_out/write_encodings.out # File to which stdout will be written
#SBATCH --error=script_out/write_encodings.err

# Note the different syntax for command substitution and variable referencing in csh
python ../write_encodings.py