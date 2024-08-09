#!/bin/bash

amber_container_path="/scratch/av144592/container/amber22_tools23_cuda.sif"
bind_path="/projects/av144592"

ligand_name="nmp422"
protein_name="pparg"

parmtop_path="/projects/av144592/new_equilibrated/nmp422/build/nmp422_pparg_hmr.prmtop"
trojectory_path="/projects/av144592/new_equilibrated/nmp422/build/9md.rst7"

cat > "md.in" <<EOF
MD explicit solvent heavy atom no rest shake dt 0.004 for 100ns
&cntrl
	imin=0, irest=0, ntx=1,
	ntt=3, temp0=300.0, gamma_ln=3, ig=-1,
	ntp=1, taup=2.0,ntb=2,
	ntc=2, ntf=2, cut=8,
	dt=0.004, nstlim=250000000,
	ntpr=250000, ntwx=250000, ntwr=250000,
/
EOF

# Function to create the Amber MD simulation script
create_amber_md_script() {
    local rep="$1"
cat > $rep-md.sh <<EOF
#!/bin/bash
#SBATCH --job-name="${rep}-${ligand_name}"
#SBATCH --partition="gpu(all)"
#SBATCH --time=99999:00:00
#SBATCH --gres=gpu:rtx_2080:1 ##reserve 1 gpu

set -x

name=$ligand_name-$protein_name-$rep

apptainer exec --nv --bind $bind_path $amber_container_path pmemd.cuda -O -i md.in -p $parmtop_path -c $trojectory_path -r $trojectory_path -o $ligand_name-$protein_name-$rep.mdout -r $ligand_name-$protein_name-$rep.rst -x $ligand_name-$protein_name-$rep.nc -inf mdinfo-$rep

EOF
    chmod +x $rep-md.sh
    sbatch $rep-md.sh
}

# Loop through the directories and create Amber MD scripts
for rep in a b c d e; do
    create_amber_md_script "$rep"
done
