#!/bin/bash

#### Edit the variables below ####

container_path="/projects/av144592/container/amber22_tools23_cuda.sif"
bind_path="/projects/av144592/"

protein_name="pparg"
protein_pdb="2Q5P_pparg.pdb"

ligand_name="mrl24"
ligand_sdf="2Q5P_mrl24.sdf"
ligand_net_molecular_charge="-1"

#### Stop editing ####

#converts the protein pdb to an amber pdb
apptainer exec --bind $bind_path $container_path pdb4amber -i "$protein_pdb" -o "$protein_name-4amber.pdb"

#converts the ligand pdv to a mol2 file
apptainer exec --bind $bind_path $container_path antechamber -i "$ligand_sdf" -fi sdf -o "$ligand_name.mol2" -fo mol2 -c bcc -nc "$ligand_net_molecular_charge" -at gaff2

#get parameters for ligand
apptainer exec --bind $bind_path $container_path parmchk2 -i "$ligand_name.mol2" -f mol2 -o "$ligand_name.frcmod"


#wirghts tleap script to determin the volume of the solvated oct box
rm leap.log
touch volume.in

cat <<EOF > volume.in
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc

MOL = loadmol2 $ligand_name.mol2
loadamberparams $ligand_name.frcmod
check MOL
saveoff MOL $ligand_name.lib
saveamberparm MOL $ligand_name.prmtop $ligand_name.rst7

loadoff $ligand_name.lib
protein = loadpdb $protein_name-4amber.pdb
complex = combine {protein MOL}

check complex
charge complex

addions complex K+ 0
addions complex Cl- 0

solvateoct complex OPCBOX 10.0

saveamberparm complex $protein_name-$ligand_name.prmtop $protein_name-$ligand_name.rst7

quit
EOF

#runs the tleap sctipt to find the volume of the solvated oct box
apptainer exec --bind $bind_path $container_path tleap -f volume.in

#extractes the volume of the oct box from the leap.log and determins how many ion need to be added to get 0.05M K+CL-
volume_a=$(grep "Volume:" leap.log | awk '{print $2}')
echo "The volume is: $volume_a"
volume_L=`echo "scale=length($volume_a)+26;$volume_a/(10^27)" | bc -l`
echo "The volume in L: $volume_L"
ion=`echo "scale=6;$volume_L*3.011*(10^22)" | bc -l`
echo "Number of ions for 150nM: $ion"
rounded_ion=$(printf "%.0f" "$ion")
echo "Number of ions for 150nM (rounded): $rounded_ion"


#wrights a tleap script build the final prmtop, rst7 and pdb
touch add_150mM_K+Cl-.in
cat <<EOF > add_150mM_K+Cl-.in
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc

MOL = loadmol2 $ligand_name.mol2
loadamberparams $ligand_name.frcmod
loadoff $ligand_name.lib

protein = loadpdb $protein_name-4amber.pdb
complex = combine {protein MOL}

check complex
charge complex

addions complex K+ 0
addions complex Cl- 0

solvateoct complex OPCBOX 10.0

addionsrand complex K+ $rounded_ion
addionsrand complex Cl- $rounded_ion

saveamberparm complex $protein_name-$ligand_name.prmtop $protein_name-$ligand_name.rst7
savepdb complex final-$protein_name-$ligand_name.pdb

quit
EOF

#runs the tleap script for the final prmtop, rst7, and pdb
apptainer exec --bind $bind_path $container_path tleap -f add_150mM_K+Cl-.in

#wrights a tleap script to get a hydrogen mass repartition prmtop file
cat <<EOF > hmassrepartition.parmed
parm $protein_name-$ligand_name.prmtop
hmassrepartition
outparm $protein_name-$ligand_name-hmr.prmtop
go
quit
EOF

apptainer exec --bind $bind_path $container_path parmed -i hmassrepartition.parmed

