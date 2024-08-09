#!/bin/bash
#SBATCH -J rmsd
#SBATCH --partition="cpu(all)"
#SBATCH --time=99999:00:00
#SBATCH --cpus-per-task=4

##################################################################################################
#Edit file paths bellow
#################################################################################################
amber_container_path="/scratch/av144592/container/amber22_tools23_cuda.sif"

project_path="/projects/av144592"

ligand="nmp422"

protein="pparg"

parmtop_path="/projects/av144592/new_equilibrated/nmp422/build/nmp422_pparg_hmr.prmtop"

ref_trajectory_path="/projects/av144592/new_equilibrated/nmp422/build/9md.rst7"

pdb_file="/projects/av144592/new_equilibrated/nmp422/build/nmp422_pparg.pdb"

##################################################################################################

last_protein_res_num="not found"
ligand_res_num="not found"
ter_num=1

array=("HIS" "HID" "HIE" "CYS" "CYX" "ASP" "GLU" "LYS" "ASH" "GLH" "LYN" "ACE" "NHE" "NME" "ARG" "SER" "THR" "ASN" "GLN" "GLY" "PRO" "ALA" "VAL" "ILE" "LEU" "MET" "PHE" "TYR" "TRP")

while [ "$last_protein_res_num" == "not found" ]; do
        if [[ "${array[*]}" =~ "$fourth_column" ]]; then
                line_number=$(grep -n '^TER' "$pdb_file" | sed -n "${ter_num}p" | cut -d':' -f1)
                ((line_number++))
                fourth_column=$(awk 'NR=='"$line_number"' {print $4}' "$pdb_file")
                last_protein_res_num="not found"
                ((ter_num++))
        else
                ((line_number-=2))
                last_protein_res_num=$(awk 'NR=='"$line_number"' {print $5}' "$pdb_file")
        fi
done

ter_num=1
fourth_column="reset"

while [ "$ligand_res_num" == "not found" ]; do
        if [[ "$fourth_column" != "MOL" ]]; then
                line_number=$(grep -n '^TER' "$pdb_file" | sed -n "${ter_num}p" | cut -d':' -f1)
                ((line_number++))
                fourth_column=$(awk 'NR=='"$line_number"' {print $4}' "$pdb_file")
                ligand_res_num="not found"
                ((ter_num++))
        fi
        if [[ "$fourth_column" == "MOL" ]]; then

                ligand_res_num=$(awk 'NR=='"$line_number"' {print $5}' "$pdb_file")
        fi
done

echo "last protein res num $last_protein_res_num"
echo "ligand res num $ligand_res_num"

####################################################################################################################

all="1-$last_protein_res_num"

cat > "rmsd_rmsf.cpptraj" <<EOF
parm $parmtop_path
reference $ref_trajectory_path parm $parmtop_path [$ref_trajectory_path]
for rep in a,b,c,d,e
	trajin ../../$ligand-$protein-\$rep.nc
	autoimage
	rms $protein-\$rep-bb :$all@CA,C,O,N ref [$ref_trajectory_path] out $protein-bb-rmsd-$ligand.dat mass
	rms $ligand-\$rep :$ligand_res_num&!@H= ref [$ref_trajectory_path] nofit out $ligand-rmsd-$protein.dat mass 
	atomicfluct $protein-\$rep-rmsf out $protein-rmsf-$ligand.dat :$all byres
	atomicfluct $ligand-\$rep-rmsf  out $ligand-rmsf-$protein.dat :$ligand_res_num byatom
	go
	clear trajin
done
quit

EOF

apptainer exec --bind "/projects/av144592/" "$amber_container_path" cpptraj -i rmsd_rmsf.cpptraj

exit 0
