#!/bin/bash
#SBATCH --job-name="relax"
#SBATCH --partition="gpu(all)"
#SBATCH --time=99999:00:00
#SBATCH --gres=gpu:rtx_2080:1 ##reserve 1 gpu

# Minimization for explicit solvent with a protein ligand complex

#### Edit the variables below ####

path_to_amber_container="/projects/av144592/container/amber22_tools23_cuda.sif"

bind_path="/projects/av144592"

prmtop_file="pparg-mrl24.prmtop"

inpcrd_file="pparg-mrl24.rst7"

pdb_file="final-pparg-mrl24.pdb"

#### Stop editing ####

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

# 5,000 steps of steepest descent then 10,000 steps of conjugate gradient with 500 kcal·mol-1·Å-2 restraints on the heavy atoms of the protein
cat > 1min.in <<EOF
minimization of solvent
&cntrl
	imin = 1, maxcyc = 15000,
	ncyc = 5000, ntx = 1,
	ntwe = 0, ntwr = 500, ntpr = 50
	ntc = 2, ntf = 2, ntb = 1, ntp = 0,
	cut = 10.0,
	ntr = 1, restraintmask = ':1-$last_protein_res_num & !@H=',
	restraint_wt = 500.,
	ioutfm = 1, ntxo = 2,
/
EOF

# 1,000 steps of steepest descent then 500 steps of conjugate gradient with 100 kcal·mol-1·Å-2 restraints on the heavy atoms of the protein
cat > 2min.in <<EOF
minimization of solvent
&cntrl
        imin = 1, maxcyc = 1500,
        ncyc = 1000, ntx = 1,
        ntwe = 0, ntwr = 500, ntpr = 50
        ntc = 2, ntf = 2, ntb = 1, ntp = 0,
        cut = 10.0,
        ntr = 1, restraintmask = ':1-$last_protein_res_num & !@H=',
	restraint_wt = 100.,
        ioutfm = 1, ntxo = 2,
/
EOF

# 1,000 steps of steepest descent then 500 steps of conjugate gradient with 1 kcal·mol-1·Å-2 restraints on the heavy atoms of the protein
cat > 3min.in <<EOF
minimization of solvent
&cntrl
        imin = 1, maxcyc = 1500,
        ncyc = 1000, ntx = 1,
        ntwe = 0, ntwr = 500, ntpr = 50
        ntc = 2, ntf = 2, ntb = 1, ntp = 0,
        cut = 10.0,
        ntr = 1, restraintmask = ':1-$last_protein_res_num & !@H=',
        restraint_wt = 1.,
        ioutfm = 1, ntxo = 2,
/
EOF

# 1,000 steps of steepest descent then 500 steps of conjugate gradient with no restraints on the heavy atoms of the protein
cat > 4min.in <<EOF
minimization of solvent
&cntrl
        imin = 1, maxcyc = 1500,
        ncyc = 1000, ntx = 1,
        ntwe = 0, ntwr = 500, ntpr = 50
        ntc = 2, ntf = 2, ntb = 1, ntp = 0,
        cut = 10.0,
        ioutfm = 1, ntxo = 2,
/
EOF


# heat up the system to 300K with 50 kcal·mol-1·Å-2 restraints on heavy atoms
cat > 5md.in <<EOF
heat up system to 300K over 1ns using 50 kcal·mol-1·Å-2 restraints on heavy atoms
&cntrl
	imin = 0, nstlim = 1000000, dt=0.001,
	irest = 0, ntx = 1, ig = -1,
	tempi = 10.0, temp0 = 300.0,
	ntc = 2, ntf = 2, tol = 0.00001,
	ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
	cut = 8.0, iwrap = 0,
	ntt = 3, gamma_ln = 1., ntb = 1, ntp = 0,
	nscm = 0,
	ntr = 1, restraintmask = ':1-$ligand_res_num & !@H=', restraint_wt = 50.0,
	nmropt = 1,
	ioutfm = 1, ntxo = 2,
/
&wt TYPE="TEMP0", istep1=0, istep2=1000000, value1=10., value2=300., /
&wt TYPE="END", /
EOF

# relax the system at a constant pressure for 1ns using 50 kcal·mol-1·Å-2 restraints on heavy atoms
cat > 6md.in <<EOF
relax the system at a constant pressure for 1ns using 50 kcal·mol-1·Å-2 restraints on heavy atoms
&cntrl
        imin = 0, nstlim = 1000000, dt=0.001,
        irest = 1, ntx = 5, ig = -1,
        temp0 = 300.0,
        ntc = 2, ntf = 2, tol = 0.00001,
        ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
        cut = 8.0, iwrap = 0,
        ntt = 3, gamma_ln = 1., ntb = 2, ntp = 1, barostat = 2
        nscm = 0,
        ntr = 1, restraintmask = ':1-$ligand_res_num & !@H=', restraint_wt = 50.0,
        ioutfm = 1, ntxo = 2,
/
EOF

# relax the system at a constant pressure for 1ns using 10 kcal·mol-1·Å-2 restraints on heavy atoms
cat > 7md.in <<EOF
relax the system at a constant pressure for 1ns using 10 kcal·mol-1·Å-2 restraints on heavy atoms
&cntrl
	imin = 0, nstlim = 1000000, dt=0.001,
        irest = 1, ntx = 5, ig = -1,
        temp0 = 300.0,
        ntc = 2, ntf = 2, tol = 0.00001,
        ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
        cut = 8.0, iwrap = 0,
        ntt = 3, gamma_ln = 1., ntb = 2, ntp = 1, barostat = 2
        nscm = 0,
        ntr = 1, restraintmask = ':1-$ligand_res_num & !@H=', restraint_wt = 10.0,
        ioutfm = 1, ntxo = 2,
/
EOF

# relax the system at a constant pressure for 1ns using 2 kcal·mol-1·Å-2 restraints on heavy atoms
cat > 8md.in <<EOF
relax the system at a constant pressure for 1ns using 2 kcal·mol-1·Å-2 restraints on heavy atoms
&cntrl
	imin = 0, nstlim = 1000000, dt=0.001,
        irest = 1, ntx = 5, ig = -1,
        temp0 = 300.0,
        ntc = 2, ntf = 2, tol = 0.00001,
        ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
        cut = 8.0, iwrap = 0,
        ntt = 3, gamma_ln = 1., ntb = 2, ntp = 1, barostat = 2
        nscm = 0,
        ntr = 1, restraintmask = ':1-$ligand_res_num & !@H=', restraint_wt = 2.0,
        ioutfm = 1, ntxo = 2,
/
EOF

# relax the system at a constant pressure for 1ns using 0 kcal·mol-1·Å-2 restraints on heavy atoms
cat > 9md.in <<EOF
relax the system at a constant pressure for 1ns using 0 kcal·mol-1·Å-2 restraints on heavy atoms
&cntrl
	imin = 0, nstlim = 1000000, dt=0.001,
	irest = 1, ntx = 5, ig = -1,
        temp0 = 300.0,
        ntc = 2, ntf = 2, tol = 0.00001,
        ntwx = 10000, ntwe = 0, ntwr = 1000, ntpr = 1000,
        cut = 8.0, iwrap = 0,
        ntt = 3, gamma_ln = 1., ntb = 2, ntp = 1, barostat = 2
        nscm = 0,
        ioutfm = 1, ntxo = 2,
/
EOF

START="`date +%s.%N`"

# Minimization Phase
for RUN in 1min 2min 3min 4min ; do
 echo "------------------------"
 echo "Minimization phase: $RUN"
 echo "------------------------"
 if [[ ! -f $RUN.rst7 ]]; then
     echo "File -- $RUN.rst7 -- does not exists. Running job..."
     apptainer exec --nv --bind $bind_path $path_to_amber_container pmemd.cuda -O -i $RUN.in -p $prmtop_file -c $inpcrd_file -ref $inpcrd_file -o $RUN.out -x $RUN.nc -r $RUN.rst7 -inf $RUN.info
 else
     echo "File -- $RUN.rst7 -- exists.  Checking the next step."
 fi
 echo ""
 inpcrd_file="$RUN.rst7"
done

# Equilibration phase - reference coords are last coords from minimize phase
REF=$inpcrd_file
for RUN in 5md 6md 7md 8md 9md ; do

 echo "------------------------"
 echo "Equilibration phase: $RUN"
 echo "------------------------"
 if [[ ! -f $RUN.rst7 ]]; then
     echo "File -- $RUN.rst7 -- does not exists. Running job..."
     apptainer exec --nv --bind $bind_path $path_to_amber_container pmemd.cuda -O -i $RUN.in -p $prmtop_file -c $inpcrd_file -ref $REF -o $RUN.out -x $RUN.nc -r $RUN.rst7 -inf $RUN.mdinfo
 else
     echo "File -- $RUN.rst7 -- exists.  Checking the next step."
 fi
  echo ""
  inpcrd_file="$RUN.rst7"

done

sed -i 's/0.3000000E+02/0.0000000E+00/g' step9.rst7

STOP="`date +%s.%N`"
TIMING=`echo "scale=4; $STOP - $START;" | bc`
echo "$TIMING seconds."
echo ""

exit 0
