#!/bin/bash
# script to convert SMOG input files from GROMACS to LAMMPS
# L. Coronel, A. Suma and C. Micheletti
# v 1.0 May 9, 2020
#
#

echo "_____________________________________________________"
echo
echo "          GROMACS to LAMMPS conversion tool       "
echo "        for SMOG input scripts and data files     "
echo
echo "Version 1.0 - May 9, 2020 "
echo
echo "Please cite in publications as:"
echo "A. Suma, L. Coronel, G. Bussi and C. Micheletti"
echo "\"Directional translocation resistance of Zika xrRNA\"."
echo "_____________________________________________________"
echo
rootname=$1
echo "Processing files " ${rootname}.{gro,top}

# Read Box vectors
box_x=$(tail -1 ${rootname}.gro | awk '{print($1)}')
box_y=$(tail -1 ${rootname}.gro | awk '{print($2)}')
box_z=$(tail -1 ${rootname}.gro | awk '{print($3)}')


awk -v box_x=$box_x -v box_y=$box_y  -v box_z=$box_z 'BEGIN{k=0;k1=-1;nangles=-1;ndp=0;ndi=0}{
#atomtypes not used
#atoms
if($1=="Structure-Based")par=1
if($2=="defaults")par=0
if($2=="pairs")par=2
if($2=="bonds")par=3
if($2=="angles")par=4
if($2=="dihedrals")par=5
#read coordinates
if(par==1){
	if(FNR==2)natoms=$1
	if(NF==7 && FNR>=3){
		ind[$4]=$4
		x[$4]=$5
		y[$4]=$6
		z[$4]=$7
		#consider mol type
		#print ind[$4],x[$4],y[$4],z[$4]
	}
}
#read pairs
if(par==2){
	if(NF==5){
		k++
		pid1[k]=$1
		pid2[k]=$2
		qsc6[k] = $4
		qsc12[k] = $5
	}
	npair=k
}

#read bond
if(par==3){
	if(NF==5){
		k1++
		bid1[k1]=$1
		bid2[k1]=$2
		r0[k1]=$4
		kb[k1]=$5
	}
	nbonds=k1
}

#read angles
if(par==4){
	if(NF==6){
		nangles++
		aid1[nangles]=$1
		aid2[nangles]=$2
		aid3[nangles]=$3
		th0[nangles]=$5
		Ka[nangles]=$6
	}
}
#read dihedrals
if(par==5){
#LAMMPS proper dihedral angles
if(NF==8){
if($8==1){
ndp++
prop_id1[ndp]=$1
prop_id2[ndp]=$2
prop_id3[ndp]=$3
prop_id4[ndp]=$4
prop_phi1[ndp]=$6
prop_Kd[ndp]=$7

}
if($8==3){
prop_phi3[ndp]=$6
prop_Kd3[ndp]=$7

}
}
#LAMMPS improper dihedral angles
if(NF==7){
ndi++

iprop_id1[ndi]=$1
iprop_id2[ndi]=$2
iprop_id3[ndi]=$3
iprop_id4[ndi]=$4
chi0[ndi]=$6
Kchi[ndi]=$7
}
}


}END{
ENERGY_RESCALING_FACTOR=0.239;
DISTANCE_RESCALING_FACTOR=10.0;
DISTANCE_RESCALING_12=DISTANCE_RESCALING_FACTOR**12
DISTANCE_RESCALING_6=DISTANCE_RESCALING_FACTOR**6
print ""
print ""
print natoms,"atoms"
print nbonds,"bonds"
print nangles,"angles"
print ndp+ndi,"dihedrals"
print ""
print natoms,"atom types"
print nbonds,"bond types"
print nangles,"angle types"
print ndp+ndi,"dihedral types"
print ""
print 0, box_x*DISTANCE_RESCALING_FACTOR, " xlo xhi"
print 0, box_y*DISTANCE_RESCALING_FACTOR, " ylo yhi"
print 0, box_z*DISTANCE_RESCALING_FACTOR, " zlo zhi"
print ""
print "Masses"
print "# masses in g/mol (amu) - uniform default value = 16.00"
for(i=1;i<=natoms;i++)print i,"16.0"
print ""
print "Atoms"
print "#atom-ID molecule-ID atom-type q x y z"
for(i=1;i<=natoms;i++)print i,1,i,0,x[i]*DISTANCE_RESCALING_FACTOR,y[i]*DISTANCE_RESCALING_FACTOR,z[i]*DISTANCE_RESCALING_FACTOR;
print ""
print "Bond Coeffs"
print ""
# We are converting from Gromacs to LAMMPS, K has to be multiplied by 0.5
for(i=1;i<=nbonds;i++)print i,0.5*kb[i]*ENERGY_RESCALING_FACTOR/(DISTANCE_RESCALING_FACTOR*DISTANCE_RESCALING_FACTOR),r0[i]*DISTANCE_RESCALING_FACTOR; ## NF=3 bond_id, K, r0
print ""
print "Bonds"
print ""
for(i=1;i<=nbonds;i++)print i,i,bid1[i],bid2[i]
print ""
print "Angle Coeffs"
print ""
for(i=1;i<=nangles;i++)print i,0.5*Ka[i]*ENERGY_RESCALING_FACTOR,th0[i]
print ""
print "Angles"
print ""
for(i=1;i<=nangles;i++)print i,i,aid1[i],aid2[i],aid3[i]
print ""
print ""
print "Dihedral Coeffs"
print "#id num_dihedral_angles(2) K1 coeff1  phi1 K3 coeff3 phi3"
for(i=1;i<=ndp;i++)print i,2,prop_Kd[i]*ENERGY_RESCALING_FACTOR,1,prop_phi1[i],prop_Kd3[i]*ENERGY_RESCALING_FACTOR,3,prop_phi3[i]
for(i=1;i<=ndi;i++)print i+ndp,1,Kchi[i]*ENERGY_RESCALING_FACTOR,1,chi0[i]+180;
print ""
print "Dihedrals"
print ""
for(i=1;i<=ndp;i++)print i,i,prop_id1[i],prop_id2[i],prop_id3[i],prop_id4[i]
for(i=1;i<=ndi;i++)print i+ndp,i+ndp,iprop_id1[i],iprop_id2[i],iprop_id3[i],iprop_id4[i]

}' ${rootname}.gro ${rootname}.top  > data_file.lammps

echo "Created file data_file.lammps"

cat ${rootname}.top | awk  'BEGIN{k=0;}{
	#atomtypes not needed
	#atoms
	if($2=="defaults")par=0
	if($2=="pairs")par=2
	if($2=="bonds")par=3
	#read pairs
	if(par==2){
	if(NF==5){
		k++
		pid1[k]=$1
		pid2[k]=$2
		qsc6[k] = $4
		qsc12[k] = $5
		#sigmaij[k]=(2*qsc12/qsc6)**(1./6.)
		#epsc[k]=qsc6/2./(sigmaij[k]**6.)
	}
	npair=k
}
}END{
	ENERGY_RESCALING_FACTOR=0.239;
	DISTANCE_RESCALING_FACTOR=10.0;
	DISTANCE_RESCALING_12=DISTANCE_RESCALING_FACTOR**12
	DISTANCE_RESCALING_6=DISTANCE_RESCALING_FACTOR**6
print ""
print ""
print "units real"
print ""
print "boundary p p p"
print "#femtoseconds"
print "timestep 2"
print ""
print "neighbor    1.0 bin"
print "neigh_modify every 1 delay 0 check yes"
print ""
print "atom_style      full"
print "bond_style	harmonic"
print "angle_style	harmonic"
#print "dihedral_style	hybrid fourier quadratic" XXX commented out 9/5/2020
print "dihedral_style	fourier"
# print "improper_style  harmonic"
print ""
print "variable		cut_off equal 2.5 "
print "variable		ep_NC_gro equal 0.01"
print "variable		ep_unit_conv equal 0.239"
print "variable		ep_NC equal ${ep_NC_gro}*${ep_unit_conv}"
print "variable		C12 equal ${ep_NC}*${cut_off}^12"
print "pair_style	lj/smog  15"
print ""
print "read_data	data_file.lammps"
print "#read_restart restart_file.1000000"
print "restart 1000000 restart_file.*"
print ""
print "pair_coeff      * * ${C12} 0."
# write the interactions
for(i=1; i<=npair;i++) {#
	a=pid1[i];
	b=pid2[i];
	if(pid1[i]>pid2[i]) {a=pid2[i]; b=pid1[i];  }
	print "pair_coeff	",a,b,qsc12[i]*ENERGY_RESCALING_FACTOR*DISTANCE_RESCALING_12,qsc6[i]*ENERGY_RESCALING_FACTOR*DISTANCE_RESCALING_6;
}
print ""
print "variable	temp equal 95.0"
print ""
print "variable        rand equal floor(random(1,99999999,1))"
print "velocity        all create ${temp} ${rand} dist gaussian"
print ""
print "# fixes for Langevin simulations"
print "fix             1 all nve"
print "fix             2 all langevin ${temp} ${temp} 2000 ${rand}"
print ""
print "# Remove center of mass displacements"
print "#fix 		3 all momentum 100 linear 1 1 1 angular"
print ""
print "dump            1 all custom 1000 dump.lammpstrj id type x y z"
print "dump_modify 1 sort id   append yes"
print "dump_modify     1 sort id"
print "#thermo_style	custom step temp etotal"
print "thermo 		10000"
print ""
print "#minimize      0.0 1.0e-8 1000 100000"
print "run		1000000"
}' > input_script.lammps
echo "Created file input_script.lammps"
echo "Done."
echo
