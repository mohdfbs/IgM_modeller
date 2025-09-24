# IgM_modeller
This repository hosts structural templates, input files, and scripts for modelling full-length glycosylated IgM antibodies and preparing them for coarse-grained molecular dynamics (MD) simulations.

Overview
IgM antibodies, with up to 12 antigen-binding fragments (Fabs), offer significant advantages over IgG—including higher multivalency, avidity, and complement activation. Despite their promise in therapeutics and diagnostics, IgM-based antibody development has been hindered by the lack of full-length, high-resolution structural models.

This repository provides a modular computational pipeline to construct and simulate glycosylated IgM structures suitable for rational design and multiscale simulations.

Pipeline components:
1) Model building
   We provide an atomic resolution structure of Cetuximab IgM, which can be used as a template for modelling IgM for any antibody using the SWISS-MODEL webserver.
   a) On the SWISS-MODEL webserver, select the "User Template" mode.
   b) Paste the sequence of the antibody to be modelled based on the template provided (swissmodel_target_sequence.fasta). You only need to change the sequences of the heavy and light chain sequences, as the J-chain sequence in provided.
   c) Select "Add Template File" and upload the template structure provided (IgM-template.pdb). SWISS-MODEL will verify the structure.
   d) Once verified, add Project Title and Email, if needed.
   e) Finally, click on Build Model.
   f) Once finished, you should be able to download a PDB file of the IgM structure of your antibody.

2) Coarse-graining
   a) Use the provided martinize.py script as well as the dssp executable to convert your atomic IgM structure of a coarse-grained representation. Make sure to add elastic network and disulfide bonds appropriately. As there are disulfide bonds betwen the different chains of the IgM, you should also merge the topology file into a single topology file.
   Below is an example command:
      python martinize.py -f swissmodel_target_structure.pdb -o system.top -x CG_IgM.pdb -dssp ./dssp -p backbone -ff martini22 -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0 -elastic -name IgM -cys 0.5 -merge A,B,C,D,E,F,G,H,K,I,J,K,L,M,N,O,P,Q,R,S,T,U

3) Glycans addition
   a) We provide a library of common glycans (structure and topolgy files compatible with coarse-grained Martini 2.2 forcefield) and a python script to add these glycans to the coarse-grained structure of the IgM model. Go to the "Martini_glycosylator_v3" folder and open the README file for more info.
   b) Change the input.dat file to match the position and type of glycans for your antibody.
   c) Using GROMACS, renumber the residues to start with 1 and label the chain to A. This is necessary to match the input.dat file.
      gmx editconf -f CG_IgM.pdb -o CG_IgM_renumber.pdb -resnr 1 -label A
   c) Rename the topology file to reflect that all residues are now in chain A.
      cp  IgM_A+IgM_B+IgM_C+IgM_D+IgM_E+IgM_F+IgM_G+IgM_H+IgM_I+IgM_J+IgM_K+IgM_L+IgM_M+IgM_N+IgM_O+IgM_P+IgM_Q+IgM_R+IgM_S+IgM_T+IgM_U.itp CG_IgM_renumber_A.itp
   d) Use the python script to add glycans.
      python martini_glycosylator_v3.py -d 1.2 -c CG_IgM_renumber.pdb -f input.dat
      You might need to change the -d value to avoid steric clashes with surrounding protein residues.

4) Preparation for coarse-grained MD simulations
   a) Put the glycosylated IgM model into a box and add water as well as neutralising ions.
      gmx editconf -f CG_IgM_renumber_chain_A_gly.pdb -o IgM_box.gro -d 1 -c -bt cubic
      gmx solvate -cs water.gro -cp IgM_box.gro -o IgM_solv.gro -p system.top
      gmx grompp -f forcefield_files/ions.mdp -c IgM_solv.gro -p system.top -o IgM_ion -maxwarn 1
      echo 15 | gmx genion -s IgM_ion.tpr -o IgM_ion.gro -p system.top -pname NA+ -nname CL- -conc 0.15 -neutral
   b) Run energy minimisation.
      gmx grompp -f forcefield_files/minimization.mdp -c IgM_ion.gro -r IgM_ion.gro -p system.top -o minimization -maxwarn 1
      gmx mdrun -v -deffnm minimization
   c) Build an index file.
      gmx make_ndx -f minimization.gro <<EOF
      1|13|14
      name 17 proteingly
      15|16
      name 18 solv
      q
      EOF
   d) Run a multi-time step equilibrations.
      gmx grompp -f forcefield_files/equilibration_1fs.mdp -c minimization.gro -r minimization.gro -p system.top -o equilibration_1fs -n index.ndx -maxwarn 1
      gmx mdrun -v -deffnm equilibration_1fs

      gmx grompp -f forcefield_files/equilibration_5fs.mdp -c equilibration_1fs.gro -r equilibration_1fs.gro -p system.top -o equilibration_5fs -n index.ndx -maxwarn 1
      gmx mdrun -v -deffnm equilibration_5fs

      gmx grompp -f forcefield_files/equilibration_10fs.mdp -c equilibration_5fs.gro -r equilibration_5fs.gro -p system.top -o equilibration_10fs -n index.ndx -maxwarn 1
      gmx mdrun -v -deffnm equilibration_10fs

Now, you will have a system containing your IgM antibody that is ready for coarse-grained MD simulations.


   
