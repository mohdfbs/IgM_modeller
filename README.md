# IgM Modeller

This repository hosts structural templates, input files, and scripts for modelling **full-length glycosylated IgM antibodies** and preparing them for **coarse-grained molecular dynamics (MD) simulations**.

---

## 🔬 Overview

IgM antibodies, with up to **12 antigen-binding fragments (Fabs)**, offer significant advantages over IgG—including higher multivalency, avidity, and complement activation. Despite their promise, IgM-based antibody development has been hindered by the **lack of full-length, high-resolution structural models**.  

This repository provides a **modular computational pipeline** to construct and simulate glycosylated IgM structures suitable for **rational design and multiscale simulations**.

---

# 📌 Pipeline

We provide the workflow in **four main steps**:  

- [Step 1: Model Building](#step-1-model-building)  
- [Step 2: Coarse-Graining](#step-2-coarse-graining)  
- [Step 3: Glycan Addition](#step-3-glycan-addition)  
- [Step 4: Preparation for MD Simulations](#step-4-preparation-for-md-simulations)  

---

<details>
<summary>🔹 Step 1: Model Building</summary>

We provide an **atomic resolution Cetuximab IgM template** for modelling.  

1. Go to [SWISS-MODEL](https://swissmodel.expasy.org) → select **User Template** mode.  
2. Paste the antibody sequence based on `swissmodel_target_sequence.fasta`.  
   - Only modify heavy & light chains.  
   - J-chain is already provided.  
3. Upload template structure: `IgM-template.pdb`.  
4. (Optional) Add Project Title & Email.  
5. Click **Build Model**.  
6. Download the final **PDB IgM model**.  

</details>

---

<details>
<summary>🔹 Step 2: Coarse-Graining</summary>

Use **martinize.py** and **DSSP** to convert your atomic model into a CG representation.  

💻 Example command:
```bash
python martinize.py \
   -f swissmodel_target_structure.pdb \
   -o system.top \
   -x CG_IgM.pdb \
   -dssp ./dssp \
   -p backbone \
   -ff martini22 \
   -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0 \
   -elastic \
   -name IgM \
   -cys 0.5 \
   -merge A,B,C,D,E,F,G,H,K,I,J,K,L,M,N,O,P,Q,R,S,T,U
```

⚠️ Remember to:  
- Add **elastic network** & **disulfide bonds** appropriately.  
- Merge topologies across all IgM chains into a **single topology file**.  

</details>

---

<details>
<summary>🔹 Step 3: Glycan Addition</summary>

We provide a **Martini-compatible glycan library** and `martini_glycosylator_v3.py`.  

1. Edit `input.dat` to define glycan types & positions.  
2. Renumber residues & chain ID with GROMACS:  
   ```bash
   gmx editconf -f CG_IgM.pdb -o CG_IgM_renumber.pdb -resnr 1 -label A
   ```  
3. Rename topology file:  
   ```bash
   cp IgM_A+...+IgM_U.itp CG_IgM_renumber_A.itp
   ```  
4. Add glycans with script:  
   ```bash
   python martini_glycosylator_v3.py \
      -d 1.2 \
      -c CG_IgM_renumber.pdb \
      -f input.dat
   ```  

👉 Adjust the `-d` value to avoid steric clashes.  

</details>

---

<details>
<summary>🔹 Step 4: Preparation for MD Simulations</summary>

1. **Solvate & Ionize**  
   ```bash
   gmx editconf -f CG_IgM_renumber_chain_A_gly.pdb -o IgM_box.gro -d 1 -c -bt cubic
   gmx solvate -cs water.gro -cp IgM_box.gro -o IgM_solv.gro -p system.top
   gmx grompp -f forcefield_files/ions.mdp -c IgM_solv.gro -p system.top -o IgM_ion -maxwarn 1
   echo 15 | gmx genion -s IgM_ion.tpr -o IgM_ion.gro -p system.top -pname NA+ -nname CL- -conc 0.15 -neutral
   ```

2. **Energy Minimization**  
   ```bash
   gmx grompp -f forcefield_files/minimization.mdp -c IgM_ion.gro -r IgM_ion.gro -p system.top -o minimization -maxwarn 1
   gmx mdrun -v -deffnm minimization
   ```

3. **Index File Creation**  
   ```bash
   gmx make_ndx -f minimization.gro <<EOF
   1|13|14
   name 17 proteingly
   15|16
   name 18 solv
   q
   EOF
   ```

4. **Equilibrations**  
   ```bash
   gmx grompp -f forcefield_files/equilibration_1fs.mdp -c minimization.gro -r minimization.gro -p system.top -o equilibration_1fs -n index.ndx -maxwarn 1
   gmx mdrun -v -deffnm equilibration_1fs

   gmx grompp -f forcefield_files/equilibration_5fs.mdp -c equilibration_1fs.gro -r equilibration_1fs.gro -p system.top -o equilibration_5fs -n index.ndx -maxwarn 1
   gmx mdrun -v -deffnm equilibration_5fs

   gmx grompp -f forcefield_files/equilibration_10fs.mdp -c equilibration_5fs.gro -r equilibration_5fs.gro -p system.top -o equilibration_10fs -n index.ndx -maxwarn 1
   gmx mdrun -v -deffnm equilibration_10fs
   ```

✅ You now have a **glycosylated IgM system** ready for CG MD simulations.  

</details>

---
