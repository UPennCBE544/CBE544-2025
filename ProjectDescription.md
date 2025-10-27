---
layout: page
mathjax: true
permalink: /Project/
---
## Course Project Logistics ##

1. [Introduction](#intro)
2. [Motivation](#MO)
3. [Terms to Know](#Terms)
4. [Plan](#Plan)
5. [Final Presentation Rubric](https://github.com/UPennCBE544/CBE5440-2024/blob/main/Final_Project/CBE%205440%20-%20Final%20Project%20-%20Presentation%20Rubric%20and%20Expectations%20Sheet.pdf)
6. [Final Report Rubric](https://github.com/UPennCBE544/CBE5440-2024/blob/main/Final_Project/CBE%205440%20-%20Final%20Project%20-%20Report.pdf)

Turn in your final report on CANVAS or email to:

```
alevoj@seas.upenn.edu, rthatch@seas.upenn.edu
```

<a name='intro'></a>

## Project Introduction ##
Material: MXenes are an emerging family of two-dimensional transition metal carbides and nitrides with the general formula M<sub>n+1</sub>X<sub>n</sub>T<sub>x</sub> (n = 1-4), where M is an early transition metal (Sc, Ti, Zr, V, Nb, Ta, etc), X is C or N, and T is a surface group such as O, F, Cl, and OH. Its main properties are: scalable synthesis, solution processability, large surface-to-volume ratio, high metallic conductivity, and mechanical strength. Due to these properties, applications for energy storage, composites, and optoelectronics are being explored.

Goal: The main scientific goal is to design a potential catalytic reaction from functionalzied Nb<sub>2</sub>C MXenes. First, we will look at the Cl-terminated Nb<sub>2</sub>C MXenes and adsorb cycloalkane, cycloalkene, cycloalkyne from 1~10 carbon loop and find a thermodynamically optimal adsorption configuration by calculating the adsorption energy. And we will run DOS, bader charge calculations to analyze the system and probe into the electronic/thermodynamic properties of the system to design a possible pathway of catalytic reaction. 

<a name='MO'></a>

### Motivation ###

- **Synthesis of 2D Organometallic Material**: MXene is a 2 dimentional transitional metal carbide or nitride which allows the material to be easily tunable with surface modification to have drastically  different properties. A deeper understanding of thermodynamic stability and adsorption mechanism supports the advancement of sustainable and ecnomically attractive option for active electrocatalysts.

- **Adsorption Configuration and Optimization**: By modeling different surface configuration of the adsorption, we gain insight into how environmental conditions (vacancy, dehydrogenation, temperature, etc.) and configuration influence the adsorption stability. This could lead to optimization strategies that improve catalyst activity and stability by tuning vacancy rates and ligand environment.

- **Catalyst Design for Functional Efficiency**: Selecting catalysts with the right functional properties is crucial for enabling manufacturable technologies. This includes considerations for cost efficiency, performance, durability, raw material abundance, and environmental impact. This project will specifically lend itself toward oulining a framework for studying activity and stability of functionalized MXenes. 

- **Experimental Comparison and Data Validation**: By aligning the surface configuration, adsorption, and catalytic behaviors observed in simulations with those measured experimentally, we can validate our models and refine predictive accuracy. This comparison helps identify discrepancies, improve our understanding of real-world conditions, and build confidence in the applicability of the computational framework for guiding experimental design and optimizing catalyst performance.

<a name='Terms'></a>

## Terminology ##
- **Surface Facet Miller Index**: The Miller index (hkl) describes the orientation of a crystal plane by indicating how it cuts through the unit cell of the crystal structure

- **Termination**: Termination refers to the specific arrangement and coordination of lattice metals and their corresponding oxygens at a surface. Different terminations, whether metal-terminated or oxygen-terminated, significantly affect surface properties, including reactivity and stability, and can result from various preparation methods.

- **Symmetric Surface**: A symmetric surface is characterized by an arrangement of atoms that is the same on both sides of the surface plane, exhibiting inversion symmetry. They may not always be stoichiometric, meaning the atomic ratios (e.g., metal to oxygen) can differ from the bulk composition.

- **Asymmetric Surface**: An asymmetric surface tends to lack inversion symmetry of atoms on both sides, but not always. We will define the clean (lacking adsorbates) asymmetric surfaces to be stoichiometric, maintaining a specific metal-to-oxygen ratio overall.

- **Pourbaix Diagram**: Also known as a thermodynamic stability diagram. It plots free energy relative to some reference surface (usually the bulk cleaved configuration) vs some variable, in our case we will plot potential ($U$) on the x axis. Free energy is calcuated for each surface via the following equation. $\Delta G = E_{\text{slab}}^{\text{DFT}} - E_{\text{reference}}^{\text{DFT}} + \sum_A \left( \Delta n_A \cdot (\text{ref}_A^{\text{DFT}} + \text{correction}_A(T, P, \mu, pH, U)) \right)$ where A is any element in our system, the correction term is a function of selected environmental conditions.

- **Surface Energy ($\gamma_{(hkl)}$)**: This can be thought of as the energy for the specific modeled surface to exist. It indicates the thermodynamic stability of a surface. Surface energy is influenced by surface geometry, termination, and adsorbates. $\gamma_{(hkl)}(pH,U) = -\frac{E_{slab} + N \cdot E_{bulk} - \sum_{i} \delta n_{i} \cdot \mu_{i}^{\text{REF}}(pH,U)}{2A}$. Note: $E_{bulk}$ is calculated by generating 4 of the same surface facet+termination each with different amounts of "bulk" and plotting total potential energy vs system size and getting the slope.

- **Wulff's Theorem**: This theorem states that the shape of a crystal/nanoparticle at equilibrium is determined by the minimization of its total surface energy. The equilibrium shape (also known as the **Wulff Construction**) is the convex hull of the surfaces with the lowest energies.

<a name='Plan'></a>

## Plan ##

For each symmetric Surface we need to go through the following path.

<img width="968" alt="Screenshot 2024-10-22 at 8 00 01 PM" src="https://github.com/user-attachments/assets/8ef50116-698b-476e-b0db-ed8eb3596ee4">

### Detailed plan: ###

We will break into groups of 2 students and each group will be assigned a series of surfaces they will be responsible for.
   
       a. {(2,2,1),(3,3,1),(4,4,1),(0,1,2)} - Group 1

       b. {(1,2,2),(1,3,3),(1,4,4),(0,2,1)} - Group 2

       c. {(2,1,1),(3,1,1),(4,1,1),(2,1,0)} - Group 3

       d. {(1,1,5),(1,5,5),(1,1,2),(1,1,3)} - Group 4

       Groups: 
         (1) Conor & Erika
         (2) Daniel & George
         (3) Shiqiang & Jinyu
         (4) Yeri & Achala

      Nanoparticle Conditions: 
         Group 1 - pH = 1, U = 0
         Group 2 - pH = 1, U = 1.3
         Group 3 - pH = 1, U = 1.6


<img width="274" alt="Screenshot 2024-11-20 at 2 58 56 PM" src="https://github.com/user-attachments/assets/527de129-a4fa-4f4b-b8a6-6bfef077d897">

0. Prep Work - we need a rutile RuO<sub>2</sub> bulk that is optimized

2. Generate Asymmetric Surface Facet - Provided by Rachel

3. Adsorb on the Asymmetric Surface at up to 3 adsorbate coverages -
      
      FIRST. Please go into your personal scratch and delete everything: ``` rm -r /anvil/scratch/x-YOURUSERNAME/Final_Project_2024 ```
   
      a. Update the Adsorption script in your home:

         - First remove the old scripts folder: rm -r ~/scripts_Final_Project

         - Copy over the adsorption script to your home: cp -r /anvil/projects/x-eve210010/scripts/scripts_Final_Project ~/
   
      b. Open the reference trajectory file: ```ag /anvil/projects/x-eve210010/REFERENCES/dopedSurface/ruo2/YOUR_FACET/clean/No_defect/0%_doped/PBE/relax/init.traj```

         In terminal: control+z, then type bg
   
      c. Determine which metal/oxygen atom/s you want to adsorb onto. If you have 4 or more metal atoms at the surface, you will need to adsorb on 1, approximately half, and all of the surface sites. We tend to adsorb on metal atoms, however, if there is an oxygen on top of the metal atom, we will instead adsorb on that oxygen.
   
      d. Modify the Adsorption file (in the last few lines) to have YOUR_FACET and the indices of the atoms you want to adsorb on top of: ```nano ~/scripts_Final_Project/asymmetric_slabs/Adsorptions_Asymmetric.py ```
      <img width="652" alt="Screenshot 2024-11-20 at 1 29 57 PM" src="https://github.com/user-attachments/assets/a0329922-7904-48bc-9184-1c28667a6898">

      e. Run the Adsorption python script: ```python ~/scripts_Final_Project/asymmetric_slabs/Adsorptions_Asymmetric.py ```

      f. Check on the ```init.traj``` files that were generated and listed after running the Adsorption python script. Do this by opening the files using ```ag FILE_PATH```

      g. **Do the other coverages that need to be generated (reminder: we want to do 1 adsorbate, approximately half adsorbed, and all adsorbed).** The grey and white table above shows the max number of sites you can adsorb onto.
   
      h. **Repeate b-f for your second facet.**

      i. Now copy over clean/No_defect surface from shared scratch to your directory structure for BOTH facets:

         1. Make sure to specify YOUR_USERNAME and your YOUR_FACET in this line: mkdir -p /anvil/scratch/YOUR_USERNAME/Final_Project_2024/dopedSurface/ruo2/YOUR_FACET/clean/No_defect/0%_doped/PBE
   
         2. Copy the init.traj from the shared scratch (again make sure to use YOUR_FACET and YOUR_USERNAME): cp /anvil/projects/x-eve210010/REFERENCES/dopedSurface/ruo2/YOUR_FACET/clean/No_defect/0%_doped/PBE/relax/init.traj /anvil/scratch/YOUR_USERNAME/Final_Project_2024/dopedSurface/ruo2/YOUR_FACET/clean/No_defect/0%_doped/PBE
   
      i. Submit the jobs: (jobs submit one every 10 seconds, so be patient)

            Run this to confirm the directories you want to submit from are listed: python ~/scripts_Final_Project/asymmetric_slabs/submit_asymmetric.py

            Modify your own submission script so that submit = True: nano ~/scripts_Final_Project/asymmetric_slabs/submit_asymmetric.py

            Now run the script again: python ~/scripts_Final_Project/asymmetric_slabs/submit_asymmetric.py

            After submission, modify your own submission script so that submit = False: nano ~/scripts_Final_Project/asymmetric_slabs/submit_asymmetric.py

      j. As jobs are completing you can run this python script that will make ```opt.traj``` files out of finished jobs: ```python /anvil/projects/x-eve210010/scripts/scripts_Final_Project/asymmetric_slabs/LOGTRAJ.py``` 
   
5. Generate Pourbaix Diagram -
   
      a. Modify the Pourbaix Script to use your own FACET: ```nano ~/scripts_Final_Project/asymmetric_slabs/Pourbaix/Pourbaix_Electrochemical.py```
   
      b. Run the script: ```python ~/scripts_Final_Project/asymmetric_slabs/Pourbaix/Pourbaix_Electrochemical.py```

      c. Now a pourbaix diagram was generated and saved in ```/anvil/projects/x-eve210010/scripts/scripts_Final_Project/asymmetric_slabs/Pourbaix/Pourbaix_Diagrams/full_diagrams/pH=1``` and the legend is coped into this folder for better viewing ```/anvil/projects/x-eve210010/scripts/scripts_Final_Project/asymmetric_slabs/Pourbaix/Pourbaix_Diagrams/Legend```

      d. We can pull the most stable configurations of the surface: ```python /anvil/projects/x-eve210010/scripts/scripts_Final_Project/asymmetric_slabs/Pourbaix/pull_most_stable.py YOUR_FACET```

   
6. Map adsorbates to the symmetric facet

   a. Make a ```sym_surf``` directory: ``` mkdir -p /anvil/scratch/YOUR_USERNAME/Final_Project_2024/dopedSurface/ruo2/sym_surf/ ```
    
   b. Copy clean symmetric surfaces for both facets into your own directory structure: ```cp -r /anvil/projects/x-eve210010/REFERENCES/dopedSurface/ruo2/sym_surf/YOUR_FACET /anvil/scratch/YOUR_USERNAME/Final_Project_2024/dopedSurface/ruo2/sym_surf/ ```

   c. Copy this folder to your scripts:

         1. mv ~/scripts_Final_Project/symmetric_slabs ~/scripts_Final_Project/symmetric_slabs_old
   
         2. cp -r /anvil/projects/x-eve210010/scripts/scripts_Final_Project/symmetric_slabs/ ~/scripts_Final_Project

   d. Modify this script to have your facet and your adsorbate configuration you want to map:
   
         1. Make modifications: nano ~/scripts_Final_Project/symmetric_slabs/adsorbate_map.py

         2. Run the script: python ~/scripts_Final_Project/symmetric_slabs/adsorbate_map.py

   e. Do this for the rest of the adsorbate configurations you want to map

   f. Submit the jobs: (jobs submit one every 10 seconds, so be patient)

            Run this to confirm the directories you want to submit from are listed: python ~/scripts_Final_Project/symmetric_slabs/submit_symmetric.py

            Modify your own submission script so that submit = True: nano ~/scripts_Final_Project/symmetric_slabs/submit_symmetric.py

            Now run the script again: python ~/scripts_Final_Project/symmetric_slabs/submit_symmetric.py

            After submission, modify your own submission script so that submit = False: nano ~/scripts_Final_Project/symmetric_slabs/submit_symmetric.py

   g.  As jobs are completing you can run this python script that will make ```opt.traj``` files out of finished jobs: ```python /anvil/projects/x-eve210010/scripts/scripts_Final_Project/asymmetric_slabs/LOGTRAJ.py``` 

