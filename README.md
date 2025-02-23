peptideAntibodyCompetition.m  
A MATLAB function calculating the equilibrium binding concentrations of two binders (B1, B2)
competitively binding to a single target (T) with two distinct dissociation constants (Kd1, Kd2).
Their complexes are designated by B1T and B2T, respectively.
The total concentrations of B1, B2 and T are designated by B1tot, B2tot and Ttot, respectively.
Since the solved equation set is a cubic one, three solutions are calculated, and the last
output of the function, goodSolutionIndices, contains the index of the meaningful solution.

fitManyAniSatCurves.m  
The program can be used for analyzing anisotropy measurements in which a the binding of a small fluorescent binder to a large target is measured. The program determines the min and max anisotropies and the dissociation constants for each curve by global fitting.  
anidata - a structure array with each element having the following fields:  
        - data: an n x 2 array, 1st column - concentration, 2nd column: anisotropy  
        - id: name of the sample  
        - remark: what its name suggests  
elementForBoundAni - this dataset will be used for finding the max anisotropy, i.e., the anisotropy of the bound peptide (if other datasets don't reach saturation). If empty, anisotropy of the free and bound peptide will also be fit globally.  
ccPeptide - concentration of the fluorescent peptide that binds to the antibody


**Functions used by fitManyAniSatCurves.m**

fitBindingAnisotropy.m  
The program fits the anisotropy increase of a small peptide as it binds to a large target site.  
cFluorescentPeptide - concentration of the peptide whose anisotropy is measured  
anisotropy - anisotropy to be fitted  
cLargeBinder - concentration of the target the peptide binds to  

simpleBindingWithDepletion.m  
The program calculates a binding equilibrium assuming ligand depletion between a ligand (L) binding to a target (T) with a dissociation constant of Kd, and returns the concentrations of free target, free ligand and the target-ligand complex.

