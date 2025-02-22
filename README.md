peptideAntibodyCompetition.m
A MATLAB function calculating the equilibrium binding concentrations of two binders (B1, B2)
competitively binding to a single target (T) with two distinct dissociation constants (Kd1, Kd2).
Their complexes are designated by B1T and B2T, respectively.
The total concentrations of B1, B2 and T are designated by B1tot, B2tot and Ttot, respectively.
Since the solved equation set is a cubic one, three solutions are calculated, and the last
output of the function, goodSolutionIndices, contains the index of the meaningful solution.
