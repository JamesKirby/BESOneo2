% _________________________________________________________________________
%
% Title:        Apply Symmetry Constraints - BESOneo2
%                                                                      
% Version:      0.1 (14/06/2022)                                        
%                                                                      
% Author:       James Kirby, PhD Candidate, RMIT University             
% _________________________________________________________________________
%
% Description:  Enforce symmetry constraints
%
% Inputs:       G           input matrix
%               symX        symmetry constraint x (true/false)
%               symY        symmetry constraint y (true/false)
%
% Outputs       G           output matrix
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function [G] = applySymmetryConstraints(G, symX, symY)
    if symX
        Gnorm = flipud(G);
        G = max(G, Gnorm);
    end
    if symY
        Gnorm = fliplr(G);
        G = max(G, Gnorm);
    end
end % END function
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\