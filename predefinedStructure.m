% _________________________________________________________________________
%
% Title:        Predefined Structure Examples - BESOneo2
%                                                                      
% Version:      2.0 (14/06/2022)                                        
%                                                                      
% Author:       James Kirby, PhD Candidate, RMIT University             
% _________________________________________________________________________
%
% Description:  Several predefined structures
%               Custom definitions can be introduced by creating a named
%               selection in this function, or replacing section D in
%               BESOneo2.
%
% Inputs:       nx          number of elements along x
%               ny          number of elements along y
%               name        name of problem definition
%               fopt        set true is using fsparse function of Engblom, 
%                           S., & Lukarski, D. (2016)
%
% Outputs       F           force vector                                    sparse, dof indexed
%               free        free dofs                                       indices, dof indexed
%               fixed       fixed dofs                                      indices, dof indexed
%               act         active elements                                 indices, element indexed
%               pasV        passive void domain/s                           indices, element indexed
%               sym         symmetric boundary segments                     edge segments, node indexed
%               cyc         periodic boundary segments                      edge segments, node indexed
%               symX        Symmetry about horizontal at (ny+1)/2           true/false
%               symY        Symmetry about vertical at (nx+1)/2             true/false
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function [F, free, fixed, act, pasV, sym, cyc, symX, symY] = predefinedStructure(nx, ny, name, fopt)

    nEl     = nx*ny; 
    nDof    = 2*(ny+1)*(nx+1);
    nodeNrs = reshape(1:(1+nx)*(1+ny), 1+ny, 1+nx);
    elNrs   = reshape(1:nEl, ny, nx);
    
    switch name
        case "cantilever"
            %% Cantilever
            lcDof   = 2*nodeNrs(ny/2+1,nx+1);
            fixed   = 1:2*(ny+1);
            [pasS,pasV] = deal([],[]);
            if fopt
                F       = fsparse(lcDof',1,-1,[nDof,1]);
            else
                F       = sparse(lcDof',1,-1,nDof,1);
            end
            free    = setdiff(1:nDof,fixed);
            act     = setdiff((1:nEl)',union(pasS,pasV));
            sym     = [];
            cyc     = [];
            symX    = true;
            symY    = false;
        case "frameL"
            %% Frame reinforcement problem as presented in: 
            % "A new generation 99 line Matlab code for compliance topology 
            % optimization and its extension to 3D" by Ferrari and Sigmund, 2020, pp.11
            [lDofv,lDofh] = deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
            fixed   = [2*nodeNrs(ny+1,1)-1,2*nodeNrs(ny+1,1),2*nodeNrs(ny+1,nx+1)];
            a1      = elNrs(1:ny/50,:);
            a2      = elNrs(:,[1:nx/50,end-nx/50+1:end]);
            a3      = elNrs(2*ny/5:end,2*nx/5:end-nx/5);
            [pasS,pasV] = deal(unique([a1(:);a2(:)]),a3(:));
            if fopt
                F = fsparse(lDofv',1,-2/(nx+1),[nDof,1]) + ...
                fsparse(lDofh,1,-(1/ny:-1/(ny).^2:0)',[nDof,1]);
            else
                F = sparse(lDofv',1,-2/(nx+1),nDof,1) + ...
                sparse(lDofh,1,-(1/ny:-1/(ny).^2:0)',nDof,1);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV));  
            sym     = [];
            cyc     = [];
            symX    = false;
            symY    = false;
        case "frameR"
            %% Frame reinforcement problem as presented in: 
            % "A new generation 99 line Matlab code for compliance topology 
            % optimization and its extension to 3D" by Ferrari and Sigmund, 2020, pp.11
            [lDofv,lDofh] = deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
            fixed   = [2*nodeNrs(ny+1,1)-1,2*nodeNrs(ny+1,1),2*nodeNrs(ny+1,nx+1)];
            a1      = elNrs(1:ny/50,:);
            a2      = elNrs(:,[1:nx/50,end-nx/50+1:end]);
            a3      = elNrs(2*ny/5:end,2*nx/5:end-nx/5);
            [pasS,pasV] = deal(unique([a1(:);a2(:)]),a3(:));
            if fopt
                F = fsparse(lDofv',1,-2/(nx+1),[nDof,1]) + ...
                fsparse(lDofh,1,(1/ny:-1/(ny).^2:0)',[nDof,1]);
            else
                F = sparse(lDofv',1,-2/(nx+1),nDof,1) + ...
                sparse(lDofh,1,(1/ny:-1/(ny).^2:0)',nDof,1);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV));  
            sym     = [];
            cyc     = [];
            symX    = false;
            symY    = false;
        case "Lbracket"
            %%  L-bracket
            elNrs   = reshape(1:nEl, ny, nx);
            lDofh   = 2*nodeNrs(:,end)-1;
            fixed   = [1:2*(ny+1):nDof, 2:2*(ny+1):nDof];
            a3      = elNrs(1:2*ny/3,nx/2:end);
            [pasS,pasV] = deal([],a3(:));
            if fopt
                F = fsparse(lDofh(2*ceil(ny/3)+1)+1,1,-1,[nDof,1]);
            else
                F = sparse(lDofh(2*ceil(ny/3)+1)+1,1,-1,nDof,1);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV)); 
            sym     = [];
            cyc     = [];
            symX    = false;
            symY    = false;
        case "MBB"
            %% MBB beam
            fixed   = [1:2:2*(ny+1), nDof];
            [pasS,pasV] = deal([],[]);
            if fopt
                F = fsparse(2,1,-1,[nDof,1]);
            else
                F = sparse(2,1,-1,nDof,1);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV));   
            sym     = [nodeNrs(1:ny+1)'];
            cyc     = [];
            symX    = false;
            symY    = false;
        case "cantilever_ndr"
            %% Cantilever with circular non-design void
            lcDof   = 2*nodeNrs(ny/2+1,nx+1);
            fixed   = 1:2*(ny+1);
            pasS    = [];
            if fopt
                F       = fsparse(lcDof',1,-1,[nDof,1]);
            else
                F       = sparse(lcDof',1,-1,nDof,1);
            end
            free    = setdiff(1:nDof,fixed);
            [cx, cy, cr] = deal(nx/2, ny/2, ny/3);
            [dy,dx] = meshgrid(1:nx, 1:ny);
            ndr       = sqrt((dx-cy).^2+(dy-cx).^2) < cr;
            pasV    = find(ndr);
            act     = setdiff((1:nEl)',union(pasS,pasV));
            sym     = [];
            cyc     = [];
            symX    = true;
            symY    = false;
        case "custom"
            %% Add custom definitions here
            msgbox('Add custom structure definition in the predefinedStructure() function','Custom Structure');
    end % END switch case
end % END function
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\