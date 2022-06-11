function [F, free, act, pasV] = predefinedStructure(nx, ny, name)

    nEl     = nx*ny; 
    nDof    = 2*(ny+1)*(nx+1);
    nodeNrs = int32(reshape(1:(1+nx)*(1+ny), 1+ny, 1+nx));
    elNrs   = reshape(1:nEl, ny, nx);
    list    = dir('**/*.*');
    if any(contains({list.name}, 'fsparse')==1)
        fopt = true;
    else
        fopt = false;
    end
    
    switch name
        case "cantilever"
            %% Cantilever
            lcDof   = 2*nodeNrs(ny/2+1,nx+1);
            fixed   = 1:2*(ny+1);
            [pasS,pasV] = deal([],[]);
            if fopt
                F       = fsparse(lcDof',1,-1,[nDof,1]);
            else
                F       = sparse(lcDof',1,-1,[nDof,1]);
            end
            free    = setdiff(1:nDof,fixed);
            act     = setdiff((1:nEl)',union(pasS,pasV));
        case "frame"
            %% Frame reinforcement problem as presented in: 
            % "A new generation 99 line Matlab code for compliance topology 
            % optimization and its extension to 3D" by Ferrari and Sigmund, 2020, pp.11
            [lDofv,lDofh] = deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
            fixed   = [1,2,nDof];
            a1      = elNrs(1:ny/50,:);
            a2      = elNrs(:,[1:nx/50,end-nx/50+1:end]);
            a3      = elNrs(2*ny/5:end,2*nx/5:end-nx/5);
            [pasS,pasV] = deal(unique([a1(:);a2(:)]),a3(:));
            if fopt
                F = fsparse(lDofv',1,-2/(nx+1),[nDof,1]) + ...
                fsparse(lDofh,1,[0:1/(ny).^2:1/ny]',[nDof,1]);
            else
                F = sparse(lDofv',1,-2/(nx+1),[nDof,1]) + ...
                sparse(lDofh,1,[0:1/(ny).^2:1/ny]',[nDof,1]);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV));  
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
                F = sparse(lDofh(2*ceil(ny/3)+1)+1,1,-1,[nDof,1]);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV)); 
        case "MBB"
            %% MBB beam
            fixed   = [1:2:2*(ny+1), nDof];
            [pasS,pasV] = deal([],[]);
            if fopt
                F = fsparse(2,1,-1,[nDof,1]);
            else
                F = sparse(2,1,-1,[nDof,1]);
            end
            free    = setdiff(1:nDof,fixed);                                      
            act     = setdiff((1:nEl)',union(pasS,pasV));   
        case "cantilever_ndr"
            %% Cantilever with circular non-design void
            lcDof   = 2*nodeNrs(ny/2+1,nx+1);
            fixed   = 1:2*(ny+1);
            pasS    = [];
            if fopt
                F       = fsparse(lcDof',1,-1,[nDof,1]);
            else
                F       = sparse(lcDof',1,-1,[nDof,1]);
            end
            free    = setdiff(1:nDof,fixed);
            [cx, cy, cr] = deal(nx/2, ny/2, ny/3);
            [dy,dx] = meshgrid(1:nx, 1:ny);
            f       = sqrt((dx-cy).^2+(dy-cx).^2) < cr;
            pasV    = find(f);
            act     = setdiff((1:nEl)',union(pasS,pasV));
    end
end