%% Title:       BESOneo2 - Efficient Minimum Compliance BESO (2D)
% Version:      1.0 (04/11/2020)
% Author:       James Kirby, PhD Candidate, RMIT University
%
% Description:  BESO based on the efficient FEA implementation as presented 
%               in "A new generation 99 line Matlab code for compliance 
%               topology optimization and its extension to 3D" by Ferrari 
%               and Sigmund, 2020.
%
%               This implementation benefits from improved FEA
%               implementation including efficient stiffness matrix
%               construction, improved sparse matrix construction using
%               fsparse, and MATLABs decomposition function with cholsky
%               factorization. 
%
%               Load and boundary specification can be modified in section
%               C and passive void and solid elements (non-designable
%               regions) can be set using pasS and pasV.
%
% Notes:        This implementation uses the fsparse function as described
%               in Engblom, S., & Lukarski, D. (2016). Fast Matlab 
%               compatible sparse assembly on multicore computers. Parallel 
%               Computing, 56, 1�17. 
%               https://doi.org/10.1016/j.parco.2016.04.001. The fsparse 
%               code is available from the first authors GitHub page at: 
%               https://github.com/stefanengblom/stenglib. 
%               ___________________________________________________________
%
function [xfinal, obj] = BESOneo2(nx,ny,volfrac,er,rmin)
    bcF     = 'symmetric';                                                  % imfilter option (zero-Neumann filter boundary)
    penal   = 3;                                                            % penalty exponent for material interpolation
    xmin    = 1e-3;                                                         % material density of void elements
    maxit   = 200;                                                          % maximum number of iterations allowed
% A) Material paramaters
    E0      = 1;                                                            % Young modulus of solid
    Emin    = 1e-6;                                                         % Young modulus of "void"
    nu      = 0.3;                                                          % Poisson ratio
% B) Discretization features
    nEl     = nx*ny;                                                        % number of elements
    nodeNrs = int32(reshape(1:(1+nx)*(1+ny), 1+ny, 1+nx));                  % indexed grid node numbers (defined as int32)
    cVec    = reshape(2*nodeNrs(1:end-1, 1:end-1)+1, nEl, 1);
    cMat    = cVec+int32([0, 1, 2*ny+[2, 3, 0, 1], -2, -1]);                % connectivity matrix
    nDof    = 2*(ny+1)*(nx+1);                                              % total number of DOFs
    [sI,sII]= deal([]);
    for j = 1:8
        sI  = cat(2, sI, j:8);
        sII = cat(2, sII, repmat(j, 1, 8-j+1));
    end
    [iK,jK] = deal(cMat(:, sI)', cMat(:, sII)');
    Iar     = sort([iK(:), jK(:)], 2, 'descend'); clear iK jK               % reduced assembly indexing
    c1      = [12;3;-6;-3;-6;-3;0;3;12;3;0;-3;-6;-3;-6;12;-3;0;-3;-6;3;...
               12;3;-6;3;-6;12;3;-6;-3;12;3;0;12;-3;12];
    c2      = [-4;3;-2;9;2;-3;4;-9;-4;-9;4;-3;2;9;-2;-4;-3;4;9;2;3;-4;...
               -9;-2;3;2;-4;3;-2;9;-4;-9;4;-4;-3;-4];
    Ke      = 1/(1-nu^2)/24*(c1+nu.*c2);                                    % lower triangular portion of element K matrix
    Ke0(tril(ones(8))==1) = Ke';
    Ke0     = reshape(Ke0,8,8);
    Ke0     = Ke0+Ke0'-diag(diag(Ke0));                                     % recover full element K matrix
% C) Loads, supports & passive domains
    elNrs   = reshape(1:nEl, ny, nx);                                       % simplifies specification of pasS and pasV
    lcDof   = 2*nodeNrs(ny/2+1,nx+1);                                       % DOFs with applied load
    fixed   = 1:2*(ny+1);                                                   % fixed DOFs
    [pasS,pasV] = deal([],[]);                                              % passive solid and void elements
    F       = fsparse(lcDof',1,-1,[nDof,1]);                                % define load vector
    free    = setdiff(1:nDof,fixed);                                        % set of free DOFs
    act     = setdiff((1:nEl)',union(pasS,pasV));                           % set of active elements
% D) Prepare filter
    [dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
    h       = max(0, rmin-sqrt(dx.^2+dy.^2));                               % discretized gaussian filter kernel
% E) Allocate & initialize other parameters
    [x,dsK] = deal(ones(nEl, 1), zeros(nEl, 1));                            % initialize vectors
    x(pasV) = xmin;                                                         % set x = xmin on pasV set
    [ch, loop, vol, U] = deal(1, 0, sum(x,'all')/nEl, zeros(nDof,1));
    while ch > 1e-6 && loop<maxit
        loop    = loop+1;                                                   % update loop counter
        vol     = max(vol*(1-er), volfrac);
        if loop > 1; olddc = dc; end
% 1) Setup & solve equilibrium equations
        sK      = Emin+x.^penal*(E0-Emin);                                  % stiffness interpolation
        dsK(act)= penal.*x(act).^(penal-1);                                 % derivative of stiffness interpolation
        sK      = reshape(Ke(:)*sK', length(Ke)*nEl, 1);
        K       = fsparse(Iar(:, 1), Iar(:, 2), sK, [nDof, nDof]);          % assemble stiffness matrix
        U(free) = decomposition(K(free, free), 'chol', 'lower')\F(free);    % solve equilibrium system
        obj(loop)= 0.5*F'*U;
% 2) Compute Sensitivities
        dc      = dsK.*sum((U(cMat)*Ke0).*U(cMat), 2);                      % derivative of compliance
        dc      = imfilter(reshape(dc, ny, nx), h, bcF, 'conv');            % filter objective sensitivity
        dc      = dc(:);
        if loop>1; dc = mean([dc,olddc], 2); end
% 3) BESO design update
        l = [min(dc), max(dc)];                                             % initial estimate for lambda
        while (l(2)-l(1))/l(2) > 1e-5
            th = (l(1)+l(2))/2.0;
            x(act) = max(1e-3, sign(dc(act)-th));
            if (sum(x) -nEl*vol) > 0; l(1) = th; else, l( 2 ) = th; end
        end
% 4) Print current results and plot design
        if loop > 10
            ch = abs(sum(obj(loop-9:loop-5))-sum(obj(loop-4:loop)))/sum(obj(loop-4:loop));
        end
        fprintf( 'It.:%7i Obj:%7.4f Vol:%7.3f ch.:%0.2e \n', ...
            loop, obj(loop), vol, ch); 
        colormap(gray); imagesc(reshape(x<1, ny, nx));
        axis equal tight off; drawnow;
    end
    xfinal = reshape(x, ny, nx);
end    

%% Cantilever example
%
%     lcDof   = 2*nodeNrs(ny/2+1,nx+1);
%     fixed   = 1:2*(ny+1);
%     [pasS,pasV] = deal([],[]);
%     F       = fsparse(lcDof',1,-1,[nDof,1]);
%     free    = setdiff(1:nDof,fixed);
%     act     = setdiff((1:nEl)',union(pasS,pasV));

%% Frame reinforcement problem as presented in: 
% "A new generation 99 line Matlab code for compliance topology 
% optimization and its extension to 3D" by Ferrari and Sigmund, 2020, pp.11
%
%     elNrs   = reshape(1:nEl, ny, nx);
%     [lDofv,lDofh] = deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
%     fixed   = [1,2,nDof];
%     a1      = elNrs(1:ny/50,:);
%     a2      = elNrs(:,[1:nx/50,end-nx/50+1:end]);
%     a3      = elNrs(2*ny/5:end,2*nx/5:end-nx/5);
%     [pasS,pasV] = deal(unique([a1(:);a2(:)]),a3(:));
%     F = fsparse(lDofv',1,-2/(nx+1),[nDof,1]) + ...
%         fsparse(lDofh,1,[0:1/(ny).^2:1/ny]',[nDof,1]);
%     free    = setdiff(1:nDof,fixed);                                      
%     act     = setdiff((1:nEl)',union(pasS,pasV));                        

%%  L-bracket example
%
%     elNrs   = reshape(1:nEl, ny, nx);
%     [lDofv,lDofh] = deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
%     fixed   = [1:2*(ny+1):nDof, 2:2*(ny+1):nDof];
%     a3      = elNrs(1:2*ny/3,nx/2:end);
%     [pasS,pasV] = deal([],a3(:));
%     F = fsparse(lDofh(2*ceil(ny/3)+1)+1,1,-1,[nDof,1]);
%     free    = setdiff(1:nDof,fixed);                                      
%     act     = setdiff((1:nEl)',union(pasS,pasV)); 

%%  MBB beam example
%
%     fixed   = [1:2:2*(ny+1), nDof];
%     [pasS,pasV] = deal([],[]);
%     F = fsparse(2,1,-1,[nDof,1]);
%     free    = setdiff(1:nDof,fixed);                                      
%     act     = setdiff((1:nEl)',union(pasS,pasV)); 

%% Specification of circular non-design void
%
%     [cx, cy, cr] = deal(nx/2, ny/2, ny/3);
%     [dy,dx] = meshgrid(1:nx, 1:ny);
%     f = sqrt((dx-cy).^2+(dy-cx).^2) < cr;
