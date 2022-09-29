% _________________________________________________________________________
%
% Title:        BESOneo2 - Efficient Minimum Compliance BESO (2D)
%                                                                      
% Version:      2.0 (14/06/2022)                                        
%                                                                      
% Author:       James Kirby, PhD Candidate, RMIT University             
% _________________________________________________________________________
%
% Description:  Minimum compliance, volume constrained topology 
%               optimization using Bi-directional Evolutionary Topology 
%               Optimization (BESO).
%               FEA method is adapted from "A new generation 99 line Matlab 
%               code for compliance Topology Optimization and its extension 
%               to 3D" by Ferrari, F., & Sigmund, O. (2020)
%
% Notes:        This implementation uses the fsparse function as described
%               in Engblom, S., & Lukarski, D. (2016). Fast Matlab 
%               compatible sparse assembly on multicore computers. Parallel 
%               Computing, 56, 1–17. 
%               https://doi.org/10.1016/j.parco.2016.04.001. The fsparse 
%               code is available from the first authors GitHub page at: 
%               https://github.com/stefanengblom/stenglib. 
%
% Inputs:       nx          Number of elements along x
%               ny          Number of elements along y
%               volfrac     Volume constraint (0 <= volfrac <= 1)
%               er          Evolutionary rate (typically 0.02)
%               rmin        Filter radius
%               definition  Inbuilt examples ["cantilever", "frame",
%                           "Lbracket", "MBB", "cantilever_ndr"]
%
% Outputs       xfinal      Optimised design
%               obj         Mean Compliance for optimised design
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function [xfinal, obj] = BESOneo2(nx, ny, volfrac, er, rmin, definition)
% 0) Preliminaries
    if nargin < 6
        definition = "custom";                                              % remove if defining own structural definition in section D
    elseif nargin < 5
        error("Not enough input arguments");
    end
    list    = dir('**/*.*');                                                % check if fsparse function exists in path
    if any(contains({list.name}, 'fsparse')==1)
        fopt = true;
    else
        fopt = false;                                                       % if not, use MATLABs inbuilt sparse function
    end
% A) Optimisation parameters    
    penal   = 3;                                                            % penalty exponent for material interpolation
    xmin    = 1e-3;                                                         % material density of void elements
    maxit   = 200;                                                          % maximum number of iterations allowed
% B) Material paramaters
    E0      = 1;                                                            % Young modulus of solid
    Emin    = 1e-6;                                                         % Young modulus of "void"
    nu      = 0.3;                                                          % Poisson ratio
% C) Discretization features
    nEL     = nx*ny;                                                        % number of elements
    nodeNrs = int32(reshape(1:(1+nx)*(1+ny), 1+ny, 1+nx));                  % indexed grid node numbers (defined as int32)
    cVec    = reshape(2*nodeNrs(1:end-1, 1:end-1)+1, nEL, 1);
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
    Ke0     = reshape(Ke0, 8, 8);
    Ke0     = Ke0 + Ke0' - diag(diag(Ke0));                                 % recover full element K matrix
% D) Loads, supports & passive domains
   [F, free, fixed, act, pasV, sym, cyc, symX, symY] = ...                  % inbuilt structural definitions
                           predefinedStructure(nx, ny, definition, fopt);   
% E) Prepare filter
    EC      = setEdgeConditions(nx, ny, fixed, F, sym, cyc);                % edge conditions defined based on definitions in [fixed, F, sym, cyc]
    Rmin    = ceil(rmin);
    sigma   = rmin/3;
    [dx,dy] = meshgrid(-Rmin:Rmin, -Rmin:Rmin);
    h       = exp(-(dx.^2+dy.^2)/(2*sigma^2));                              % construct Gaussian kernel
    h       = h./sum(h, 'all');                                             % normalise the kernel
% F) Allocate & initialize other parameters
    [x,dsK] = deal(ones(nEL, 1), zeros(nEL, 1));                            % initialize vectors
    x(pasV) = xmin;                                                         % set x = xmin on pasV set
    [ch, loop, vol, U] = deal(1, 0, sum(x,'all')/nEL, zeros(nDof,1));
    while ch > 1e-5 && loop<maxit
        loop    = loop + 1;                                                 % update loop counter
        vol     = max(vol*(1-er), volfrac);                                 % target volume for current loop iteration
        if loop > 1; olddc = dc; end                                        % save previous iterations sensitivities
% 1) Setup & solve equilibrium equations
        sK      = Emin + x.^penal*(E0-Emin);                                % stiffness interpolation
        dsK     = penal.*x.^(penal-1);                                      % derivative of stiffness interpolation
        sK      = reshape(Ke(:)*sK', length(Ke)*nEL, 1);
        if fopt
            K = fsparse(Iar(:, 1), Iar(:, 2), sK, [nDof, nDof]);            % assemble stiffness matrix
        else
            K = sparse(Iar(:, 1), Iar(:, 2), sK, nDof, nDof); 
        end
        U(free) = decomposition(K(free, free), 'chol', 'lower')\F(free);    % solve equilibrium system
        obj(loop)= 0.5*F'*U;                                                % compute mean compliance
% 2) Compute Sensitivities
        dc      = dsK.*sum((U(cMat)*Ke0).*U(cMat), 2);                      % derivative of compliance
        dc      = applyEdgeConditions(nx, ny, dc, Rmin, EC);                % apply edge conditions to sensitivity field
        dc      = conv2(dc, h, 'same');                                     % filter objective sensitivity
        dc      = dc(Rmin+1:end-Rmin, Rmin+1:end-Rmin);                     % extract original domain
        dc      = applySymmetryConstraints(dc, symX, symY);                 % enforce symmetries if specified in [symX, symY]
        dc      = dc(:);
        if loop>1; dc = mean([dc, olddc], 2); end
% 3) BESO design update
        s = [min(dc), max(dc)];                                             % initial estimate for lambda
        while (s(2)-s(1))/s(2) > 1e-7                                       % BESO design update
            sbar = mean([s(1), s(2)]);
            x(act) = max(xmin, sign(dc(act)-sbar));
            if (sum(x) -nEL*vol) > 0
                s(1) = sbar; 
            else
                s(2) = sbar; 
            end
        end
% 4) Print current results and plot design
        if loop > 10
            ch = abs(sum(obj(loop-9:loop-5))-sum(obj(loop-4:loop)))/...
                                                     sum(obj(loop-4:loop));
        end
        fprintf( 'It.:%7i Obj:%7.4f Vol:%7.3f ch.:%0.2e \n', ...
            loop, obj(loop), vol, ch); 
        colormap(gray); imagesc(reshape(x<1, ny, nx));
        axis equal tight off; drawnow;
    end
    xfinal = reshape(x, ny, nx);
end % END function
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\