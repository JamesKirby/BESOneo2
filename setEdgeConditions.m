% _________________________________________________________________________
%
% Title:        Set Edge Conditions - BESOneo2
%                                                                      
% Version:      0.1 (14/06/2022)                                        
%                                                                      
% Author:       James Kirby, PhD Candidate, RMIT University             
% _________________________________________________________________________
%
% Description:  Set edge conditions in the following order: periodic and
%               symmetric -> applied loads replicate -> applied
%               displacement constraints replicate
%
% Inputs:       ny          number of elements along x
%               nx          number of elements along y
%               fixed       fixed dofs                                      indices, dof indexed
%               F           force vector                                    sparse, dof indexed
%               sym         symmetric boundary segments                     edge segments, node indexed
%               cyc         periodic boundary segments                      edge segments, node indexed
%
% Outputs       EC          edge conditions                                 (ny+2) by (nx+2) matrix
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function EC = setEdgeConditions(nx, ny, fixed, F, sym, cyc)

    EC          = zeros(ny+2, nx+2);
    nodeNrs     = reshape(1:(1+nx)*(1+ny), 1+ny, 1+nx);
    repL        = unique(ceil(find(F)./2));
    repD        = unique(ceil(fixed(:)./2));
    upperEdge   = nodeNrs(1,:);
    lowerEdge   = nodeNrs(end,:);
    leftEdge    = nodeNrs(:,1);
    rightEdge   = nodeNrs(:,end);

    %% Set Periodic boundaries
    if ~isempty(cyc)
        % Check corner conditions first
        if ismember(nodeNrs(1,1), cyc) && ~ismember(nodeNrs(2,1), cyc) && ~ismember(nodeNrs(1,2), cyc)
            EC(1,1) = 1; EC(2,1) = 1; EC(1,2) = 1;
        elseif ismember(nodeNrs(1,1), cyc) && ismember(nodeNrs(2,1), cyc) && ismember(nodeNrs(1,2), cyc)
            EC(1,1) = 1;
        end
        if ismember(nodeNrs(end,1), cyc) && ~ismember(nodeNrs(end-1,1), cyc) && ~ismember(nodeNrs(end,2), cyc)
            EC(end,1) = 1; EC(end-1,1) = 1; EC(end,2) = 1;
        elseif ismember(nodeNrs(end,1), cyc) && ismember(nodeNrs(end-1,1), cyc) && ismember(nodeNrs(end,2), cyc)
            EC(end,1) = 1;
        end
        if ismember(nodeNrs(1,end), cyc) && ~ismember(nodeNrs(1, end-1), cyc) && ~ismember(nodeNrs(2,end), cyc)
            EC(1,end) = 1; EC(2,end) = 1; EC(1,end-1) = 1;
        elseif ismember(nodeNrs(1,end), cyc) && ismember(nodeNrs(1, end-1), cyc) && ismember(nodeNrs(2,end), cyc)
            EC(1,end) = 1;
        end
        if ismember(nodeNrs(end,end), cyc) && ~ismember(nodeNrs(end-1,end), cyc) && ~ismember(nodeNrs(end,end-1), cyc)
            EC(end,end) = 1; EC(end-1,end) = 1; EC(end,end-1) = 1;
        elseif ismember(nodeNrs(end,end), cyc) && ismember(nodeNrs(end-1,end), cyc) && ismember(nodeNrs(end,end-1), cyc)
            EC(end,end) = 1;
        end
        
        % Upper boundary
        UE      = ismember(upperEdge, cyc);
        ECsub   = EC(1, 2:end-1);
        UEECid  = ECsub > 0;
        UECC    = bwconncomp(UE(2:end-1));
        for ii = 1:UECC.NumObjects
            if numel(UECC.PixelIdxList{ii}) == 1
                UE2     = zeros(1,nx);
                UE2(UECC.PixelIdxList{ii}:UECC.PixelIdxList{ii}+1) = 1;
                UE2(UEECid) = ECsub(UEECid);
                EC(1, 2:end-1) = UE2;
                ECsub   = EC(1, 2:end-1);
                UEECid  = ECsub > 0;
            end
        end
        UE      = double(UE(1:end-1) & UE(2:end))*1;
        UE(UEECid) = ECsub(UEECid);
        EC(1, 2:end-1) = UE;
        
        % Lower boundary
        LE      = ismember(lowerEdge, cyc);
        ECsub   = EC(end, 2:end-1);
        LEECid  = ECsub > 0;
        LECC    = bwconncomp(LE(2:end-1));
        for ii = 1:LECC.NumObjects
            if numel(LECC.PixelIdxList{ii}) == 1
                LE2     = zeros(1,nx);
                LE2(LECC.PixelIdxList{ii}:LECC.PixelIdxList{ii}+1) = 1;
                LE2(LEECid) = ECsub(LEECid);
                EC(end, 2:end-1) = LE2;
                ECsub   = EC(end, 2:end-1);
                LEECid  = ECsub > 0;
            end
        end
        LE      = double(LE(1:end-1) & LE(2:end))*1;
        LE(LEECid) = ECsub(LEECid);
        EC(end, 2:end-1) = LE;
        
        % LHS boundary
        LHE     = ismember(leftEdge, cyc);
        ECsub   = EC(2:end-1, 1);
        LHEECid = ECsub > 0;
        LHECC   = bwconncomp(LHE(2:end-1));
        for ii = 1:LHECC.NumObjects
            if numel(LHECC.PixelIdxList{ii}) == 1
                LHE2    = zeros(ny,1);
                LHE2(LHECC.PixelIdxList{ii}:LHECC.PixelIdxList{ii}+1) = 1;
                LHE2(LHEECid) = ECsub(LHEECid);
                EC(2:end-1, 1) = LHE2;
                ECsub   = EC(2:end-1, 1);
                LHEECid = ECsub > 0;
            end
        end
        LHE     = double(LHE(1:end-1) & LHE(2:end))*1;
        LHE(LHEECid) = ECsub(LHEECid);
        EC(2:end-1, 1) = LHE;
        
        % RHS boundary
        RHE     = ismember(rightEdge, cyc);
        ECsub   = EC(2:end-1, end);
        RHEECid = ECsub > 0;
        RHECC   = bwconncomp(RHE(2:end-1));
        for ii = 1:RHECC.NumObjects
            if numel(RHECC.PixelIdxList{ii}) == 1
                RHE2    = zeros(ny,1);
                RHE2(RHECC.PixelIdxList{ii}:RHECC.PixelIdxList{ii}+1) = 1;
                RHE2(RHEECid) = ECsub(RHEECid);
                EC(2:end-1, end) = RHE2;
                ECsub   = EC(2:end-1, end);
                RHEECid = ECsub > 0;
            end
        end
        RHE     = double(RHE(1:end-1) & RHE(2:end))*1;
        RHE(RHEECid) = ECsub(RHEECid);
        EC(2:end-1, end) = RHE;
    end

    %% Set Symmetric boundaries
    if ~isempty(sym)
        % Check corner conditions first
        if ismember(nodeNrs(1,1), sym) && ~ismember(nodeNrs(2,1), sym) && ~ismember(nodeNrs(1,2), sym)
            if EC(1,1) == 0; EC(1,1) = 2; end
            if EC(2,1) == 0; EC(2,1) = 2; end
            if EC(1,2) == 0; EC(1,2) = 2; end
        elseif ismember(nodeNrs(1,1), sym) && ismember(nodeNrs(2,1), sym) && ismember(nodeNrs(1,2), sym)
            if EC(1,1) == 0; EC(1,1) = 2; end
        end
        if ismember(nodeNrs(end,1), sym) && ~ismember(nodeNrs(end-1,1), sym) && ~ismember(nodeNrs(end,2), sym)
            if EC(end,1) == 0; EC(end,1) = 2; end
            if EC(end-1,1) == 0; EC(end-1,1) = 2; end
            if EC(end,2) == 0; EC(end,2) = 2; end
        elseif ismember(nodeNrs(end,1), sym) && ismember(nodeNrs(end-1,1), sym) && ismember(nodeNrs(end,2), sym)
            if EC(end,1) == 0; EC(end,1) = 2; end
        end
        if ismember(nodeNrs(1,end), sym) && ~ismember(nodeNrs(1, end-1), sym) && ~ismember(nodeNrs(2,end), sym)
            if EC(1,end) == 0; EC(1,end) = 2; end
            if EC(2,end) == 0; EC(2,end) = 2; end
            if EC(1,end-1) == 0; EC(1,end-1) = 2; end
        elseif ismember(nodeNrs(1,end), sym) && ismember(nodeNrs(1, end-1), sym) && ismember(nodeNrs(2,end), sym)
            if EC(1,end) == 0; EC(1,end) = 2; end
        end
        if ismember(nodeNrs(end,end), sym) && ~ismember(nodeNrs(end-1,end), sym) && ~ismember(nodeNrs(end,end-1), sym)
            if EC(end,end) == 0; EC(end,end) = 2; end
            if EC(end-1,end) == 0; EC(end-1,end) = 2; end
            if EC(end,end-1) == 0; EC(end,end-1) = 2; end
        elseif ismember(nodeNrs(end,end), sym) && ismember(nodeNrs(end-1,end), sym) && ismember(nodeNrs(end,end-1), sym)
            if EC(end,end) == 0; EC(end,end) = 2; end
        end
        
        % Upper boundary
        UE      = ismember(upperEdge, sym);
        ECsub   = EC(1, 2:end-1);
        UEECid  = ECsub > 0;
        UECC    = bwconncomp(UE(2:end-1));
        for ii = 1:UECC.NumObjects
            if numel(UECC.PixelIdxList{ii}) == 1
                UE2     = zeros(1,nx);
                UE2(UECC.PixelIdxList{ii}:UECC.PixelIdxList{ii}+1) = 2;
                UE2(UEECid) = ECsub(UEECid);
                EC(1, 2:end-1) = UE2;
                ECsub   = EC(1, 2:end-1);
                UEECid  = ECsub > 0;
            end
        end
        UE      = double(UE(1:end-1) & UE(2:end))*2;
        UE(UEECid) = ECsub(UEECid);
        EC(1, 2:end-1) = UE;
        
        % Lower boundary
        LE      = ismember(lowerEdge, sym);
        ECsub   = EC(end, 2:end-1);
        LEECid  = ECsub > 0;
        LECC    = bwconncomp(LE(2:end-1));
        for ii = 1:LECC.NumObjects
            if numel(LECC.PixelIdxList{ii}) == 1
                LE2     = zeros(1,nx);
                LE2(LECC.PixelIdxList{ii}:LECC.PixelIdxList{ii}+1) = 2;
                LE2(LEECid) = ECsub(LEECid);
                EC(end, 2:end-1) = LE2;
                ECsub   = EC(end, 2:end-1);
                LEECid  = ECsub > 0;
            end
        end
        LE      = double(LE(1:end-1) & LE(2:end))*2;
        LE(LEECid) = ECsub(LEECid);
        EC(end, 2:end-1) = LE;
        
        % LHS boundary
        LHE     = ismember(leftEdge, sym);
        ECsub   = EC(2:end-1, 1);
        LHEECid = ECsub > 0;
        LHECC   = bwconncomp(LHE(2:end-1));
        for ii = 1:LHECC.NumObjects
            if numel(LHECC.PixelIdxList{ii}) == 1
                LHE2    = zeros(ny,1);
                LHE2(LHECC.PixelIdxList{ii}:LHECC.PixelIdxList{ii}+1) = 2;
                LHE2(LHEECid) = ECsub(LHEECid);
                EC(2:end-1, 1) = LHE2;
                ECsub   = EC(2:end-1, 1);
                LHEECid = ECsub > 0;
            end
        end
        LHE     = double(LHE(1:end-1) & LHE(2:end))*2;
        LHE(LHEECid) = ECsub(LHEECid);
        EC(2:end-1, 1) = LHE;
        
        % RHS boundary
        RHE     = ismember(rightEdge, sym);
        ECsub   = EC(2:end-1, end);
        RHEECid = ECsub > 0;
        RHECC   = bwconncomp(RHE(2:end-1));
        for ii = 1:RHECC.NumObjects
            if numel(RHECC.PixelIdxList{ii}) == 1
                RHE2    = zeros(ny,1);
                RHE2(RHECC.PixelIdxList{ii}:RHECC.PixelIdxList{ii}+1) = 2;
                RHE2(RHEECid) = ECsub(RHEECid);
                EC(2:end-1, end) = RHE2;
                ECsub   = EC(2:end-1, end);
                RHEECid = ECsub > 0;
            end
        end
        RHE     = double(RHE(1:end-1) & RHE(2:end))*2;
        RHE(RHEECid) = ECsub(RHEECid);
        EC(2:end-1, end) = RHE;
    end

    %% Set Replicate boundaries
        % Check corner conditions for loads
        if ismember(nodeNrs(1,1), repL) && ~ismember(nodeNrs(2,1), repL) && ~ismember(nodeNrs(1,2), repL)
            if EC(1,1) == 0; EC(1,1) = 3; end
            if EC(2,1) == 0; EC(2,1) = 3; end
            if EC(1,2) == 0; EC(1,2) = 3; end
        elseif ismember(nodeNrs(1,1), repL) && ismember(nodeNrs(2,1), repL) && ismember(nodeNrs(1,2), repL)
            if EC(1,1) == 0; EC(1,1) = 3; end
        end
        if ismember(nodeNrs(end,1), repL) && ~ismember(nodeNrs(end-1,1), repL) && ~ismember(nodeNrs(end,2), repL)
            if EC(end,1) == 0; EC(end,1) = 3; end
            if EC(end-1,1) == 0; EC(end-1,1) = 3; end
            if EC(end,2) == 0; EC(end,2) = 3; end
        elseif ismember(nodeNrs(end,1), repL) && ismember(nodeNrs(end-1,1), repL) && ismember(nodeNrs(end,2), repL)
            if EC(end,1) == 0; EC(end,1) = 3; end
        end
        if ismember(nodeNrs(1,end), repL) && ~ismember(nodeNrs(1, end-1), repL) && ~ismember(nodeNrs(2,end), repL)
            if EC(1,end) == 0; EC(1,end) = 3; end
            if EC(2,end) == 0; EC(2,end) = 3; end
            if EC(1,end-1) == 0; EC(1,end-1) = 3; end
        elseif ismember(nodeNrs(1,end), repL) && ismember(nodeNrs(1, end-1), repL) && ismember(nodeNrs(2,end), repL)
            if EC(1,end) == 0; EC(1,end) = 3; end
        end
        if ismember(nodeNrs(end,end), repL) && ~ismember(nodeNrs(end-1,end), repL) && ~ismember(nodeNrs(end,end-1), repL)
            if EC(end,end) == 0; EC(end,end) = 3; end
            if EC(end-1,end) == 0; EC(end-1,end) = 3; end
            if EC(end,end-1) == 0; EC(end,end-1) = 3; end
        elseif ismember(nodeNrs(end,end), repL) && ismember(nodeNrs(end-1,end), repL) && ismember(nodeNrs(end,end-1), repL)
            if EC(end,end) == 0; EC(end,end) = 3; end
        end
        % Check corner conditions for displacement boundaries
        if ismember(nodeNrs(1,1), repD) && ~ismember(nodeNrs(2,1), repD) && ~ismember(nodeNrs(1,2), repD)
            if EC(1,1) == 0; EC(1,1) = 3; end
            if EC(2,1) == 0; EC(2,1) = 3; end
            if EC(1,2) == 0; EC(1,2) = 3; end
        elseif ismember(nodeNrs(1,1), repD) && ismember(nodeNrs(2,1), repD) && ismember(nodeNrs(1,2), repD)
            if EC(1,1) == 0; EC(1,1) = 3; end
        end
        if ismember(nodeNrs(end,1), repD) && ~ismember(nodeNrs(end-1,1), repD) && ~ismember(nodeNrs(end,2), repD)
            if EC(end,1) == 0; EC(end,1) = 3; end
            if EC(end-1,1) == 0; EC(end-1,1) = 3; end
            if EC(end,2) == 0; EC(end,2) = 3; end
        elseif ismember(nodeNrs(end,1), repD) && ismember(nodeNrs(end-1,1), repD) && ismember(nodeNrs(end,2), repD)
            if EC(end,1) == 0; EC(end,1) = 3; end
        end
        if ismember(nodeNrs(1,end), repD) && ~ismember(nodeNrs(1, end-1), repD) && ~ismember(nodeNrs(2,end), repD)
            if EC(1,end) == 0; EC(1,end) = 3; end
            if EC(2,end) == 0; EC(2,end) = 3; end
            if EC(1,end-1) == 0; EC(1,end-1) = 3; end
        elseif ismember(nodeNrs(1,end), repD) && ismember(nodeNrs(1, end-1), repD) && ismember(nodeNrs(2,end), repD)
            if EC(1,end) == 0; EC(1,end) = 3; end
        end
        if ismember(nodeNrs(end,end), repD) && ~ismember(nodeNrs(end-1,end), repD) && ~ismember(nodeNrs(end,end-1), repD)
            if EC(end,end) == 0; EC(end,end) = 3; end
            if EC(end-1,end) == 0; EC(end-1,end) = 3; end
            if EC(end,end-1) == 0; EC(end,end-1) = 3; end
        elseif ismember(nodeNrs(end,end), repD) && ismember(nodeNrs(end-1,end), repD) && ismember(nodeNrs(end,end-1), repD)
            if EC(end,end) == 0; EC(end,end) = 3; end
        end
        
        % Check remaining conditions for loads
        % Upper boundary
        UE      = ismember(upperEdge, repL);
        ECsub   = EC(1, 2:end-1);
        UEECid  = ECsub > 0;
        UECC    = bwconncomp(UE(2:end-1));
        for ii = 1:UECC.NumObjects
            if numel(UECC.PixelIdxList{ii}) == 1
                UE2     = zeros(1,nx);
                UE2(UECC.PixelIdxList{ii}:UECC.PixelIdxList{ii}+1) = 3;
                UE2(UEECid) = ECsub(UEECid);
                EC(1, 2:end-1) = UE2;
                ECsub   = EC(1, 2:end-1);
                UEECid  = ECsub > 0;
            else
                LL      = numel(UECC.PixelIdxList{ii});
                LE      = LL+1;
                step    = [1:4:LL/2, 2:4:LL/2]; 
                apply   = zeros(1, LE);
                apply(step) = 1;
                apply   = apply + fliplr(apply);
                apply(apply == 1) = 3;
                apply(apply == 0) = 4;
                EC(1, 2:end-1) = apply;
                ECsub   = EC(1, 2:end-1);
                UEECid  = ECsub > 0;
            end
        end
        UE      = double(UE(1:end-1) & UE(2:end))*3;
        UE(UEECid) = ECsub(UEECid);
        EC(1, 2:end-1) = UE;
        
        % Lower boundary
        LE      = ismember(lowerEdge, repL);
        ECsub   = EC(end, 2:end-1);
        LEECid  = ECsub > 0;
        LECC    = bwconncomp(LE(2:end-1));
        for ii = 1:LECC.NumObjects
            if numel(LECC.PixelIdxList{ii}) == 1
                LE2     = zeros(1,nx);
                LE2(LECC.PixelIdxList{ii}:LECC.PixelIdxList{ii}+1) = 3;
                LE2(LEECid) = ECsub(LEECid);
                EC(end, 2:end-1) = LE2;
                ECsub   = EC(end, 2:end-1);
                LEECid  = ECsub > 0;
            else
                LL      = numel(LECC.PixelIdxList{ii});
                LE      = LL+1;
                step    = [1:4:LL/2, 2:4:LL/2]; 
                apply   = zeros(1, LE);
                apply(step) = 1;
                apply   = apply + fliplr(apply);
                apply(apply == 1) = 3;
                apply(apply == 0) = 4;
                EC(end, 2:end-1) = apply;
                ECsub   = EC(end, 2:end-1);
                LEECid  = ECsub > 0;    
            end
        end
        LE      = double(LE(1:end-1) & LE(2:end))*3;
        LE(LEECid) = ECsub(LEECid);
        EC(end, 2:end-1) = LE;
        
        % LHS boundary
        LHE     = ismember(leftEdge, repL);
        ECsub   = EC(2:end-1, 1);
        LHEECid = ECsub > 0;
        LHECC   = bwconncomp(LHE(2:end-1));
        for ii = 1:LHECC.NumObjects
            if numel(LHECC.PixelIdxList{ii}) == 1
                LHE2    = zeros(ny,1);
                LHE2(LHECC.PixelIdxList{ii}:LHECC.PixelIdxList{ii}+1) = 3;
                LHE2(LHEECid) = ECsub(LHEECid);
                EC(2:end-1, 1) = LHE2;
                ECsub   = EC(2:end-1, 1);
                LHEECid = ECsub > 0;
            else
                LL      = numel(LHECC.PixelIdxList{ii});
                LE      = LL+1;
                step    = [1:4:LL/2, 2:4:LL/2]'; 
                apply   = zeros(LE, 1);
                apply(step) = 1;
                apply   = apply + flipud(apply);
                apply(apply == 1) = 3;
                apply(apply == 0) = 4;
                EC(2:end-1, 1) = apply;
                ECsub   = EC(2:end-1, 1);
                LHEECid = ECsub > 0;
            end
        end
        LHE     = double(LHE(1:end-1) & LHE(2:end))*3;
        LHE(LHEECid) = ECsub(LHEECid);
        EC(2:end-1, 1) = LHE;
        
        % RHS boundary
        RHE     = ismember(rightEdge, repL);
        ECsub   = EC(2:end-1, end);
        RHEECid = ECsub > 0;
        RHECC   = bwconncomp(RHE(2:end-1));
        for ii = 1:RHECC.NumObjects
            if numel(RHECC.PixelIdxList{ii}) == 1
                RHE2    = zeros(ny,1);
                RHE2(RHECC.PixelIdxList{ii}:RHECC.PixelIdxList{ii}+1) = 3;
                RHE2(RHEECid) = ECsub(RHEECid);
                EC(2:end-1, end) = RHE2;
                ECsub   = EC(2:end-1, end);
                RHEECid = ECsub > 0;
            else
                LL      = numel(RHECC.PixelIdxList{ii});
                LE      = LL+1;
                step    = [1:4:LL/2, 2:4:LL/2]'; 
                apply   = zeros(LE, 1);
                apply(step) = 1;
                apply   = apply + flipud(apply);
                apply(apply == 1) = 3;
                apply(apply == 0) = 4;
                EC(2:end-1, end) = apply;
                ECsub   = EC(2:end-1, end);
                RHEECid = ECsub > 0;    
            end
        end
        RHE     = double(RHE(1:end-1) & RHE(2:end))*3;
        RHE(RHEECid) = ECsub(RHEECid);
        EC(2:end-1, end) = RHE;

        % Check remaining conditions for displacement boundaries
        % Upper boundary
        UE      = ismember(upperEdge, repD);
        ECsub   = EC(1, 2:end-1);
        UEECid  = ECsub > 0;
        UECC    = bwconncomp(UE(2:end-1));
        for ii = 1:UECC.NumObjects
            if numel(UECC.PixelIdxList{ii}) == 1
                UE2     = zeros(1,nx);
                UE2(UECC.PixelIdxList{ii}:UECC.PixelIdxList{ii}+1) = 3;
                UE2(UEECid) = ECsub(UEECid);
                EC(1, 2:end-1) = UE2;
                ECsub   = EC(1, 2:end-1);
                UEECid  = ECsub > 0;
            end
        end
        UE      = double(UE(1:end-1) & UE(2:end))*3;
        UE(UEECid) = ECsub(UEECid);
        EC(1, 2:end-1) = UE;
        
        % Lower boundary
        LE      = ismember(lowerEdge, repD);
        ECsub   = EC(end, 2:end-1);
        LEECid  = ECsub > 0;
        LECC    = bwconncomp(LE(2:end-1));
        for ii = 1:LECC.NumObjects
            if numel(LECC.PixelIdxList{ii}) == 1
                LE2     = zeros(1,nx);
                LE2(LECC.PixelIdxList{ii}:LECC.PixelIdxList{ii}+1) = 3;
                LE2(LEECid) = ECsub(LEECid);
                EC(end, 2:end-1) = LE2;
                ECsub   = EC(end, 2:end-1);
                LEECid  = ECsub > 0;
            end
        end
        LE      = double(LE(1:end-1) & LE(2:end))*3;
        LE(LEECid) = ECsub(LEECid);
        EC(end, 2:end-1) = LE;
        
        % LHS boundary
        LHE     = ismember(leftEdge, repD);
        ECsub   = EC(2:end-1, 1);
        LHEECid = ECsub > 0;
        LHECC   = bwconncomp(LHE(2:end-1));
        for ii = 1:LHECC.NumObjects
            if numel(LHECC.PixelIdxList{ii}) == 1
                LHE2    = zeros(ny,1);
                LHE2(LHECC.PixelIdxList{ii}:LHECC.PixelIdxList{ii}+1) = 3;
                LHE2(LHEECid) = ECsub(LHEECid);
                EC(2:end-1, 1) = LHE2;
                ECsub   = EC(2:end-1, 1);
                LHEECid = ECsub > 0;
            end
        end
        LHE     = double(LHE(1:end-1) & LHE(2:end))*3;
        LHE(LHEECid) = ECsub(LHEECid);
        EC(2:end-1, 1) = LHE;
        
        % RHS boundary
        RHE     = ismember(rightEdge, repD);
        ECsub   = EC(2:end-1, end);
        RHEECid = ECsub > 0;
        RHECC   = bwconncomp(RHE(2:end-1));
        for ii = 1:RHECC.NumObjects
            if numel(RHECC.PixelIdxList{ii}) == 1
                RHE2    = zeros(ny,1);
                RHE2(RHECC.PixelIdxList{ii}:RHECC.PixelIdxList{ii}+1) = 3;
                RHE2(RHEECid) = ECsub(RHEECid);
                EC(2:end-1, end) = RHE2;
                ECsub   = EC(2:end-1, end);
                RHEECid = ECsub > 0;
            end
        end
        RHE     = double(RHE(1:end-1) & RHE(2:end))*3;
        RHE(RHEECid) = ECsub(RHEECid);
        EC(2:end-1, end) = RHE;

        EC(EC==4) = 0;

end % END function
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\