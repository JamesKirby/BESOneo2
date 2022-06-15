% _________________________________________________________________________
%
% Title:        Apply Edge Conditions - BESOneo2
%                                                                      
% Version:      0.1 (14/06/2022)                                        
%                                                                      
% Author:       James Kirby, PhD Candidate, RMIT University             
% _________________________________________________________________________
%
% Description:  Apply edge conditions defined by an (ny+2) by (nx+2) matrix
%               where each of the corners (1,1),(1,end),(end,1),(end,end)
%               define edge conditions at the padded corners
%               (1:r,1:r),(1:r,end-r+1:end),(end-r+1:end,1:r),
%               (end-r+1:end,end-r+1:end)
%               and where the remaining border elements define edge
%               conditions for linear boundary segments
%               0: free/zero edge condition
%               1: periodic edge condition
%               2: symmetric edge condition
%               3: replicate edge condition
%
% Inputs:       ny          number of elements along x
%               nx          number of elements along y
%               dc          matrix on which edge conditions are applied
%               r           pad distance
%               EC          edge condition matrix
%
% Outputs       DC          output matrix of size (ny+2r) by (nx+2r)
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function DC = applyEdgeConditions(nx, ny, dc, r, EC)

    dc = reshape(dc, ny, nx);
    DC = padarray(dc, [r r], 0);

    % Apply any corner edge conditions
    % Upper left
    switch EC(1,1)
        case 1 % periodic
            DC(1:r, 1:r) = dc(end-r+1:end, end-r+1:end);
        case 2 % symmetric
            DC(1:r, 1:r) = rot90(fliplr(dc(1:r, 1:r)),-1);
        case 3 % replicate
            DC(1:r, 1:r) = dc(1);
        % otherwise remains zero
    end
    %Lower left
    switch EC(end,1)
        case 1 % periodic
            DC(end-r+1:end, 1:r) = dc(1:r, end-r+1:end);
        case 2 % symmetric
            DC(end-r+1:end, 1:r) = rot90(fliplr(dc(end-r+1, 1:r)),-1);
        case 3 % replicate
            DC(end-r+1:end, 1:r) = dc(end,1);
        % otherwise remains zero
    end
    % Lower right
    switch EC(end,end)
        case 1 % periodic
            DC(end-r+1:end, end-r+1:end) = dc(1:r, 1:r);
        case 2 % symmetric
            DC(end-r+1:end, end-r+1:end) = rot90(fliplr(dc(end-r+1:end, end-r+1:end)),-1);
        case 3 % replicate
            DC(end-r+1:end, end-r+1:end) = dc(end,end);
        % otherwise remains zero
    end
    % Upper right
    switch EC(1,end)
        case 1 % periodic
            DC(1:r, end-r+1:end) = dc(end-r+1:end, 1:r);
        case 2 % symmetric
            DC(1:r, end-r+1:end) = rot90(fliplr(dc(1:r, end-r+1:end)),-1);
        case 3 % replicate
            DC(1:r, end-r+1:end) = dc(1,end);
        % otherwise remains zero
    end

    % Apply edge conditions to linear boundary segments
    upper = EC(1, 2:end-1);
    lower = EC(end, 2:end-1);
    lhs = EC(2:end-1, 1);
    rhs = EC(2:end-1, end);

    for ii = 1:nx
        switch upper(ii)
            case 1 % periodic
                DC(1:r, r+ii) = dc(end-r+1:end, ii);
            case 2 % symmetric
                DC(1:r, r+ii) = flipud(dc(1:r, ii));
            case 3 % replicate
                DC(1:r, r+ii) = dc(1, ii);
            otherwise % otherwise remains zero
                DC(1:r, r+ii) = 0;
        end
        switch lower(ii)
            case 1 % periodic
                DC(end-r+1:end, r+ii) = dc(1:r, ii);
            case 2 % symmetric
                DC(end-r+1:end, r+ii) = flipud(dc(end-r+1:end, ii));
            case 3 % replicate
                DC(end-r+1:end, r+ii) = dc(end, ii);
            otherwise % otherwise remains zero
                DC(end-r+1:end, r+ii) = 0;
        end
    end
    for ii = 1:ny
        switch lhs(ii)
            case 1 % periodic
                DC(r+ii, 1:r) = dc(ii, end-r+1:end);
            case 2 % symmetric
                DC(r+ii, 1:r) = fliplr(dc(ii, 1:r));
            case 3 % replicate
                DC(r+ii, 1:r) = dc(ii, 1);
            otherwise % otherwise remains zero
                DC(r+ii, 1:r) = 0;
        end
        switch rhs(ii)
            case 1 % periodic
                DC(r+ii, end-r+1:end) = dc(ii, 1:r);
            case 2 % symmetric
                DC(r+ii, end-r+1:end) = fliplr(dc(ii, end-r+1:end));
            case 3 % replicate
                DC(r+ii, end-r+1:end) = dc(ii, end);
            otherwise % otherwise remains zero
                DC(r+ii, end-r+1:end) = 0;
        end
    end

end % END function
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\