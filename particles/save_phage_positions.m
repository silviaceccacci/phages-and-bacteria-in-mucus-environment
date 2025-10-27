function phage_positions_over_time = save_phage_positions(phages, step_idx, phage_positions_over_time)
% savePhagePositions
% ---------------------------------------------------------------
% Stores the positions of all phages at a given time step
% into a preallocated matrix phage_positions_over_time.
%
% Usage:
%   phage_positions_over_time = savePhagePositions(phages, step_idx, phage_positions_over_time)
%
% Inputs:
%   phages                  : struct array with field .position (2×1 vector)
%   step_idx                : current time step index (integer)
%   phage_positions_over_time : preallocated matrix of size [num_steps, 2*num_phages]
%
% Output:
%   phage_positions_over_time : updated matrix with current positions stored
%
% Each row of phage_positions_over_time stores:
%   [x1, y1, x2, y2, ..., xN, yN]
% where phage i’s position corresponds to columns (2*i - 1 : 2*i).
%
% ---------------------------------------------------------------

    num_phages_total = size(phage_positions_over_time,2)/2;
    positionsP = NaN(num_phages_total,2);

    for p = 1:length(phages)
        id = phages(p).id;
        positionsP(id,:) = phages(p).position(:)';
    end

    % Flatten and assign to the current time step
    phage_positions_over_time(step_idx, :) = reshape(positionsP.', 1, []);
end
