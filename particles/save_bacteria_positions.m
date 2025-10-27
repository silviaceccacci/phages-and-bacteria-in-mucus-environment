function bact_positions_over_time = save_bacteria_positions(bacteria, step_idx, bact_positions_over_time)
% saveBacteriaPositions
% ---------------------------------------------------------------
% Stores the positions of all bacteria at a given time step
% into a preallocated matrix bact_positions_over_time.
%
% Usage:
%   bact_positions_over_time = saveBacteriaPositions(bacteria, step_idx, bact_positions_over_time)
%
% Inputs:
%   bacteria                : struct array with field .position (2×1 vector)
%   step_idx                : current time step index (integer)
%   bact_positions_over_time: preallocated matrix of size [num_steps, 2*num_bacteria]
%
% Output:
%   bact_positions_over_time: updated matrix with current positions stored
%
% Each row of bact_positions_over_time stores:
%   [x1, y1, x2, y2, ..., xN, yN]
% where bacterium i’s position corresponds to columns (2*i - 1 : 2*i).
%
% ---------------------------------------------------------------

    % Extract all bacteria positions as an N×2 matrix
    positionsB = reshape([bacteria.position], 2, []).';  % each row = [x, y]

    % Flatten and assign to the current time step
    bact_positions_over_time(step_idx, :) = reshape(positionsB.', 1, []);
end
