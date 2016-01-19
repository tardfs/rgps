function L1_state = nlib_L1_baseline(L1_state, measurementTime, sat_list, Eph, obsA, obsB)

if isempty(L1_state)
    % init L1_state
    % sat_list - is the array of GPS sattelite IDs
    L1_state.sat_list = sat_list ;
    
    % get rough A & B positions
    [A_pos, ~, ~] = recpo_ls(obsA,(1:length(obsA)), measurementTime, Eph ) ;
    [B_pos, ~, ~] = recpo_ls(obsB,(1:length(obsB)), measurementTime, Eph ) ;
    
    % initialize state vector
    % x(1) - ECEF delta X
    % x(2) - ECEF delta Y
    % x(3) - ECEF delta Z
    x = zeros(6,1) ;
    x(1:3) = A_pos - B_pos ;    
end
