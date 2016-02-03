function L1 = nlib_L1_baseline(fout, L1, measurementTime, sat_list, Eph, prMesA, prMesB, cpMesA, cpMesB)
% L1 - is L1 state structure

easyLib = getFullPath('..\\easy') ;

fprintf(fout,'<L1_BASELINE><time: %d>\n', measurementTime ) ;

numSat = length(sat_list) ; % number of sattelites available 
if numSat<5
    fprintf(fout, '<Error!> [nlib_L1_baseline]: not enought sattelites.\n') ;
    return ;
end

% [sat_list, sat_reindex] = sort(sat_list) ;
% Eph = Eph(:,sat_reindex) ;
% obsA = obsA(sat_reindex) ;
% obsB = obsB(sat_reindex) ;

if isempty(L1)

    fprintf(fout, '\t<INIT L1 STATE>\n' ) ;
    
    % filter parameters 
    dt = 2 ;                % time step, sec
    
    % initial values 
    sigma2_ecef =     5 ;   % $\sigma^2_{\Delta X}$ - initial variance of relative position
    sigma2_dN =  50 ;       % $\sigma^2_{\Delta\nabla\phi^{12}_{AB}}$ - initial 
                            % variance of double difference carrier phase
                            % measurements
    qp = 15 ;
    qN = 1.1e-2 ;
    
    ra = 10.24 ;            % code variance
    rb = 0.0994 ;           % phase variance
    rc = 5.12 ;             % code covariance
    rd = 0.0497 ;           % phase covariance
    
    
    % init L1_state
    % sat_list - is the array of GPS sattelite IDs
    L1.sat_list = sat_list ;
    
    % choose base sattelite index in the list
    L1.base_sv = sat_list(1) ;
    fprintf(fout,'\t<BASE SATTELITE ID:><%d>\n', L1_state.base_sv ) ;
    
    % get rough A & B positions
    addpath(easyLib) ;
    [A_pos, ~, ~] = recpo_ls(prMesA(:),(1:numSat)', measurementTime, Eph ) ;
    [B_pos, ~, ~] = recpo_ls(prMesB(:),(1:numSat)', measurementTime, Eph ) ;
    rmpath(easyLib) ;
    
    % numSat determines size x_len of state vector, x_len = 3(ECEF) +
    % (numSat-1) double differences
    x_len = 3+ (numSat-1) ;
    
    % initialize state vector
    % x(1)..x(3) are relative ECEF coordinates
    % x(1) - ECEF delta X
    % x(2) - ECEF delta Y
    % x(3) - ECEF delta Z
    x = zeros(x_len,1) ;
    x(1:3) = A_pos - B_pos ;

    % x(4)..x(end) are double-difference carrier-phase ambiguity terms
    % x(4) - $\Delta\nabla\phi^{12}_{AB}$
    % x(5) - $\Delta\nabla\phi^{13}_{AB}$
    % x(6) - $\Delta\nabla\phi^{14}_{AB}$
    % x(7) - $\Delta\nabla\phi^{15}_{AB}$
    % ...
    k = 4 ; % index of the first double difference at state vector
    dphi_AB_1 = cpMesA(1) - cpMesB(1) ;
    for n=2:numSat
        x(k) = dphi_AB_1 - (cpMesA(n)-cpMesB(n)) ;
        k = k + 1 ;
    end
    
    % State vector covariance matrix P(x_len,x_len)
    P = diag([ones(3,1)*sigma2_ecef; ones(numSat-1,1)*sigma2_dN] ) ;
    
    % Discrete noise covariance matrix
    Qd = diag( [ones(3,1)*qp; ones(numSat-1,1)*qN]*dt ) ;
    
    % Measurement covariance matrix, R = E(z'z), $z=($\Delta\nabla\rho^{12}_{AB}$)$ - measurements vector
    R = zeros((numSat-1)*2, (numSat-1)*2) ; % oservation
    E_matrix = (ones(numSat-1)-diag(ones(numSat-1,1))) ;
    D_matrix = diag(ones(numSat-1,1)) ;
    R(1:(numSat-1),1:(numSat-1)) = E_matrix*rc + D_matrix*ra ;
    R(numSat:end,numSat:end) = E_matrix*rd + D_matrix*rb ;
    
    
end
