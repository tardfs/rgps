function L1 = nlib_L1_baseline(fout, L1, measurementTime, sat_list, Eph, prMesA, prMesB, cpMesA, cpMesB)
% L1 - is L1 state structure

easyLib = getFullPath('..\\easy') ;

fprintf(fout,'<L1_BASELINE><time: %d>\n', measurementTime ) ;

M = length(sat_list) ; % number of sattelites available 
if M<5
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
    
    ra = 10.24 ;            % code variances
    rb = 
    
    
    % init L1_state
    % sat_list - is the array of GPS sattelite IDs
    L1.sat_list = sat_list ;
    
    % choose base sattelite index in the list
    L1.base_sv = sat_list(1) ;
    fprintf(fout,'\t<BASE SATTELITE ID:><%d>\n', L1_state.base_sv ) ;
    
    % get rough A & B positions
    addpath(easyLib) ;
    [A_pos, ~, ~] = recpo_ls(prMesA(:),(1:M)', measurementTime, Eph ) ;
    [B_pos, ~, ~] = recpo_ls(prMesB(:),(1:M)', measurementTime, Eph ) ;
    rmpath(easyLib) ;
    
    % initialize state vector
    % x(1)..x(3) are relative ECEF coordinates
    % x(1) - ECEF delta X
    % x(2) - ECEF delta Y
    % x(3) - ECEF delta Z
    x_len = 7 ;
    x = zeros(x_len,1) ;
    x(1:3) = A_pos - B_pos ;

    % x(4)..x(7) are double-difference carrier-phase ambiguity terms
    % $\Delta\nabla\phi^{jk}_{AB}$
    % x(4) - $\Delta\nabla\phi^{12}_{AB}$
    % x(5) - $\Delta\nabla\phi^{13}_{AB}$
    % x(6) - $\Delta\nabla\phi^{14}_{AB}$
    % x(7) - $\Delta\nabla\phi^{15}_{AB}$
    k = 4 ;
    dphi_AB_1 = cpMesA(1) - cpMesB(1) ;
    for n=2:5
        x(k) = dphi_AB_1 - (cpMesA(n)-cpMesB(n)) ;
        k = k + 1 ;
    end
    
    % State vector covariance matrix
    P = diag([sigma2_ecef sigma2_ecef sigma2_ecef sigma2_dN sigma2_dN sigma2_dN sigma2_dN]) ;
    
    % Discrete noise matrix
    Qd = diag([qp qp qp qN qN qN qN]*dt) ;
    
    % Measurement covariance matrix
    
    
end
