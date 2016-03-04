function L1 = nlib_L1_baseline(fout, L1, settings, measurementTime, sat_list, Eph, prMesA, prMesB, cpMesA, cpMesB)

% some physical constants
v_light = 299792458 ;	     % vacuum speed of light m/s
f1 = 154*10.23e6 ;		     % L1 frequency, Hz
lambda1 = v_light/f1 ;	     % wavelength on L1:  .19029367  m

% L1 - is L1 state structure

easyLib = getFullPath('..\\easy') ;

fprintf(fout,'<L1_BASELINE><time: %d>\n', measurementTime ) ;

numSat = length(sat_list) ; % number of sattelites available 
if numSat<settings.minSatNum
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
    qp = 30 ;
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
    fprintf(fout,'\t<BASE SATTELITE ID:><%d>\n', L1.base_sv ) ;
    
    % get rough A & B positions
    addpath(easyLib) ;
    [A_pos, ~, ~, A_basic_obs] = recpo_ls(prMesA(:),(1:numSat)', measurementTime, Eph ) ;
    [B_pos, ~, ~, ~] = recpo_ls(prMesB(:),(1:numSat)', measurementTime, Eph ) ;
    rmpath(easyLib) ;
    
    satPos = A_basic_obs(:,1:3) ;
    
    % numSat determines size x_len of state vector, x_len = 3(ECEF) +
    % (numSat-1) double differences
    x_len = 3+ (numSat-1) ;
    
    % initialize state vector
    % x(1)..x(3) are relative ECEF coordinates
    % x(1) - ECEF delta X
    % x(2) - ECEF delta Y
    % x(3) - ECEF delta Z
    x = zeros(x_len,1) ;
    x(1:3) = A_pos(1:3) - B_pos(1:3) ;

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
    
    % Fundamental matrix
    F = eye(x_len) ;
    
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
    
    % compute Esv matrix
    % each row is union vector e corresponded to n-th sattelite
    Esv = satPos - repmat(A_pos(1:3).',numSat,1) ;
    for s=1:numSat
        Esv(s,:) = Esv(s,:)/norm(Esv(s,:),2) ;
    end
    
    % Observation marix
    H = zeros((numSat-1)*2,x_len) ;
    for n=1:numSat-1
        H(n,1:3) = Esv(1,:) - Esv(n+1,:) ;
        H(n+(numSat-1),1:3) = H(n,1:3)/lambda1 ;
        H(n+(numSat-1),n+3) = 1 ;
    end    
    
    L1.x = x ;
    L1.F = F ;
    L1.P = P ;
    L1.Qd = Qd ;
    L1.R = R ;
    L1.H = H ;
    
else
    % checking for sat list
    rbldFlag = 0 ;
    if length(L1.sat_list)~=length(sat_list)
        rbldFlag = 1 ;
    end
    maplist = zeros(size(L1.sat_list)) ;    
    for n=1:length(L1.sat_list)
        k = find(sat_list==L1.sat_list(n),1) ;
        if isempty(k)
            disp(L1.sat_list) ;
            disp(sat_list) ;
            error('can`t resolve sattelite list\n') ;            
        end
        maplist(n) = k ;
        if maplist(n)~=n
            rbldFlag = 1 ;
        end
    end
    if rbldFlag
        % rebuild sattelite list and corresponded vectors
        sat_list = L1.sat_list ;
        numSat = length(sat_list) ; % number of sattelites available 
        Eph = Eph(:,maplist) ;
        Eph(1,:) = 1:numSat ;
        prMesA = prMesA(maplist) ;
        prMesB = prMesB(maplist) ;
        cpMesA = cpMesA(maplist) ;
        cpMesB = cpMesB(maplist) ;
    end
    
    addpath(easyLib) ;
    [A_pos, ~, ~, A_basic_obs] = recpo_ls(prMesA(:),(1:numSat)', measurementTime, Eph ) ;
    rmpath(easyLib) ;
    satPos = A_basic_obs(:,1:3) ;
    % compute Esv matrix
    % each row is union vector e corresponded to n-th sattelite
    Esv = satPos - repmat(A_pos(1:3).',numSat,1) ;
    for s=1:numSat
        Esv(s,:) = Esv(s,:)/norm(Esv(s,:),2) ;
    end
    
end

x = L1.x ;
F = L1.F ;
P = L1.P ;
Qd = L1.Qd ;
R = L1.R ;
H = L1.H ;

% prediction stage
x = F*x ;
P = F*P*F' + Qd ;

% ambiguity resolution
addpath(easyLib) ;
%[a,~,~,~] = lambda(x(4:end), P(4:end,4:end)) ;
%x(1:3) = x(1:3) - P(1:3,4:end)*inv(P(4:end,4:end))*(x(4:end)-a(:,1)) ;
%x(4:end) = a(:,1) ;
rmpath(easyLib) ;

% get measurement vector z
zrho = zeros(numSat-1,1) ;
zphi = zeros(numSat-1,1) ;
drho_AB_1 = prMesA(1) - prMesB(1) ;
dphi_AB_1 = cpMesA(1) - cpMesB(1) ;
for n=2:numSat
    zrho(n-1) = drho_AB_1 - (prMesA(n)-prMesB(n)) ;
    zphi(n-1) = dphi_AB_1 - (cpMesA(n)-cpMesB(n)) ;
end
z = [zrho;zphi] ;
L1.z = z ;

% correction stage
K = P*H'*inv(H*P*H'+R) ;
x = x + K*(z - H*x) ;
P = P - K*H*P ;

% save context
L1.x = x ;
L1.F = F ;
L1.P = P ;
L1.Qd = Qd ;
L1.R = R ;
L1.H = H ;
