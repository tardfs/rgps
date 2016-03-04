% Baseline solver using L1 measurments only
function [baseline_data, x_data, P_data, z_data] = nlib_L1_baseline_solver(settings, fout, measurmentsA, measurmentsB)

BL = 0 ; % number of baselines computed

N = length(measurmentsA) ; % number of receiver A measurements
K = length(measurmentsB) ; % number of receiver B measurements

baseline_data = zeros(5,min(N,K)) ;
% 1 - time
% 2:4 - baseline x,y,z
% 5 - reserved (GDOP)

x_data = zeros(7,min(N,K)) ;
P_data = zeros(7,min(N,K)) ;
z_data = zeros(8,min(N,K)) ;


L1_state = [] ;

% compute time step for Kalman filter
% dt_A = zeros(N-1,1) ;
% dt_B = zeros(N-1,1) ;
% for n=2:N
%     dt_A(n-1) = measurmentsA{n}{1}.msrTow - measurmentsA{n-1}{1}.msrTow ;
%     dt_B(n-1) = measurmentsB{n}{1}.msrTow - measurmentsB{n-1}{1}.msrTow ;
% end

% get Tow vectors for A & B receivers
recvA_Tow_v = zeros(N,1) ;
for n=1:N
    recvA_msr = measurmentsA{n} ;
    % get rcvA measurement time
    recvA_Tow_v(n) = recvA_msr{1}.msrTow ;
end
recvB_Tow_v = zeros(K,1) ;
for k=1:K
    recvB_msr = measurmentsB{k} ;
    % get rcvA measurement time
    recvB_Tow_v(k) = recvB_msr{1}.msrTow ;
end

% A & time synchronization
AB_n2k = zeros(N,1) ;
AB_time_err = zeros(N,1) ;
for n=1:N
    [time_dif,k] = min(abs(recvB_Tow_v-recvA_Tow_v(n))) ;
    AB_time_err(n) = time_dif ;
    AB_n2k(n) = k ;
end

if settings.check_for_time_sync
    figure(20) ;
    set(gcf,'NumberTitle','off') ;
    set(gcf,'Name', 'A-B Time mismatch' ) ; 
    hold off, 
    plot((recvA_Tow_v-recvA_Tow_v(1))/60, AB_time_err*settings.v_light,'LineWidth',2) ;
    set(gca,'FontSize',14) ;
    grid on ;
    xlabel('Time, sec') ;
    ylabel('Light travel, m') ;
    title('A-B Time mismatch in terms of correpsonded light travel') ;
    c=input('Proceed with relative coordinates: [Y/N]?','s') ;
    if isempty(c)
        c = 'Y' ;
    end
    if c~='Y' && c~='y'
        return ;
    end
end

for n=1:N
    recvA_msr = measurmentsA{n} ;
    % get rcvA measurement time
    TowA = recvA_msr{1}.msrTow ;
    % find corresponded measurment for recvB
    k = -1 ;
    if AB_time_err(n)<=settings.timeSyncTol
        k = AB_n2k(n) ;
    end
%     for p=1:K
%         recvB_msr = measurmentsB{p} ;
%         TowB = recvB_msr{1}.msrTow ;
% %        if round(TowA)==round(TowB)
%         if TowA==TowB
%             %fprintf(fout,'Delta t:%f\n', TowA-TowB ) ;
%             k = p ;
%             break ;
%         end
%     end
    if k>0
        % corresponding measurment found for receiver 2
        recvB_msr = measurmentsB{k} ;
        sat_list = [] ;
        prMesA = [] ;
        prMesB = [] ;
        cpMesA = [] ;
        cpMesB = [] ;
        Eph = [] ;
        NUMSAT = 0 ;
        % proceed with sattelite list
        for g=1:length(recvA_msr)
            s1 = -1 ;
            for s=1:length(recvB_msr)
                if recvA_msr{g}.svId==recvB_msr{s}.svId
                    % found the same sattelite for receiver B
                    s1 = s ;
                    break ;
                end
            end
            if s1>0 && NUMSAT<settings.maxSatNum
                NUMSAT = NUMSAT + 1 ;
                sat_list(NUMSAT) = recvA_msr{g}.svId ;
                prMesA(NUMSAT) = recvA_msr{g}.prMes ;
                prMesB(NUMSAT) = recvB_msr{s1}.prMes ;
                cpMesA(NUMSAT) = recvA_msr{g}.cpMes ;
                cpMesB(NUMSAT) = recvB_msr{s1}.cpMes ;
                Eph(:,NUMSAT) = eph2easy(recvA_msr{g}.s_eph, NUMSAT ) ;
            end
        end
        if NUMSAT>=settings.minSatNum
            BL = BL + 1 ;
            
            measurementTime = TowA ;
            
            L1_state = nlib_L1_baseline(fout, L1_state, settings, measurementTime, sat_list, Eph, prMesA, prMesB, cpMesA, cpMesB ) ;

            baseline_data(1, BL) = measurementTime ;
            baseline_data(2:4, BL) = L1_state.x(1:3) ;
            
            x_data(1:length(L1_state.x), BL)  = L1_state.x ;
            P_data(1:length(L1_state.x), BL) = diag(L1_state.P) ;

            z_data(:, BL) = L1_state.z(1:8) ;
        end
        
    end
end

baseline_data = baseline_data(:,1:BL) ;
x_data = x_data(:,1:BL) ;
P_data = P_data(:,1:BL) ;
z_data = z_data(:,1:BL) ;
