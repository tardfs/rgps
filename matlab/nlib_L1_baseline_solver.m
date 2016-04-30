% Baseline solver using L1 measurments only
function [baseline_data, x_data, P_data, z_data, H_data, phi_resid_data] = nlib_L1_baseline_solver(settings, fout, measurmentsA, measurmentsB, mean_ecef_A, easy_ecef_A)

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
H_data = zeros(32,min(N,K)) ;
phi_resid_data = zeros(12,min(N,K)) ;


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
    %plot((recvA_Tow_v-recvA_Tow_v(1))/60, AB_time_err*settings.v_light,'LineWidth',2) ;
    plot(AB_time_err*settings.v_light,'LineWidth',2) ;
    set(gca,'FontSize',14) ;
    grid on ;
    xlabel('Epoch #') ;
    %ylabel(sprintf('Light travel, m (1ms:%3.2f m)', settings.v_light*1e-3 )) ;
    ylabel('ms') ;
    title(sprintf('A-B Time mismatch: %s', settings.fnameA ),'interpreter','none') ;
    
    figure(21) ;
    hold off , plot(recvA_Tow_v(2:end)-recvA_Tow_v(1:end-1), 'LineWidth',2) ;
    hold on ,  plot(recvB_Tow_v(2:end)-recvB_Tow_v(1:end-1), 'r-','LineWidth',2) ;
    title(sprintf('A & B Time step: %s', settings.fnameA ),'interpreter','none') ;
    legend('A time step', 'B time step') ;
    set(gca,'FontSize',14) ;
    grid on ;
    
    c=input('Proceed with relative coordinates: [Y/N]?','s') ;
    if isempty(c)
        c = 'Y' ;
    end
    if c~='Y' && c~='y'
        return ;
    end
    
    startIndex = input('Index from:') ;
    if startIndex<1 || startIndex>N
        startIndex = 1 ;
    end
    endIndex = input('Index to:') ;
    if endIndex<1 || endIndex>N
        endIndex = N ;
    end
    
    figure(22),
    clf,
    set(gcf,'NumberTitle','off') ;
    set(gcf,'Name', '[RGPS] Sattelite Sky Plot' ) ;

    axes() ; set(gca,'Position', [.03 .52 .42 .42] ) ;
    [~, ~, ~, hlMapA, hl_sat_list_A] = nlib_plot_sats(settings, measurmentsA(startIndex:endIndex)) ;
    title( sprintf('%s,range: %d..%d', settings.fnameA, startIndex, endIndex ),'interpreter','none','FontSize',14)  ;
    axes() ; set(gca,'Position', [.03 .04 .44 .42] ) ;
    %[hlMapA, hl_sat_list_A] = nlib_highlight_smap(settings, sMapA, ElA) ;
    nlib_sat_timeline(hlMapA)
    title(sprintf('%1d,',hl_sat_list_A),'FontSize',14) ;

    axes() ; set(gca,'Position', [.53 .52 .42 .42] ) ;    
    [~, ~, ~, hlMapB, hl_sat_list_B] = nlib_plot_sats(settings, measurmentsB(AB_n2k(startIndex:endIndex))) ;
    title( sprintf('%s,range: %d..%d', settings.fnameB, AB_n2k(startIndex), AB_n2k(endIndex) ),'interpreter','none','FontSize',14 ) ;
    axes() ; set(gca,'Position', [.53 .04 .44 .42] ) ;
    %[hlMapB, hl_sat_list_B] = nlib_highlight_smap(sMapB, ElB) ;
    nlib_sat_timeline(hlMapB)
    title(sprintf('%1d,',hl_sat_list_B),'FontSize',14) ;
    
    mean_ecef_A = nlib_mean_ecef(easy_ecef_A(:,startIndex:endIndex)) ;
    
    sat_list = intersect( hl_sat_list_A, hl_sat_list_B ) ;
    
    fprintf('Proposed sat list ') ; fprintf('%1d,', sat_list) ;
    c = input('Accept [Y/N]?','s') ;
    if isempty(c)
        c = 'Y' ;
    end
    if c~='Y' && c~='y'
        return ;
    end
    settings.enableSvId = sat_list ;
    
end

for n=startIndex:endIndex    
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
            if nnz(recvA_msr{g}.svId==settings.enableSvId)>0
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
        end
        if NUMSAT>=settings.minSatNum
            BL = BL + 1 ;
            
            measurementTime = TowA ;
            
            L1_state = nlib_L1_baseline(fout, L1_state, settings, measurementTime, sat_list, Eph, prMesA, prMesB, cpMesA, cpMesB, mean_ecef_A ) ;

            baseline_data(1, BL) = measurementTime ;
            %baseline_data(2:4, BL) = L1_state.x(1:3) ;
            baseline_data(2:4, BL) = L1_state.b ;
            
            x_data(1:length(L1_state.x), BL)  = L1_state.x ;
            P_data(1:length(L1_state.x), BL) = diag(L1_state.P) ;

            z_data(:, BL) = L1_state.z(1:8) ;
            hData = L1_state.H(1:NUMSAT-1,1:3) ;
            H_data(1:(NUMSAT-1)*3, BL ) = hData(:) ;
            phi_resid_data(1:NUMSAT-1, BL) = L1_state.phi_resid ;
            
            if settings.drawStateCovariances
                figure(23) ;
                set(gcf,'NumberTitle','off') ;
                set(gcf,'Name', sprintf('State vector covariance matrix %d',n) ) ;
                contourf(L1_state.P) ;
                title('State vector covariance matrix') ;
                set(gca,'XTickLabel',{'x','y','z','N12','N13','N14','N15','N16','N17','N18','N19'}) ;
                set(gca,'YTickLabel',{'x','y','z','N12','N13','N14','N15','N16','N17','N18','N19'}) ;
                axis square ;
                drawnow ;
            end
            
        end
        
    end
end

baseline_data = baseline_data(:,1:BL) ;
x_data = x_data(:,1:BL) ;
P_data = P_data(:,1:BL) ;
z_data = z_data(:,1:BL) ;
H_data = H_data(1:(NUMSAT-1)*3,1:BL) ;
phi_resid_data = phi_resid_data(1:NUMSAT-1,1:BL) ;
