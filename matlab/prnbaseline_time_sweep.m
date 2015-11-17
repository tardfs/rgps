%function prnbaseline_time_sweep()
clc ;

fout = fopen('receiver_out.txt','w+t') ;
htmlFile = 'RS_matv_50mm.html'
recv1File = '..\data\RS_matv_50mm_01.mat' ;
recv2File = '..\data\RS_matv_50mm_02.mat' ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
ubx_ecef_A = ubxEcef ;
ubx_geodetic_A = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver B data...') ;
load(recv2File) ;
recv2_measurments = measurments_queue ;
ubx_ecef_B = ubxEcef ;
ubx_geodetic_B = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ;
fprintf('\n') ;

%get A data
gnss_ecef_A = nlib_gnss_ecef_solver(recv1_measurments) ;
easy_ecef_A = nlib_easy_ecef_solver(recv1_measurments) ;
% averaging poistion
ubx_ecef_mA = nlib_mean_ecef(ubx_ecef_A) ;
gnss_ecef_mA = nlib_mean_ecef(gnss_ecef_A) ;
easy_ecef_mA = nlib_mean_ecef(easy_ecef_A) ;
% geodetic position
ubx_geodetic_mA = nlib_ecef2geodetic(ubx_ecef_mA) ;
gnss_geodetic_A = nlib_ecef2geodetic(gnss_ecef_A) ;
easy_geodetic_A = nlib_ecef2geodetic(easy_ecef_A) ;
gnss_geodetic_mA = nlib_ecef2geodetic(gnss_ecef_mA) ;
easy_geodetic_mA = nlib_ecef2geodetic(easy_ecef_mA) ;


%get B data
gnss_ecef_B = nlib_gnss_ecef_solver(recv2_measurments) ;
easy_ecef_B = nlib_easy_ecef_solver(recv2_measurments) ;
% averaging poistion
ubx_ecef_mB = nlib_mean_ecef(ubx_ecef_B) ;
gnss_ecef_mB = nlib_mean_ecef(gnss_ecef_B) ;
easy_ecef_mB = nlib_mean_ecef(easy_ecef_B) ;
% geodetic position
ubx_geodetic_mB = nlib_ecef2geodetic(ubx_ecef_mB) ;
gnss_geodetic_B = nlib_ecef2geodetic(gnss_ecef_B) ;
easy_geodetic_B = nlib_ecef2geodetic(easy_ecef_B) ;
gnss_geodetic_mB = nlib_ecef2geodetic(gnss_ecef_mB) ;
easy_geodetic_mB = nlib_ecef2geodetic(easy_ecef_mB) ;

% baseline
ubx_baseline_m = ubx_ecef_mB - ubx_ecef_mA ;
gnss_baseline_m = gnss_ecef_mB - gnss_ecef_mA ;
easy_baseline_m = easy_ecef_mB - easy_ecef_mA ;
ubx_baseline = nlib_coarse_baseline_solver(ubx_ecef_A, ubx_ecef_B) ;
easy_baseline = nlib_coarse_baseline_solver(easy_ecef_A, easy_ecef_B) ;

% baseline statistics
%stat_num_ubx_bl = size(ubx_baseline,2) ;
%stat_num_easy_bl = size(easy_baseline,2) ;
%stat_num_prng_bl = size(prng_baseline,2) ;

figure(1) ;
hold off ;

Colorv = [0.1 0.1 0.8] ;
n = 0 ;
x_var = zeros(10000,1) ;
x_mean = zeros(10000,1) ;
for t_delta=-6e-3:25e-5:9e-3
    
    prng_baseline = nlib_pseudorange_baseline_solver( easy_ecef_mA, recv1_measurments, recv2_measurments, t_delta ) ;
    prng_baseline_m = nlib_mean_ecef( prng_baseline ) ;
    
    n = n+1 ;
    
    x_mean(n) = mean(prng_baseline(2,:)) ;
    x_var(n) = var(prng_baseline(2,:)) ;


    % plot data
    fprintf('\n') ;
    %fprintf('<ubx_baseline_m:> <ECEF>%f,%f,%f <R>%f\n', ubx_baseline_m(2),ubx_baseline_m(3), ubx_baseline_m(4), norm(ubx_baseline_m(2:4) )) ;
    %fprintf('<gnss_baseline_m:> <ECEF>%f,%f,%f <R>%f\n', gnss_baseline_m(2),gnss_baseline_m(3), gnss_baseline_m(4), norm(gnss_baseline_m(2:4) )) ;
    %fprintf('<easy_baseline_m:> <ECEF>%f,%f,%f <R>%f\n', easy_baseline_m(2),easy_baseline_m(3), easy_baseline_m(4), norm(easy_baseline_m(2:4) )) ;
    fprintf('<prng_time_corr0:> %f\n', t_delta ) ;
    fprintf('<prng_baseline_m:> <ECEF>%f,%f,%f <R>%f\n', prng_baseline_m(2),prng_baseline_m(3),prng_baseline_m(4), norm(prng_baseline_m(2:4) )) ;
    fprintf('<prng_x_mean_val:> <%f>\n', x_mean(n) ) ;
    fprintf('<prng_x_variance:> <%f>\n', x_var(n) ) ;

    %hold off ;
    if size(prng_baseline,2)>1
        plot((prng_baseline(1,:)-prng_baseline(1,1))/60,prng_baseline(2,:),'Color',Colorv,'LineWidth',2) ;
        hold on ;
        %plot((prng_baseline(1,:)-prng_baseline(1,1))/60,prng_baseline(3,:),'Color',[0.1 0.7 0.1],'LineWidth',2) ;
        %hold on ;
    end
    
    Colorv(1) = Colorv(1)+0.05 ;
    Colorv(2) = Colorv(2)+0.05 ;
    Colorv(3) = Colorv(3)-0.04 ;
    Colorv(Colorv>1.0) = 1.0 ;
    Colorv(Colorv<0.0) = 0.0 ;
    
    drawnow ;

end

x_var = x_var(1:n) ;
x_mean = x_mean(1:n) ;

%hold off ;
% if size(prng_baseline,2)>1
%     plot((prng_baseline(1,:)-prng_baseline(1,1))/60,nlib_ecef_norm2(prng_baseline),'Color',[0.1 0.1 0.8],'LineWidth',2) ;
%     hold on ;
% end
% if size(ubx_baseline,2)>1
%     plot((ubx_baseline(1,:)-ubx_baseline(1,1))/60,nlib_ecef_norm2(ubx_baseline),'Color',[0.1 0.7 0.1],'LineWidth',2) ;
% end
% hold on ;
% if size(easy_baseline,2)>1
%     plot((easy_baseline(1,:)-easy_baseline(1,1))/60,nlib_ecef_norm2(easy_baseline),'Color',[0.7 0.2 0.2],'LineWidth',2) ;
% end

grid on ;
set(gca,'FontSize',14) ;
xlabel('sec') ;
ylabel('baseline, m') ;
title('Baseline') ;

figure(2) ;
hold off ;
plot(-6e-3:25e-5:9e-3,x_mean,'Color',[0.1 0.1 0.8],'LineWidth',2) ;
hold on ;
grid on ;
%plot(-6e-3:25e-5:9e-3,sqrt(x_var),'Color',[0.1 0.7 0.1],'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('time correction, sec') ;
%ylabel('ECEF X, m') ;
title('X mean coord') ;

fclose(fout) ;

