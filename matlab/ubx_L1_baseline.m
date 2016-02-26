clc ;

fout = fopen('ubx_L1_out.txt','w+t') ;

recv1File = '..\data\RS_matv_50mm_01.mat' ;
recv2File = '..\data\RS_matv_50mm_02.mat' ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
load(recv1File) ;
measurments_A = measurments_queue ;
ubx_ecef_A = ubxEcef ;
ubx_geodetic_A = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver B data...') ;
load(recv2File) ;
measurments_B = measurments_queue ;
ubx_ecef_B = ubxEcef ;
ubx_geodetic_B = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ;
fprintf('\n') ;

cbaseline_data = nlib_coarse_baseline_solver(ubx_ecef_A, ubx_ecef_B) ;
[baseline_data, x_data, P_data] = nlib_L1_baseline_solver(fout, measurments_A, measurments_B) ;

figure(1) ;
hold off ,
plot((cbaseline_data(1,1:1000)-cbaseline_data(1,1))/60, cbaseline_data(2,1:1000)) ;
hold on,
plot((baseline_data(1,1:1000)-baseline_data(1,1))/60, baseline_data(2,1:1000),'r-') ;
xlabel('sec') ;

figure(2) ;
hold off, plot(P_data(1:3,:).', 'LineWidth',2 ), 
set(gca,'FontSize',14) ;
legend('\sigma_x','\sigma_y','\sigma_z') ;
grid on ;
ylabel('m^2') ;
title([recv1File,', ', recv2File]) ;

figure(3) ;
f1 = 154*10.23e6 ;		     % L1 frequency, Hz
v_light = 299792458 ;	     % vacuum speed of light m/s
lambda1 = v_light/f1 ;	     % wavelength on L1:  .19029367  m
hold off, plot(P_data(4:end,:).'*lambda1, 'LineWidth',2 ), 
set(gca,'FontSize',14) ;
legend(...
    '\sigma^2_{\Delta\nabla\phi^{12}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{13}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{14}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{15}_{AB}}' ...
    ) ;
ylabel('m^2') ;
grid on ;
title([recv1File,', ', recv2File]) ;


fclose(fout) ;
