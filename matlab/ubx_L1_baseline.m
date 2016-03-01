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
[baseline_data, x_data, P_data, z_data] = nlib_L1_baseline_solver(fout, measurments_A, measurments_B) ;

figure(1) ;
hold off ,
plot((cbaseline_data(1,:)-cbaseline_data(1,1))/60, cbaseline_data(2,:)) ;
hold on,
plot((baseline_data(1,:)-baseline_data(1,1))/60, baseline_data(2,:),'r-') ;
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

figure(4) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Code Double Difference' ) ;
hold off, plot(z_data(1:4,:)') ;
set(gca,'FontSize',14) ;
xlabel('epoch #') ;
ylabel('\Delta\nabla\rho^{ij_{AB}}') ;
title(sprintf('Std(\\Delta\\nabla\\rho^{ij_{AB}}):%5.3f m', mean(std(z_data(1:4,:),[],2)))) ;
legend('\Delta\nabla\rho^{12}_{AB}','\Delta\nabla\rho^{13}_{AB}',...
    '\Delta\nabla\rho^{14}_{AB}','\Delta\nabla\rho^{15}_{AB}') ;
grid on ;

figure(5) ;
minutesScale = (baseline_data(1,:)-baseline_data(1,1))/60 ;
visData = z_data(8,:).' ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Phase Double Differences' ) ;
hold off, plot(minutesScale,visData*lambda1) ;
set(gca,'FontSize',14) ;
xlabel('minutes') ;
ylabel('\Delta\nabla\phi^{ij_{AB}}, m') ;
std_phi = mean(std(visData,[],1)) ;
title(sprintf('Std(\\Delta\\nabla\\phi^{ij_{AB}}):%5.3f cycles (%5.3f m)', std_phi, std_phi*lambda1 )) ;
legend('\Delta\nabla\phi^{15}_{AB}','\Delta\nabla\phi^{13}_{AB}',...
    '\Delta\nabla\phi^{14}_{AB}','\Delta\nabla\phi^{15}_{AB}') ;
grid on ;



fclose(fout) ;
