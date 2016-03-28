clc ;

settings.easyLib = getFullPath('..\\easy') ;
settings.gnssLib = getFullPath('..\\softgnss') ;
settings.recv1File = '..\data\RS_matv_50mm_01.mat' ;
settings.recv2File = '..\data\RS_matv_50mm_02.mat' ;
settings.v_light = 299792458 ;	     % vacuum speed of light m/s
settings.check_for_time_sync = 1 ;
settings.timeSyncTol = 1e-3 ; % Synchronization tolerance
settings.minSatNum = 5 ;
settings.maxSatNum = 20 ;
settings.enableSvId = [1,4,8,14,19,22,32] ;
[~,settings.fnameA] = fileparts( settings.recv1File ) ;
[~,settings.fnameB] = fileparts( settings.recv2File ) ;


fout = fopen('ubx_L1_out.txt','w+t') ;

fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
load(settings.recv1File) ;
measurments_A = measurments_queue ;
ubx_ecef_A = ubxEcef ;
ubx_geodetic_A = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver B data...') ;
load(settings.recv2File) ;
measurments_B = measurments_queue ;
ubx_ecef_B = ubxEcef ;
ubx_geodetic_B = ubxGeodetic ;
fprintf(repmat('\b',1,160)) ;
fprintf('\n') ;

% get precise A position
mean_ecef_A = nlib_mean_ecef(nlib_easy_ecef_solver( measurments_A )) ;

cbaseline_data = nlib_coarse_baseline_solver(ubx_ecef_A, ubx_ecef_B) ;
[baseline_data, x_data, P_data, z_data] = nlib_L1_baseline_solver(settings, fout, measurments_A, measurments_B) ;

figure(1) ;
hold off ,
plot((cbaseline_data(1,:)-cbaseline_data(1,1))/60, cbaseline_data(2,:)) ;
hold on,
plot((baseline_data(1,:)-baseline_data(1,1))/60, baseline_data(2,:),'r-') ;
hold off, 
plot((baseline_data(1,:)-baseline_data(1,1))/60, nlib_ecef_norm2(baseline_data),'r-') ;
hold off, 
plot(baseline_data(2:4,:).','LineWidth',2) ;
set(gca,'FontSize',14) ;
ylabel('m') ;
grid on ;
title(['ECEF x-y-z: ',settings.fnameA,'<->', settings.fnameB],'interpreter','none') ;
xlabel('sec') ;

figure(2) ;
hold off, plot(P_data(1:3,:).', 'LineWidth',2 ), 
set(gca,'FontSize',14) ;
legend('\sigma_x','\sigma_y','\sigma_z') ;
grid on ;
ylabel('m^2') ;
title([settings.recv1File,', ', settings.recv2File], 'interpreter', 'none' ) ;

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
title([settings.recv1File,', ', settings.recv2File]) ;

figure(4) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Code Double Difference' ) ;
numDif = size(z_data,1)/2 ;
hold off, plot(z_data(1:numDif,:)','LineWidth',1) ;
set(gca,'FontSize',14) ;
xlabel('epoch #') ;
ylabel('\Delta\nabla\rho^{ij_{AB}}') ;
title(sprintf('Std(\\Delta\\nabla\\rho^{ij_{AB}}):%5.3f m', mean(std(z_data(1:4,:),[],2)))) ;
legend('\Delta\nabla\rho^{12}_{AB}','\Delta\nabla\rho^{13}_{AB}',...
    '\Delta\nabla\rho^{14}_{AB}','\Delta\nabla\rho^{15}_{AB}') ;
grid on ;

figure(5) ;
minutesScale = (baseline_data(1,:)-baseline_data(1,1))/60 ;
numDif = size(z_data,1)/2 ;
visData = z_data(numDif+1:end,:).' ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Phase Double Differences' ) ;
hold off, plot(minutesScale,visData*lambda1,'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('minutes') ;
ylabel('\Delta\nabla\phi^{ij_{AB}}, m') ;
std_phi = mean(std(visData,[],1)) ;
title(sprintf('Std(\\Delta\\nabla\\phi^{ij_{AB}}):%5.3f cycles (%5.3f m)', std_phi, std_phi*lambda1 )) ;
legend('\Delta\nabla\phi^{12}_{AB}','\Delta\nabla\phi^{13}_{AB}',...
    '\Delta\nabla\phi^{14}_{AB}','\Delta\nabla\phi^{15}_{AB}') ;
grid on ;

figure(6) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'UTM relative coordinates' ) ;
addpath(settings.easyLib) ;
utm_data = nlib_ecef_to_utm(mean_ecef_A, baseline_data) ;
rmpath(settings.easyLib) ;
plot( utm_data(2,:), utm_data(3,:), 'LineWidth',2 ) ;
set(gca,'FontSize',14) ;
xlabel('UTM east, m') ;
ylabel('UTM north, m') ;
grid on ;
title(['UTM East-North: ',settings.fnameA,'-', settings.fnameB],'interpreter','none') ;

fclose(fout) ;
