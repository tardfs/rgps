clc ;
opengl software ;

settings.easyLib = getFullPath('..\\easy') ;
settings.gnssLib = getFullPath('..\\softgnss') ;
%settings.recv1File = '..\data\RS_matv_50mm_01.mat' ;
%settings.recv2File = '..\data\RS_matv_50mm_01.mat' ;
settings.recv1File = '..\data\SITE247J.mat' ;
settings.recv2File = '..\data\SITE24~1.mat' ;
[~,settings.fnameA] = fileparts( settings.recv1File ) ;
[~,settings.fnameB] = fileparts( settings.recv2File ) ;
settings.v_light = 299792458 ;	     % vacuum speed of light m/s
settings.f1 = 154*10.23e6 ;		     % L1 frequency, Hz
settings.lambda1 = settings.v_light/settings.f1 ;   % wavelength on L1:  .19029367  m
settings.check_for_time_sync = 1 ;
settings.useSolBasedPrMes = 1 ;
settings.timeSyncTol = 10e-3 ; % Synchronization tolerance
settings.minSatElevation = 15 ;
settings.minSatNum = 5 ;
settings.maxSatNum = 20 ;
settings.drawStateCovariances = 1 ;
%settings.enableSvId = [1,4,8,14,19,22,32] ; % matv_50mm 1040..1140
settings.enableSvId = [1,4,13,20,24,25] ; % matv_50mm 900..1200
%settings.enableSvId = [1,4,8,10,11,18,32] ; % matv_50mm 40..140
%settings.enableSvId = [1,4,8,14,22,32] ; % matv_50mm 800..1000
%settings.enableSvId = [4,8,11,14,19,22,32] ; % RS_matv_1400mm_680mm 50..150
%settings.enableSvId = [1,3,4,11,14,19,32] ; % RS_matv_1400mm_680mm 1460..1600
%settings.enableSvId = [15,16,18,21,22,27] ; % RS_matv_2000mm 570..720

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
fprintf(repmat('\b',1,160)) ; fprintf('get precise A position...') ;
easy_ecef_A = nlib_easy_ecef_solver( measurments_A ) ;
mean_ecef_A = nlib_mean_ecef(easy_ecef_A) ;
%mean_ecef_A(2) = mean_ecef_A(2) - 10 ;
%mean_ecef_A(3) = mean_ecef_A(3) - 10 ;
%mean_ecef_A(4) = mean_ecef_A(4) - 10 ;
fprintf(repmat('\b',1,160)) ;
fprintf('\n') ;

fprintf(repmat('\b',1,160)) ; fprintf('get precise B position...') ;
easy_ecef_B = nlib_easy_ecef_solver( measurments_B ) ;
fprintf(repmat('\b',1,160)) ;
fprintf('\n') ;


%cbaseline_data = nlib_coarse_baseline_solver(ubx_ecef_A, ubx_ecef_B) ;
[baseline_data, x_data, P_data, z_data, H_data, phi_resid_data ] = nlib_L1_baseline_solver(settings, fout, measurments_A, measurments_B, mean_ecef_A, easy_ecef_A, easy_ecef_B ) ;

figure(1) ;
% hold off ,
% plot((cbaseline_data(1,:)-cbaseline_data(1,1))/60, cbaseline_data(2,:)) ;
% hold on,
% plot((baseline_data(1,:)-baseline_data(1,1))/60, baseline_data(2,:),'r-') ;
hold off, 
plot(baseline_data(2:4,:).','LineWidth',2) ;
set(gca,'FontSize',14) ;
ylabel('m') ;
grid on ;
title(['ECEF x-y-z: ',settings.fnameA,'-', settings.fnameB],'interpreter','none') ;
xlabel('sec') ;
legend('x','y','z') ;

figure(2) ;
hold off, plot(P_data(1:3,:).', 'LineWidth',2 ), 
set(gca,'FontSize',14) ;
legend('\sigma_x^2','\sigma_y^2','\sigma_z^2') ;
grid on ;
ylabel('m^2','interpreter','tex') ;
xlabel('Epoch #','interpreter','none') ;
title([settings.recv1File,', ', settings.recv2File], 'interpreter', 'none' ) ;

figure(3) ;
f1 = 154*10.23e6 ;		     % L1 frequency, Hz
v_light = 299792458 ;	     % vacuum speed of light m/s
lambda1 = v_light/f1 ;	     % wavelength on L1:  .19029367  m
hold off, plot(P_data(4:end,:).'*settings.lambda1, 'LineWidth',2 ), 
set(gca,'FontSize',14) ;
legend(...
    '\sigma^2_{\Delta\nabla\phi^{12}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{13}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{14}_{AB}}', ...
    '\sigma^2_{\Delta\nabla\phi^{15}_{AB}}' ...
    ) ;
ylabel('m^2','interpreter','tex') ;
xlabel('Epoch #','interpreter','none') ;
grid on ;
title([settings.recv1File,', ', settings.recv2File], 'interpreter', 'none') ;

figure(4) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Measured (Input) Code Double Difference' ) ;
numDif = size(z_data,1)/2 ;
hold off, plot(z_data(1:numDif,:)','LineWidth',1) ;
set(gca,'FontSize',14) ;
xlabel('epoch #','interpreter','none') ;
ylabel('m') ;
title(sprintf('Std(\\Delta\\nabla\\rho^{ij_{AB}}):%5.3f m', mean(std(z_data(1:numDif,:),[],2))),'interpreter','tex') ;
legend('\Delta\nabla\rho^{12}_{AB}','\Delta\nabla\rho^{13}_{AB}',...
    '\Delta\nabla\rho^{14}_{AB}','\Delta\nabla\rho^{15}_{AB}','\Delta\nabla\rho^{16}_{AB}') ;
grid on ;

figure(5) ;
minutesScale = (baseline_data(1,:)-baseline_data(1,1))/60 ;
numDif = size(z_data,1)/2 ;
visData = z_data(numDif+1:end,:).' ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Measured (Input) Phase Double Differences' ) ;
hold off, plot(minutesScale,visData*settings.lambda1,'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('minutes') ;
ylabel('\Delta\nabla\phi^{ij_{AB}}, m','interpreter','tex') ;
std_phi = mean(std(visData,[],1)) ;
title(sprintf('Std(\\Delta\\nabla\\phi^{ij_{AB}}):%5.3f cycles (%5.3f m)', std_phi, std_phi*settings.lambda1 ),'interpreter','tex') ;
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

figure(7) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Estimated distance using ECEF' ) ;
hold off, 
plot((baseline_data(1,:)-baseline_data(1,1))/60, nlib_ecef_norm2(baseline_data),'r-','LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('time, sec') ;
ylabel('ECEF Distance, m') ;
grid on ;
title(['ECEF D: ',settings.fnameA,'-', settings.fnameB],'interpreter','none') ;

figure(8) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Estimated distance using UTM' ) ;
plot( (baseline_data(1,:)-baseline_data(1,1))/60, sqrt(utm_data(2,:).^2+utm_data(3,:).^2), 'LineWidth',2 ) ;
set(gca,'FontSize',14) ;
xlabel('time, sec') ;
ylabel('UTM Distance, m') ;
grid on ;
title(['UTM D: ',settings.fnameA,'-', settings.fnameB],'interpreter','none') ;

figure(9) ;
minutesScale = (baseline_data(1,:)-baseline_data(1,1))/60 ;
visData = x_data(4:end,:).' ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Filtered & resolved Phase Double Differences' ) ;
hold off, plot(minutesScale,visData*lambda1,'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('minutes') ;
ylabel('m') ;
std_phi = mean(std(visData,[],1)) ;
title(sprintf('Std(\\Delta\\nabla\\phi^{ij_{AB}}):%5.3f cycles (%5.3f m)', std_phi, std_phi*lambda1 )) ;
legend('\Delta\nabla\phi^{12}_{AB}','\Delta\nabla\phi^{13}_{AB}',...
    '\Delta\nabla\phi^{14}_{AB}','\Delta\nabla\phi^{15}_{AB}') ;
grid on ;

figure(10) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'Phase residuals' ) ;
visData = phi_resid_data.' ;
hold off, plot( visData ,'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('Epoch #') ;
ylabel('m') ;
title('\Delta\nabla\phi^{ij}_{AB}-N^{ij}_{AB}', 'interpreter', 'tex') ;
grid on ;

figure(11) ;
set(gcf,'NumberTitle','off') ;
set(gcf,'Name', 'H_Data' ) ;
visData = H_data.' ;
hold off, plot( visData ,'LineWidth',2) ;
set(gca,'FontSize',14) ;
xlabel('Epoch #') ;
ylabel('m') ;
grid on ;

fclose(fout) ;
