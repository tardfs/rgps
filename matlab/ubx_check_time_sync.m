settings.easyLib = getFullPath('..\\easy') ;
settings.gnssLib = getFullPath('..\\softgnss') ;
settings.recv1File = '..\data\RS_twin_300mm_01.mat' ;
settings.recv2File = '..\data\RS_twin_300mm_02.mat' ;
[~,settings.fnameA] = fileparts( settings.recv1File ) ;
[~,settings.fnameB] = fileparts( settings.recv2File ) ;
settings.v_light = 299792458 ;	     % vacuum speed of light m/s
settings.f1 = 154*10.23e6 ;		     % L1 frequency, Hz
settings.lambda1 = settings.v_light/settings.f1 ;   % wavelength on L1:  .19029367  m

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

% get precise A position
fprintf(repmat('\b',1,160)) ; fprintf('get precise B position...') ;
easy_ecef_B = nlib_easy_ecef_solver( measurments_B ) ;



hold off ,
n = 27 ;
N = 2700 ;
plot( (easy_ecef_A(1,1:N) - easy_ecef_A(5,1:N)/settings.v_light)-(easy_ecef_B(1,1+n:N+n)- easy_ecef_B(5,1+n:N+n)/settings.v_light), 'LineWidth', 1 ) ;

hold on,
plot( (easy_ecef_A(1,1:N))-(easy_ecef_B(1,1+n:N+n)),'-', 'LineWidth', 1 ) ;

legend('''sattelite'' time difference from receivers point of view','RXM-RAWX rcvTow difference') ;
grid on ;
xlabel('seconds') ;
ylabel('Time Mismatch, seconds') ;
set(gca,'FontSize',14) ;


