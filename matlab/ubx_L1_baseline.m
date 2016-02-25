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
baseline_data = nlib_L1_baseline_solver(fout, measurments_A, measurments_B) ;

hold off ,
plot((cbaseline_data(1,1:1000)-cbaseline_data(1,1))/60, cbaseline_data(2,1:1000)) ;
hold on,
plot((baseline_data(1,1:1000)-baseline_data(1,1))/60, baseline_data(2,1:1000),'r-') ;
xlabel('sec') ;

fclose(fout) ;
