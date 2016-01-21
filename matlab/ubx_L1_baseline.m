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

baseline_data = nlib_L1_baseline_solver(fout, measurments_A, measurments_B) ;

fclose(fout) ;
