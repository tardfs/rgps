function view_receiver()
clc ;
easyLib = getFullPath('..\\easy') ;

fout = fopen('receiver_out.txt','w+t') ;
recv1File = '..\data\RS_twin_30mm_01.mat' ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
ubxEcefA = ubxEcef ;
ubxGeodeticA = ubxGeodetic ;

% Prepare navigation data
N = length(recv1_measurments) ;

% get averaged ubx
ubx_mA = mean(ubxEcefA(2:4,:),2) ;
% map to geoid
addpath(easyLib) ;
[ubx_mA_phi,ubx_mA_lambda,ubx_mA_h] = togeod(6378137,298.257223563,ubx_mA(1),ubx_mA(2),ubx_mA(3)) ;
rmpath(easyLib) ;

% compute position using different methods
for idx1=1:N
    recv1_msr_set = recv1_measurments{idx1} ;
    [ux,uy,uz,gdop] = gnss_pvt_solver( fout, recv1_msr_set) ;
    gnssApos(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop(1)] ;
    [ux,uy,uz,gdop] = easy_pvt_solver( fout, recv1_msr_set) ;
    easyApos(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop] ;
end

% get averaged gnss
gnss_mA = mean(gnssApos(2:4,:),2) ;
% map to geoid
addpath(easyLib) ;
[gnss_mA_phi,gnss_mA_lambda,gnss_mA_h] = togeod(6378137,298.257223563,gnss_mA(1),gnss_mA(2),gnss_mA(3)) ;
rmpath(easyLib) ;



% get averaged easy
easy_mA = mean(easyApos(2:4,:),2) ;
% map to geoid
addpath(easyLib) ;
[easy_mA_phi,easy_mA_lambda,easy_mA_h] = togeod(6378137,298.257223563,easy_mA(1),easy_mA(2),easy_mA(3)) ;
rmpath(easyLib) ;


gnssGeodeticA = zeros(2,N) ;
easyGeodeticA = zeros(2,N) ;
addpath(easyLib) ;
for idx1=1:min(N,30)
    [gnssGeodeticA(1,idx1), gnssGeodeticA(2,idx1)] = togeod(6378137,298.257223563, gnssApos(2,idx1),gnssApos(3,idx1),gnssApos(4,idx1) ) ;
    [easyGeodeticA(1,idx1), easyGeodeticA(2,idx1)] = togeod(6378137,298.257223563, easyApos(2,idx1),easyApos(3,idx1),easyApos(4,idx1) ) ;
end
rmpath(easyLib) ;


fprintf(fout,'https://static-maps.yandex.ru/1.x/?ll=%8.6f,%8.6f&size=450,450&z=16&l=map', ubx_mA_lambda,ubx_mA_phi ) ;
fprintf(fout,'&pt=%8.6f,%8.6f,pmgns1', ubxGeodeticA(3,1),ubxGeodeticA(2,1) ) ;    

for n=1:min(N,30)
    fprintf(fout,'~%8.6f,%8.6f,pmrds', easyGeodeticA(2,n),easyGeodeticA(1,n) ) ;
end

for n=1:min(N,30)
    fprintf(fout,'~%8.6f,%8.6f,pmlbs', gnssGeodeticA(2,n),gnssGeodeticA(1,n) ) ;
end

for n=1:min(N,30)
    fprintf(fout,'~%8.6f,%8.6f,pmgns', ubxGeodeticA(3,n),ubxGeodeticA(2,n) ) ;
end

fprintf(fout,'\n\n') ;

fprintf(fout,'https://static-maps.yandex.ru/1.x/?ll=%8.6f,%8.6f&size=450,450&z=16&l=map', ubx_mA_lambda,ubx_mA_phi ) ;
fprintf(fout,'&pt=%8.6f,%8.6f,pmgnm1', ubx_mA_lambda,ubx_mA_phi ) ;
fprintf(fout,'~%8.6f,%8.6f,pmrdm2', easy_mA_lambda,easy_mA_phi ) ;
fprintf(fout,'~%8.6f,%8.6f,pmlbm3', gnss_mA_lambda,gnss_mA_phi ) ;


fclose(fout) ;
fprintf('\n') ;

function [ux,uy,uz,gdop, phi,lambda,h, el,az ] = gnss_pvt_solver( fout, measurments_queue)
    gnssLib = getFullPath('..\\softgnss') ;
    addpath(gnssLib) ;
    
    numSat = length(measurments_queue) ;
    satpos = zeros(3,numSat) ;
    obs = zeros(numSat,1) ;
    
    %sat_list = [] ;
    %gnssId = 0 ;
    %n = 0 ;
    fprintf(fout, '<GNSS PVT SOLVER SAT LIST>\n') ;
    %while length(sat_list)<numSat
    for n=1:length(measurments_queue)
        fprintf(fout, '\t<svId>%2d', measurments_queue{n}.svId ) ;
        fprintf(fout, '\t\t<SAT POS>(%7.2f,%7.2f,%7.2f)\n', measurments_queue{n}.sat_x,measurments_queue{n}.sat_y,measurments_queue{n}.sat_z) ;
        sat_distance = sqrt(measurments_queue{n}.sat_x*measurments_queue{n}.sat_x+measurments_queue{n}.sat_y*measurments_queue{n}.sat_y+measurments_queue{n}.sat_z*measurments_queue{n}.sat_z) ;
        fprintf(fout, '\t\t\t\t<SAT DISTANCE>%7.2f, <PSEUDORANGE>%7.2f\n', sat_distance, measurments_queue{n}.prMes ) ;
        satpos(1,n) = measurments_queue{n}.sat_x ;
        satpos(2,n) = measurments_queue{n}.sat_y ;
        satpos(3,n) = measurments_queue{n}.sat_z ;
        obs(n) = measurments_queue{n}.prMes + measurments_queue{n}.sat_clk_corr*299792458 ;
    end
    
    settings.c    = 299792458 ;    % The speed of light, [m/s]
    settings.useTropCorr = 1 ;     % Use troposphere correction
    [pos, el, az, dop] = leastSquarePos(satpos, obs, settings) ;
    
    [phi, lambda, h] = cart2geo(pos(1), pos(2), pos(3), 5 ) ;
    
    ux = pos(1) ;
    uy = pos(2) ;
    uz = pos(3) ;
    gdop = dop ;
    fprintf(fout,'<GNSS NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GDOP><%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        pos(1), pos(2), pos(3), gdop(1), phi, lambda, h ) ;
    rmpath(gnssLib) ;


function [ux,uy,uz,gdop, phi,lambda,h, el,az ] = easy_pvt_solver(fout,measurments_queue)
    easyLib = getFullPath('..\\easy') ;
    addpath(easyLib) ;
    
    numSat = length(measurments_queue) ;
    obs = zeros(numSat,1) ;
    
    Eph = zeros(21,length(measurments_queue)) ;
    for n=1:length(measurments_queue)
        obs(n) = measurments_queue{n}.prMes ;        
        Eph(:,n) = eph2easy(measurments_queue{n}.s_eph, n ) ; % convert to easy ephemeris
    end
    
    [easy_pos, el, gdop] = recpo_ls(obs,(1:length(obs)),measurments_queue{1}.msrTow,Eph) ;
    [phi,lambda,h] = togeod(6378137,298.257223563,easy_pos(1),easy_pos(2),easy_pos(3)) ;

    ux = easy_pos(1) ;
    uy = easy_pos(2) ;
    uz = easy_pos(3) ;
    
    fprintf(fout,'<EASY NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GDOP><%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        easy_pos(1), easy_pos(2), easy_pos(3), gdop, phi, lambda, h ) ;
    
    rmpath(easyLib) ;


