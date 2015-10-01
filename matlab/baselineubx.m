function baselineubx()
easyLib = getFullPath('..\\easy') ;

clc ;
fout = fopen('rltv_out.txt','w+t') ;
recv1File = '..\data\RS_twin11_01.mat' ;
recv2File = '..\data\RS_twin11_02.mat' ;
outFile = '..\data\rltv_out.mat' ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
ubxEcefA = ubxEcef ;
fprintf(repmat('\b',1,160)) ; fprintf('load receiver B data...') ;
load(recv2File) ;
ubxEcefB = ubxEcef ;
recv2_measurments = measurments_queue ;

if 0
    % limit dataset
    N = 1000 ;
    if length(recv1_measurments)>2*N
        recv1_measurments = recv1_measurments(1,1:N) ;
        recv2_measurments = recv2_measurments(1,1:N) ;
    end
end

N = length(recv1_measurments) ;

% position format
% X(1,:) - TOW, sec
% X(2,:) - ECEF X, m
% X(3,:) - ECEF Y, m
% X(4,:) - ECEF Z, m
% X(5,:) - ECEF DOP, m

gnssApos = zeros(5,N) ; 
gnssBpos = zeros(5,N) ;
easyApos = zeros(5,N) ;
easyBpos = zeros(5,N) ;
numBaselines = 0 ;
coarseBaseline = zeros(5,N) ;
%dcosBaseline = zeros(5,N) ;
prBaseline = zeros(5,N) ;

fprintf(repmat('\b',1,160)) ; fprintf('Initial processing...') ;


for idx1=1:N
    recv1_msr_set = recv1_measurments{idx1} ;
    [ux,uy,uz,gdop] = gnss_pvt_solver( fout, recv1_msr_set) ;
    gnssApos(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop(1)] ;
    [ux,uy,uz,gdop] = easy_pvt_solver( fout, recv1_msr_set) ;
    easyApos(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop] ;
end
% best possible receiver A position
master_pos = mean(easyApos(2:4,:),2) ;
addpath(easyLib) ;
[m_phi,m_lambda,m_h] = togeod(6378137,298.257223563,master_pos(1),master_pos(2),master_pos(3)) ;
rmpath(easyLib) ;

fprintf(fout,'<MASTER NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n\n', ...
        master_pos(1), master_pos(2), master_pos(3), m_phi, m_lambda, m_h ) ;
    
for idx2=1:length(recv2_measurments)
    recv2_msr_set = recv2_measurments{idx2} ;
    [ux,uy,uz,gdop] = gnss_pvt_solver( fout, recv2_msr_set) ;
    gnssBpos(:,idx2) = [recv2_msr_set{1}.msrTow; ux;uy;uz; gdop(1)] ;
    [ux,uy,uz,gdop] = easy_pvt_solver( fout, recv2_msr_set) ;
    easyBpos(:,idx2) = [recv2_msr_set{1}.msrTow; ux;uy;uz; gdop] ;
end
% best possible receiver B position
B_pos = mean(easyBpos(2:4,:),2) ;
addpath(easyLib) ;
[B_phi,B_lambda,B_h] = togeod(6378137,298.257223563,B_pos(1),B_pos(2),B_pos(3)) ;
rmpath(easyLib) ;
fprintf(fout,'<ROVER NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n\n', ...
        B_pos(1), B_pos(2), B_pos(3), B_phi, B_lambda, B_h ) ;

fprintf(fout,'<AVG BASELINE><%f,%f,%f>\n\n', ...
        B_pos(1)-master_pos(1), B_pos(2)-master_pos(2), B_pos(3)-master_pos(3) ) ;
    
% Baseline algorithm evaluation
n_show = 0 ;
for idx1=1:N
    
    if idx1>n_show
        fprintf(repmat('\b',1,160)) ; 
        fprintf( '<MEASURMENT IDX>%4d/%d <BASELINES:>%d', idx1, N, numBaselines ) ;
        n_show = idx1 + 20 ;
    end
        
    recv1_msr_set = recv1_measurments{idx1} ;
    % get rcv1 measurement time
    rcv1Tow = recv1_msr_set{1}.msrTow ;
    % find corresponded measurment for recv2
    idx2 = -1 ;
    for n=1:length(recv2_measurments)
        recv2_msr_set = recv2_measurments{n} ;
        rcv2Tow = recv2_msr_set{1}.msrTow ;
        if rcv1Tow==rcv2Tow
            idx2 = n ;
            break ;
        end
    end
    if idx2>0
        % corresponding measurment found for receiver 2
        numBaselines = numBaselines + 1 ;
        
        recv2_msr_set = recv2_measurments{idx2} ;        
        fprintf(fout, '<DIF:> TOW_A:%f TOW_B:%f\n', recv1_msr_set{1}.msrTow, ...
            recv2_msr_set{1}.msrTow ) ;
        
        [x_a,y_a,z_a] = easy_pvt_solver(fout,recv1_msr_set) ;
        [x_b,y_b,z_b] = easy_pvt_solver(fout,recv2_msr_set) ;
        
        coarseBaseline(1,numBaselines) = recv1_msr_set{1}.msrTow ;
        coarseBaseline(2,numBaselines) = x_b-x_a ;
        coarseBaseline(3,numBaselines) = y_b-y_a ;
        coarseBaseline(4,numBaselines) = z_b-z_a ;
        coarseBaseline(5,numBaselines) = 0 ;
        
        fprintf(fout, '<COARSE BASELINE:> <%7.2f,%7.2f,%7.2f> <R><%7.2f>\n', ...
            x_b-x_a, y_b-y_a, z_b-z_a, ...
            sqrt(coarseBaseline(2:4,numBaselines)'*coarseBaseline(2:4,numBaselines)) ) ;

        obs1 = [] ;
        obs2 = [] ;
        Eph = [] ;
        numsat = 0 ;
        % proceed with sattelite list
        for g=1:length(recv1_msr_set)
             recv1_msr = recv1_msr_set{g} ;
            % find similar sattelite for second receiver
            idx_msr = -1 ;
            for s=1:length(recv2_msr_set)
                if recv1_msr.svId==recv2_msr_set{s}.svId
                    idx_msr = s ;
                    break ;
                end
            end
            if idx_msr>0
                numsat = numsat + 1 ;
                recv2_msr = recv2_msr_set{idx_msr} ;
                obs1(numsat) = recv1_msr.prMes ;
                obs2(numsat) = recv2_msr.prMes ;
                Eph(:,numsat) = eph2easy(recv1_msr.s_eph, numsat ) ;
            end
        end
        
        addpath(easyLib) ;
        [omc, base] = baseline(master_pos,obs1(:),obs2(:),(1:length(obs1))',rcv1Tow,Eph) ;
        rmpath(easyLib) ;
        
        prBaseline(1, numBaselines) = rcv1Tow ;
        prBaseline(2, numBaselines) = base(1) ;
        prBaseline(3, numBaselines) = base(2) ;
        prBaseline(4, numBaselines) = base(3) ;
        prBaseline(5, numBaselines) = 0 ;
    end
end

coarseBaseline = coarseBaseline(:,1:numBaselines) ;
prBaseline = prBaseline(:,1:numBaselines) ;
save(outFile,'ubxEcefA','ubxEcefB','gnssApos','gnssBpos','easyApos',...
    'easyBpos','coarseBaseline','prBaseline','master_pos','B_pos') ;
fprintf('\n') ;
fclose(fout) ;

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


