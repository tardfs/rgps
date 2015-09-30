function baselineubx()
easyLib = getFullPath('..\\easy') ;

clc ;
fout = fopen('rltv_out.txt','w+t') ;
recv1File = '..\data\RS_twin12_01.mat' ;
recv2File = '..\data\RS_twin12_02.mat' ;
outFile = '..\data\rltv_out.mat' ;
load(recv1File) ;
recv1_measurments = measurments_queue ;
load(recv2File) ;
recv2_measurments = measurments_queue ;

if 1
    % limit dataset
    N = 1000 ;
    if length(recv1_measurments)>2*N
        recv1_measurments = recv1_measurments(1,1:N) ;
        recv2_measurments = recv2_measurments(1,1:N) ;
    end
end

N = length(recv1_measurments) ;

coarseApos = zeros(N,3) ;
coarseBpos = zeros(N,3) ;
coarseBaseline = zeros(N,4) ;
fineBaseline = zeros(N,4) ;

for idx1=1:N
    recv1_msr_set = recv1_measurments{idx1} ;
    [x_a,y_a,z_a] = pvt_solver(fout,recv1_msr_set) ;
    coarseApos(idx1,:) = [x_a,y_a,z_a] ;
end
% best possible receiver A position
master_pos = mean(coarseApos,1).' ;
addpath(easyLib) ;
[m_phi,m_lambda,m_h] = togeod(6378137,298.257223563,master_pos(1),master_pos(2),master_pos(3)) ;
rmpath(easyLib) ;

fprintf(fout,'<NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n\n', ...
        master_pos(1), master_pos(2), master_pos(3), m_phi, m_lambda, m_h ) ;

bases = [] ;
for idx1=1:length(recv1_measurments)
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
        recv2_msr_set = recv2_measurments{idx2} ;
        fprintf(fout, '<DIF:> TOW_A:%f TOW_B:%f\n', recv1_msr_set{1}.msrTow, ...
            recv2_msr_set{1}.msrTow ) ;
        [x_a,y_a,z_a] = pvt_solver(fout,recv1_msr_set) ;
        [x_b,y_b,z_b] = pvt_solver(fout,recv2_msr_set) ;
        
        coarseApos(idx1,:) = [x_a,y_a,z_a] ;
        coarseBpos(idx1,:) = [x_b,y_b,z_b] ;
        
        coarseBaseline(idx1,1) = x_b-x_a ;
        coarseBaseline(idx1,2) = y_b-y_a ;
        coarseBaseline(idx1,3) = z_b-z_a ;
        coarseBaseline(idx1,4) = sqrt((x_b-x_a)^2+ (y_b-y_a)^2+ (z_b-z_a)^2) ;
        
        fprintf(fout, '<COARSE BASELINE:> <%7.2f,%7.2f,%7.2f> <R><%7.2f>\n', ...
            x_b-x_a, y_b-y_a, z_b-z_a, ...
            coarseBaseline(idx1,4) ) ;

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
        
        bases = [bases base] ;
    end
end
save(outFile,'coarseBaseline','fineBaseline','coarseApos','coarseBpos','bases') ;
fclose(fout) ;


function [ux,uy,uz] = pvt_solver(fout,measurments_queue)
    easyLib = getFullPath('..\\easy') ;
    gnssLib = getFullPath('..\\softgnss') ;
    %rmpath(easyLib) ;
    %rmpath(gnssLib) ;
    
    numSat = length(measurments_queue) ;
    satpos = zeros(3,numSat) ;
    obs = zeros(numSat,1) ;
    
    %sat_list = [] ;
    %gnssId = 0 ;
    %n = 0 ;
    fprintf(fout, '<PVT SOLVER SAT LIST>\n') ;
    %while length(sat_list)<numSat
    Eph = zeros(21,length(measurments_queue)) ;
    for n=1:length(measurments_queue)
%         m_index = get_actual_measurment(measurments_queue,sat_list,gnssId) ;
%         if m_index<1
%             break ;
%         end
        fprintf(fout, '\t<svId>%2d', measurments_queue{n}.svId ) ;
        fprintf(fout, '\t\t<SAT POS>(%7.2f,%7.2f,%7.2f)\n', measurments_queue{n}.sat_x,measurments_queue{n}.sat_y,measurments_queue{n}.sat_z) ;
        sat_distance = sqrt(measurments_queue{n}.sat_x*measurments_queue{n}.sat_x+measurments_queue{n}.sat_y*measurments_queue{n}.sat_y+measurments_queue{n}.sat_z*measurments_queue{n}.sat_z) ;
        fprintf(fout, '\t\t\t\t<SAT DISTANCE>%7.2f, <PSEUDORANGE>%7.2f\n', sat_distance, measurments_queue{n}.prMes ) ;
        satpos(1,n) = measurments_queue{n}.sat_x ;
        satpos(2,n) = measurments_queue{n}.sat_y ;
        satpos(3,n) = measurments_queue{n}.sat_z ;
        obs(n) = measurments_queue{n}.prMes ;
        
        Eph(:,n) = eph2easy(measurments_queue{n}.s_eph, n ) ; % convert to easy ephemeris
    end
    
%     settings.c    = 299792458 ;    % The speed of light, [m/s]
%     settings.useTropCorr = 1 ;     % Use troposphere correction
%     addpath(gnssLib) ;
%     [pos, el, az, dop] = leastSquarePos(satpos, obs, settings) ;
%     rmpath(gnssLib) ;
    
    addpath(easyLib) ;
    [easy_pos, easy_el, easy_gdop, easy_basic_obs] = recpo_ls(obs,(1:length(obs)),measurments_queue{1}.msrTow,Eph) ;
    [phi,lambda,h] = togeod(6378137,298.257223563,easy_pos(1),easy_pos(2),easy_pos(3)) ;
    rmpath(easyLib) ;

     addpath(gnssLib) ;
     [phi, lambda, h] = cart2geo(easy_pos(1), easy_pos(2), easy_pos(3), 5 ) ;
     rmpath(gnssLib) ;
    
    ux = easy_pos(1) ;
    uy = easy_pos(2) ;
    uz = easy_pos(3) ;
%     addpath(gnssLib) ;
%     [phi, lambda, h] = cart2geo(pos(1), pos(2), pos(3), 5 ) ;
%     rmpath(gnssLib) ;
%     fprintf(fout,'<NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n\n', ...
%         easy_pos(1), easy_pos(2), easy_pos(3), phi, lambda, h ) ;


