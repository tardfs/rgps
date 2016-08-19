 function ea_rinex2mat()
% convert easy example RINEX data to raw matlab measurements files

easyLib = getFullPath('..\\easy') ;
addpath(easyLib) ;

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat') ; 
Eph = get_eph('eph.dat') ;

% We identify the master observation file and open it
%ofile1  = 'SITE247J.01O' ; 
ofile1  = 'SITE24~1.01O' ;
outFile = '..\data\SITE24~1.mat' ;
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1) ;
NoObs_types1 = size(Obs_types1,2)/2 ;

fid1 = fopen( ofile1, 'rt' ) ;

% At end of ofile2 we overwrite empty observations with NaN's to obtain 22 valid epochs
qend = 22 ;
measurments_queue = cell(qend,1) ;
for q = 1:qend
    
    [time1, dt1, sats1, eof1] = fepoch_0(fid1) ;

    NoSv1 = size(sats1,1) ;
    obsm = grabdata(fid1, NoSv1, NoObs_types1) ;
    
    C1_type = fobs_typ(Obs_types1,'C1') ;
    P1_type = fobs_typ(Obs_types1,'P1') ;
    
    instant_measurments = cell(NoSv1,1) ;
    for n_sat_obs = 1:NoSv1
        % matlab mat part    
        s_msr = make_measurment_data() ;
        s_msr.gnssId = 0 ;
        s_msr.svId = sats1(n_sat_obs) ;
        s_msr.msrTow = time1 ;
        s_msr.week = 0 ;
        s_msr.leapS = 0 ;
        s_msr.satTow = 0 ;
        
        nc = find( Eph(1,:)==s_msr.svId, 1 ) ;        
        s_msr.s_eph = easy2eph( Eph(:,nc) ) ;
        % gate to Denis Akos code
        s_msr.prMes = obsm( n_sat_obs, C1_type ) ;
        s_msr.cpMes = obsm( n_sat_obs, P1_type ) ;
        s_msr.sat_x = 0 ;
        s_msr.sat_y = 0 ;
        s_msr.sat_z = 0 ;
        s_msr.sat_clk_corr = 0 ;
        % put measurment into the queue
        instant_measurments{n_sat_obs} = s_msr ;
    end
    
    measurments_queue{q} = instant_measurments ;
    
end

ubxEcef = 0 ;
ubxEcefCount = 0 ;
ubxGeodetic = 0 ;
ubxGeodeticCount = 0 ;
save(outFile,'measurments_queue','ubxEcef','ubxEcefCount', 'ubxGeodetic', 'ubxGeodeticCount') ;

fclose( fid1 ) ;
rmpath(easyLib) ;



function s_msr = make_measurment_data()
    s_msr.gnssId = -1 ;
    s_msr.svId = -1 ;
    s_msr.msrTow = 0 ;
    s_msr.week = 0 ;
    s_msr.leapS = 0 ;
    s_msr.prMes = 0 ;
    s_msr.cpMes = 0 ;
    s_msr.satTow = 0 ;
    s_msr.sat_x = 0 ;
    s_msr.sat_y = 0 ;
    s_msr.sat_z = 0 ;
    s_msr.sat_clk_corr = 0 ;
    s_msr.s_eph = {} ;
