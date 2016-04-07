function [Az, El, sMap, hlsMap, sat_list] = nlib_plot_sats(settings,measurmentsA)

%settings.easyLib = getFullPath('..\\easy') ;
%settings.gnssLib = getFullPath('..\\softgnss') ;

%settings.recv1File = '..\data\RS_matv_1400mm_680mm_01.mat ' ;

%fprintf(repmat('\b',1,160)) ; fprintf('load receiver A data...') ;
%load(settings.recv1File) ;

%measurmentsA = measurments_queue ;

%fprintf(repmat('\b',1,160)) ;
%fprintf('\n') ;


addpath(settings.easyLib) ;

N = length(measurmentsA) ;
easy_data = zeros(5,N) ;
sMap = zeros(32,N) ;
Az = zeros(32,N) ;
El = zeros(32,N) ;
for idx1=1:N
    recv1_msr_set = measurmentsA{idx1} ;
    [ux,uy,uz,gdop, sat_map ] = easy_pvt_solver( recv1_msr_set) ;
    easy_data(:,idx1) = [recv1_msr_set{1}.msrTow; ux;uy;uz; gdop(1)] ;
    satIds = sat_map(:,1) ;
    sMap(satIds,idx1) = 1 ;
    Az(satIds,idx1) = sat_map(:,2) ;
    El(satIds,idx1) = sat_map(:,3) ;
    
%    addpath(settings.gnssLib) ;
%    skyPlot(sat_map(:,2),sat_map(:,3),satIds) ;
%    rmpath(settings.gnssLib) ;
    
end

rmpath(settings.easyLib) ;

[hlsMap, sat_list] = nlib_highlight_smap(settings,sMap, El) ;


nlib_skyplot( Az, El, sMap, sat_list ) ;

% hold off,
% for s=1:32
%     satAz = Az(s,sMap(s,:)~=0)/180*pi ;
%     satEl = El(s,sMap(s,:)~=0)/180*pi ;
%     if ~isempty(satAz)
%         h= polar( satAz, cos(satEl) ) ;
%         set(h,'LineWidth',2,'Color',[0.2 0.4 0.9]) ;
%         hold on ,
%         h = polar( satAz(end), cos(satEl(end)), 'rs' ) ;
%         set(h,'MarkerSize',14) ;
%         [xt,yt] = pol2cart(satAz(end), cos(satEl(end))) ;
%         text(xt,yt,sprintf('%1d',s),'HorizontalAlignment','Center','FontSize',12) ;
%     end
% end
% set(gca,'FontSize',14)
% title(sprintf('Карта движения спутников %s', settings.recv1File )) ;

function [ux,uy,uz,gdop,sat_map] = easy_pvt_solver(measurments_queue)

numSat = length(measurments_queue) ;
obs = zeros(numSat,1) ;

Eph = zeros(21,length(measurments_queue)) ;
for n=1:length(measurments_queue)
    obs(n) = measurments_queue{n}.prMes ;        
    Eph(:,n) = eph2easy(measurments_queue{n}.s_eph, n ) ; % convert to easy ephemeris
end

[easy_pos, el, gdop, basic_obs] = recpo_ls(obs,(1:length(obs)),measurments_queue{1}.msrTow,Eph) ;

satPos = basic_obs(:,1:3) ;

sat_map = zeros(numSat,4) ;
for n=1:numSat
    [Az, El, D] = topocent(easy_pos(1:3),satPos(n,:).'-easy_pos(1:3)) ;
    sat_map(n, 1) = measurments_queue{n}.svId ;
    sat_map(n, 2) = Az ;
    sat_map(n, 3) = El ;
end

%[phi,lambda,h] = togeod(6378137,298.257223563,easy_pos(1),easy_pos(2),easy_pos(3)) ;

ux = easy_pos(1) ;
uy = easy_pos(2) ;
uz = easy_pos(3) ;

