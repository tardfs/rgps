clc , clear all ;
easyLib = getFullPath('..\\easy') ;

resultFile = '..\data\rltv_out.mat' ;

fprintf('load baseline data...\n') ;
load(resultFile) ;

addpath(easyLib) ;
ubx_master = mean(ubxEcefA(2:4,10:20),2) ;
[m_phi,m_lambda,m_h] = togeod(6378137,298.257223563,ubx_master(1),ubx_master(2),ubx_master(3)) ;
fprintf('<GEO>%f,%f,%f\n', m_phi,m_lambda,m_h ) ;
rmpath(easyLib) ;

addpath(easyLib) ;
for idx=1:size(ubxEcefA,2)
    [m_phi,m_lambda,m_h] = togeod(6378137,298.257223563,ubxEcefA(2,idx),ubxEcefA(3,idx),ubxEcefA(4,idx)) ;
    fprintf('<GEO>%f,%f,%f\n', m_phi,m_lambda,m_h ) ;
end
rmpath(easyLib) ;

% Compare ECEF
titler = ['\Delta ECEF X';'\Delta ECEF Y';'\Delta ECEF Z'] ;
for nc=2:4
    figure(nc-1) ;
    hold off
    plot((easyApos(1,:)-easyApos(1,1))/60,easyApos(nc,:)-easyApos(nc,1),'Color',[0.2 0.2 0.9],'LineWidth',2) ;
    hold on
    plot((gnssApos(1,:)-gnssApos(1,1))/60,gnssApos(nc,:)-gnssApos(nc,1),'Color',[0.9 0.2 0.2],'LineWidth',2) ;
    hold on,plot((ubxEcefA(1,20:end)-ubxEcefA(1,20))/60,ubxEcefA(nc,20:end)-ubxEcefA(nc,20),'Color',[0.2 0.6 0.2],'LineWidth',2) ;
    grid on ;
    set(gca,'FontSize',14) ;
    title(titler(nc-1,:)) ;
    xlabel('Time, min') ;
    ylabel('\Delta, m') ;
    legend('easy','gnss','ubx') ;
    %print -depsc ..\tex\delta_ecef
    export_fig(  sprintf('..\\tex\\delta_ecef_%d.png', nc-1 )) ;
end

figure(5)
hold off
plot((easyApos(1,:)-easyApos(1,1))/60,easyApos(5,:),'Color',[0.3 0.2 0.9],'LineWidth',2) ;
hold on
plot((gnssApos(1,:)-gnssApos(1,1))/60,gnssApos(5,:),'b-.','Color',[0.9 0.2 0.3],'LineWidth',2) ;
hold on,plot((ubxEcefA(1,20:end)-ubxEcefA(1,20))/60,ubxEcefA(5,20:end),'Color',[0.3 0.6 0.2],'LineWidth',2) ;
grid on ;
set(gca,'FontSize',14) ;
title('GDOP') ;
ylabel('GDOP, m') ;
xlabel('Time, min') ;
legend('easy','gnss','ubx') ;


titler = ['Baseline X';'Baseline Y';'Baseline Z'] ;
for nc=2:4
    figure(nc-2+6) ;
    hold off
    plot((coarseBaseline(1,:)-coarseBaseline(1,1))/60,-coarseBaseline(nc,:),'Color',[0.2 0.2 0.9],'LineWidth',2) ;
    hold on
    plot((prBaseline(1,:)-prBaseline(1,1))/60,prBaseline(nc,:),'Color',[0.9 0.2 0.2],'LineWidth',2) ;
    grid on ;
    set(gca,'FontSize',14) ;
    title(titler(nc-1,:)) ;
    xlabel('Time, min') ;
    ylabel('Baseline, m') ;
    legend('Coarse','PR') ;
    %print -depsc ..\tex\delta_ecef
    export_fig(  sprintf('..\\tex\\baseline_%d.png', nc-1 )) ;
end



return ;


addpath(easyLib) ;
[m_phi,m_lambda,m_h] = togeod(6378137,298.257223563,master_pos(1),master_pos(2),master_pos(3)) ;
[B_phi,B_lambda,B_h] = togeod(6378137,298.257223563,B_pos(1),B_pos(2),B_pos(3)) ;
[ubx_phi,ubx_lambda,ubx_h] = togeod(6378137,298.257223563,ubx_master(1),ubx_master(2),ubx_master(3)) ;
[ubxB_phi,ubxB_lambda,ubxB_h] = togeod(6378137,298.257223563,ubx_B(1),ubx_B(2),ubx_B(3)) ;
rmpath(easyLib) ;

fprintf('Check averaged navigation solution\n') ;
fprintf('<MASTER NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        master_pos(1), master_pos(2), master_pos(3), m_phi, m_lambda, m_h ) ;
fprintf('<ROVER  NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        B_pos(1), B_pos(2), B_pos(3), B_phi, B_lambda, B_h ) ;
fprintf('<AVG BASELINE> %fm\n', sqrt((B_pos-master_pos)'*(B_pos-master_pos))) ;

fprintf('Check UBX navigation solution\n') ;
fprintf('<MASTER NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        ubx_master(1), ubx_master(2), ubx_master(3), ubx_phi, ubx_lambda, ubx_h ) ;
fprintf('<ROVER  NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
        ubx_B(1), ubx_B(2), ubx_B(3), ubxB_phi, ubxB_lambda, ubxB_h ) ;
fprintf('<AVG BASELINE> %fm\n', sqrt((ubx_B-ubx_master)'*(ubx_B-ubx_master))) ;


fprintf('\n') ;

