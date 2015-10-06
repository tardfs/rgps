function ubxdump()

% TOW - epoch counter, epoch length=6sec
% all TOW fields in sec!!!

gnssLib = getFullPath('..\\softgnss') ;
addpath(gnssLib) ;

sattelites = {} ;
measurments_count = 0 ;
measurments_queue = {} ;
ubxEcefCount = 0 ;
ubxEcef = zeros(5,1000) ;
ubxGeodeticCount = 0 ;
ubxGeodetic = zeros(7,1000) ;
clc ;
% clear all ;
outFile = '..\data\RS_twin_30mm_01.mat' ;
f = fopen('..\data\RS_twin_30mm_01.bin','r') ; % ..\data\rs_004.bin
fout = fopen('ubx_out.txt','w+t') ;
n_show=0 ;
n_eph_t = 0 ;
n_measurment_t = 0 ;
if f~=-1
    packet_count = 0 ;
    crc_ok = 0 ;
    longest_packet = 0 ;
    longest_crc_ok = 0 ;
    while(~feof(f))
        %if ftell(f)>2000000
        %    break ;
        %end
        if ftell(f)>n_show
            fprintf(repmat('\b',1,160)) ; fprintf( '<FILE POS>%07d <RAW MEASURMENT_t:>%10.3f <RAW MEASURMENTS COLLECTED:>%4d <UBX_ECEF:>%4d <UBX_GEODETIC:>%4d', n_show, n_measurment_t, measurments_count, ubxEcefCount, ubxGeodeticCount ) ;
            n_show = n_show+10000 ;
        end
        % Read sync char 1        
        sync_byte1 = fread(f,1,'uint8') ;
        if sync_byte1==hex2dec('B5')
            sync_byte2 = fread(f,1,'uint8') ;
            if sync_byte2==hex2dec('62')
                % looks like UBX header
                crc_data_start = ftell(f) ;
                msg_class = fread(f,1,'uint8') ;
                msg_id = fread(f,1,'uint8') ;
                payload_len = fread(f,1,'uint16',0,'ieee-le') ;

                save_pos = ftell(f) ;
                poff = save_pos ;

                packet_count = packet_count + 1 ;
                
                % compute check sum
                fseek(f,crc_data_start,'bof') ;
                crc_data_length = 4 + payload_len ;
                CK_A = 0 ;
                CK_B = 0 ;
                for crc_idx = 1:crc_data_length
                    v = fread(f,1,'uint8') ;
                    CK_A = mod(CK_A + v,256) ;
                    CK_B = mod(CK_B + CK_A,256) ;
                end
                %fseek(f,save_pos,'bof') ;
                %poff = ftell(f) ;
                field_off = ftell(f)-poff ;
                CK_1 = fread(f,1,'uint8') ;
                CK_2 = fread(f,1,'uint8') ;
                if CK_1==CK_A && CK_2==CK_B
                    crc_ok = crc_ok + 1 ;
                    if longest_crc_ok<crc_data_length
                        longest_crc_ok = crc_data_length ;
                    end
                end
                if longest_packet<crc_data_length
                    longest_packet = crc_data_length ;
                end
                
                fprintf(fout,'PKT:<%d> OFF:<%d> UBX MSG_CLASS:<%1d>\n', packet_count, crc_data_start-2, msg_class) ;
                switch msg_class
                    case 1
                        fprintf(fout,'\t<NAV> Navigation Results: Position, Speed, Time, Acceleration, Heading, DOP, SVs used\n') ;
                    case 2
                        fprintf(fout,'\t<RXM> Receiver Manager Messages: Satellite Status, RTC Status\n') ;
                    case 4
                        fprintf(fout,'\t<INF> Information Messages: Printf-Style Messages, with IDs such as Error, Warning, Notice\n') ;
                    case 5
                        fprintf(fout,'\t<ACK> Ack/Nack Messages: as replies to CFG Input Messages\n') ;
                    case 6
                        fprintf(fout,'\t<CFG> Configuration Input Messages: Set Dynamic Model, Set DOP Mask, Set Baud Rate, etc.\n') ;
                    case 9
                        fprintf(fout,'\t<UPD> Firmware Update Messages: Memory/Flash erase/write, Reboot, Flash identification, etc..\n') ;
                    case 10
                        fprintf(fout,'\t<MON> Monitoring Messages: Comunication Status, CPU Load, Stack Usage, Task Status\n') ;
                    case 11
                        fprintf(fout,'\t<AID> AssistNow Aiding Messages: Ephemeris, Almanac, other A-GPS data input\n') ;
                    case 12
                        fprintf(fout,'\t<TIM> Timing Messages: Time Pulse Output, Timemark Results\n') ;
                    case 19
                        fprintf(fout,'\t<MGA> Multi-GNSS Assistance: Assistance data for various GNSS\n') ;
                    case 33
                        fprintf(fout,'\t<LOG> Logging Messages: Log creation, deletion, info and retrieval\n') ;
                    otherwise
                        fprintf(fout,'\t<!!!> Unknown msg_class\n') ;
                end
                fprintf(fout,'UBX MSG_ID:<%1d>\n', msg_id ) ;
                fprintf(fout,'UBX PAYLOAD LENGTH:%d\n', payload_len ) ;
                
                if (CK_1==CK_A && CK_2==CK_B) %|| (msg_class==2 && msg_id==19)  % CRC is Ok - proceed with parsing
                    fseek(f,save_pos,'bof') ;
                    if msg_class==1
                        switch msg_id
                            case 1
                                fprintf(fout,'\t<NAV-POSECEF> Position Solution in ECEF\n') ;
                                field_off = ftell(f)-poff ;
                                iTOW = fread(f,1,'uint32') ;
                                fprintf(fout,'\t\t[%3d]<iTOW> %d ms\n', field_off, iTOW ) ;
                                field_off = ftell(f)-poff ;
                                ecefX = fread(f,1,'int32') ;
                                fprintf(fout,'\t\t[%3d]<ecefX> cm %d\n', field_off, ecefX ) ;
                                field_off = ftell(f)-poff ;
                                ecefY = fread(f,1,'int32') ;
                                fprintf(fout,'\t\t[%3d]<ecefY> cm %d\n', field_off, ecefY ) ;
                                field_off = ftell(f)-poff ;
                                ecefZ = fread(f,1,'int32') ;
                                fprintf(fout,'\t\t[%3d]<ecefZ> cm %d\n', field_off, ecefZ ) ;
                                field_off = ftell(f)-poff ;
                                pAcc = fread(f,1,'uint32') ;
                                fprintf(fout,'\t\t[%3d]<pAcc> cm %d\n', field_off, pAcc ) ;
                                ubxEcefCount = ubxEcefCount+1 ;
                                ubxEcef(:,ubxEcefCount) = [ iTOW*1e-3; ... % convert to sec
                                                            ecefX/100; ... % convert to m
                                                            ecefY/100; ...
                                                            ecefZ/100; ...
                                                            pAcc/100 ] ;
                            case 2                                
                                fprintf(fout,'\t<UBX-NAV-POSLLH> Geodetic Position Solution\n') ;
                                ubxGeodeticCount = ubxGeodeticCount + 1 ;
                                field_off = ftell(f)-poff ;
                                iTOW = fread(f,1,'uint32') ;
                                fprintf(fout,'\t\t[%3d]<iTOW> %d ms\n', field_off, iTOW ) ;
                                field_off = ftell(f)-poff ;
                                lon = fread(f,1,'int32')*1e-7 ;
                                fprintf(fout,'\t\t[%3d]<long> %5.3f\n', field_off, lon ) ;
                                field_off = ftell(f)-poff ;
                                lat = fread(f,1,'int32')*1e-7 ;
                                fprintf(fout,'\t\t[%3d]<lat><long>  %7.5f,%7.5f\n', field_off, lat, lon ) ;
                                field_off = ftell(f)-poff ;
                                height = fread(f,1,'int32')/100 ;
                                fprintf(fout,'\t\t[%3d]<height>  %f m\n', field_off, height ) ;
                                field_off = ftell(f)-poff ;
                                hMSL = fread(f,1,'int32')/100 ;
                                fprintf(fout,'\t\t[%3d]<hMSL>  %f m\n', field_off, hMSL ) ;
                                field_off = ftell(f)-poff ;
                                hAcc = fread(f,1,'uint32')/100 ;
                                fprintf(fout,'\t\t[%3d]<hAcc>  %f m\n', field_off, hAcc ) ;
                                field_off = ftell(f)-poff ;
                                vAcc = fread(f,1,'uint32')/100 ;
                                fprintf(fout,'\t\t[%3d]<vAcc>  %f m\n', field_off, vAcc ) ;
                                
                                ubxGeodetic(:,ubxGeodeticCount) = [ iTOW; ...
                                                                    lat; lon ; ...
                                                                    height; ...
                                                                    hMSL; ...
                                                                    hAcc; ...
                                                                    vAcc ] ;
                        end
                    elseif msg_class==2
                        switch msg_id
                            case 19
                                fprintf(fout,'\t<RXM-SFRBX> Raw Subframe Data\n') ;
                                field_off = ftell(f)-poff ;
                                gnssId = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<gnssId> %d\n', field_off, gnssId ) ;
                                field_off = ftell(f)-poff ;
                                svId = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<svId> %d\n', field_off, svId ) ;
                                field_off = ftell(f)-poff ;
                                reserved1 = fread(f,1,'uint8') ;
                                if reserved1==0
                                    fprintf(fout,'\t\t[%3d]<reserved1> %d\n', field_off, reserved1 ) ;
                                else
                                    fprintf(fout,'\t\t[%3d]<reserved1> %d - Error. Should be zero value.\n', field_off, reserved1 ) ;
                                end
                                field_off = ftell(f)-poff ;
                                freqId = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<freqId> %d\n', field_off, freqId ) ;
                                field_off = ftell(f)-poff ;
                                numWords = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<numWords> %d\n', field_off, numWords ) ;
                                field_off = ftell(f)-poff ;
                                reserved2 = fread(f,1,'uint8') ;
                                %if reserved2==0
                                    fprintf(fout,'\t\t[%3d]<reserved2> %d\n', field_off, reserved2 ) ;
                                %else
                                %    fprintf(fout,'\t\t[%3d]<reserved2> %d - Error. Should be zero value.\n', field_off, reserved2 ) ;
                                %end
                                field_off = ftell(f)-poff ;
                                version = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<version> %d\n', field_off, version ) ;
                                field_off = ftell(f)-poff ;
                                reserved3 = fread(f,1,'uint8') ;
                                %if reserved3==0
                                    fprintf(fout,'\t\t[%3d]<reserved3> %d\n', field_off, reserved3 ) ;
                                %else
                                %    fprintf(fout,'\t\t[%3d]<reserved3> %d - Error. Should be zero value.\n', field_off, reserved3 ) ;
                                %end
                                if gnssId==0
                                    fprintf(fout,'\t\tgnssId <GPS>:\n') ;

                                     s_i = get_sattelite(sattelites,gnssId, svId) ;
                                     if s_i==-1
                                         s_data = make_sattelite_data() ;
                                         s_data.gnssId = 0 ;
                                         s_data.svId = svId ;
                                         sattelites{end+1} = s_data ;
                                         s_i = length(sattelites) ;
                                     end

                                    % check parity bits
                                    raw_subframe = zeros(numWords,1) ;
                                    bin_subframe = repmat('0',numWords,30) ;
                                    hcrc_ok = zeros(numWords,1) ;
                                    %D29 = 0 ;
                                    %D30 = 0 ;
                                    for nw=1:numWords
                                        field_off = ftell(f)-poff ;
                                        raw_subframe(nw) = fread(f,1,'uint32') ;
                                        msg_word = dec2bin(raw_subframe(nw),32) ;
                                        %msg_word = msg_word(end:-1:1) ;
                                        bin_subframe(nw,:) = msg_word(3:end) ;

                                        tmp = bin2int(msg_word(1:2)) ;
                                        D29 = tmp(1) ;
                                        D30 = tmp(2) ;

                                        fprintf(fout,'\t\t\t[%3d]<word %2d> %08X %s %s D29*:%1d D30*:%1d', field_off, nw, raw_subframe(nw), msg_word(1:2), bin_subframe(nw,:), D29, D30 ) ;

                                        [data_ok,D29,D30,pbits] = check_hamming(bin_subframe(nw,:), D29, D30) ;
                                        if data_ok 
                                            fprintf(fout, ' %s Parity check ok\n', pbits ) ;
                                            hcrc_ok(nw) = 1 ;
                                        else
                                            fprintf(fout, ' %s Parity check error\n', pbits ) ;                                        
                                        end
                                    end

                                    SFID = -1 ;
                                    if hcrc_ok(1)==1
                                        on_tlm(fout, bin_subframe(1,:)) ;
                                    end
                                    if hcrc_ok(2)==1
                                        SFID = on_how(fout, bin_subframe(2,:)) ;
                                    end
                                    switch SFID
                                        case 1
                                            if hcrc_ok(3)==1
                                                sattelites{s_i} = on_sf1(fout, bin_subframe,sattelites{s_i}) ;
                                            end                                        
                                        case 2
                                            if hcrc_ok(3)==1
                                                sattelites{s_i} = on_sf2(fout, bin_subframe,sattelites{s_i}) ;
                                            end                                        
                                        case 3
                                            if hcrc_ok(3)==1
                                                sattelites{s_i} = on_sf3(fout, bin_subframe, sattelites{s_i}) ;
                                            end                                        
                                    end

                                    % TLM
    %                                 nw = 0 ;
    %                                 field_off = ftell(f)-poff ;
    %                                 gpsWord = fread(f,1,'uint32') ;
    %                                 fprintf(fout,'\t\t\t[%3d]<word %2d> %08X %s\n', field_off, nw, gpsWord, dec2bin(gpsWord,32) ) ;
    %                                 tlmBits = dec2bin(gpsWord,32) ; 
    %                                 %tlmBits= tlmBits(3:end) ;
    %                                 on_tlm(fout,tlmBits(3:end)) ;
    %                                 [data_ok,D29,D30] = check_hamming(tlmBits(3:end),0,0) ;
    %                                 if data_ok
    %                                     fprintf(fout,'\t\t\tTLM parity check Ok\n') ;                                    
    %                                 end
    %                                 if strcmpi(tlmBits(1:8),'10001011')==1
    %                                     fprintf(fout,'\t\t\t\tTLM preamble: %s.\n', tlmBits(1:8) ) ;                                    
    %                                     % HOW
    %                                     nw = 1 ;
    %                                     field_off = ftell(f)-poff ;
    %                                     gpsWord = fread(f,1,'uint32') ;
    %                                     fprintf(fout,'\t\t\t[%3d]<word %2d> %08X %s\n', field_off, nw, gpsWord, dec2bin(gpsWord,32) ) ;
    %                                     HowBits = dec2bin(gpsWord,32) ; 
    %                                     HowBits= HowBits(3:end) ;
    %                                     TOW = bin2dec(HowBits(1:17)) ;
    %                                     SFID = bin2dec(HowBits(20:22)) ;
    %                                     fprintf(fout, '\t\t\t\t<TOW> %d\n', TOW ) ;
    %                                     fprintf(fout, '\t\t\t\t<SFID> %d\n', SFID ) ;                                    
    %                                 else                                    
    %                                     fprintf(fout,'\t\t\t\t[Error]: TLM preamble missed\n') ;
    %                                 end
                                end
                            case 21
                                fprintf(fout,'\t<RXM-RAWX> Multi-GNSS Raw Measurement Data\n') ;
                                instant_measurments = {} ;
                                
                                field_off = ftell(f)-poff ;
                                rcvTow = fread(f,1,'double') ;
                                fprintf(fout,'\t\t[%3d]<rcvTow> %f seconds\n', field_off, rcvTow ) ;
                                n_measurment_t = rcvTow ;
                                field_off = ftell(f)-poff ;
                                gpsWeek = fread(f,1,'uint16') ;
                                fprintf(fout,'\t\t[%3d]<week> %d\n', field_off, gpsWeek ) ;
                                field_off = ftell(f)-poff ;
                                leapS = fread(f,1,'int8') ;
                                fprintf(fout,'\t\t[%3d]<leapS> %d\n', field_off, leapS ) ;
                                field_off = ftell(f)-poff ;
                                numMeas = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<numMeas> %d\n', field_off, numMeas ) ;
                                field_off = ftell(f)-poff ;
                                recStat = fread(f,1,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<recStat> %s\n', field_off, dec2bin(recStat,8) ) ;
                                field_off = ftell(f)-poff ;
                                reserved1 = fread(f,3,'uint8') ;
                                fprintf(fout,'\t\t[%3d]<reserved1> %d %d %d\n', field_off, reserved1(1),reserved1(2),reserved1(3) ) ;                            
                                for nmeasure=1:numMeas
                                    fprintf(fout,'\t\t\t<Measurment> %d\n', nmeasure ) ;
                                    field_off = ftell(f)-poff ;
                                    prMes = fread(f,1,'double') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<prMes/pseudorange> %f m\n', field_off, prMes ) ;
                                    field_off = ftell(f)-poff ;
                                    cpMes = fread(f,1,'double') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<cpMes/carrier phase> %f cycles\n', field_off, cpMes ) ;
                                    field_off = ftell(f)-poff ;
                                    doMes = fread(f,1,'single') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<doMes/Doppler> %f Hz\n', field_off, doMes ) ;
                                    field_off = ftell(f)-poff ;
                                    gnssId = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<gnssId> %d\n', field_off, gnssId ) ;
                                    field_off = ftell(f)-poff ;
                                    svId = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<svId> %d\n', field_off, svId ) ;
                                    field_off = ftell(f)-poff ;
                                    reserved2 = fread(f,1,'uint8') ;
                                    if reserved2==0
                                        fprintf(fout,'\t\t\t\t[%3d]<reserved2> %d\n', field_off, reserved2 ) ;
                                    else
                                        fprintf(fout,'\t\t\t\t[%3d]<reserved2> %d - Error. Should be zero value.\n', field_off, reserved2 ) ;
                                    end
                                    field_off = ftell(f)-poff ;
                                    freqId = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<freqId> %d\n', field_off, freqId ) ;
                                    field_off = ftell(f)-poff ;
                                    locktime = fread(f,1,'uint16') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<PLL locktime> %d\n', field_off, locktime ) ;
                                    field_off = ftell(f)-poff ;
                                    cno = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<cno> %d\n', field_off, cno ) ;
                                    field_off = ftell(f)-poff ;
                                    prStdev = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<prStdev> %s\n', field_off, dec2bin(prStdev,8) ) ;
                                    field_off = ftell(f)-poff ;
                                    cpStdev = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<cpStdev> %s\n', field_off, dec2bin(cpStdev,8) ) ;
                                    field_off = ftell(f)-poff ;
                                    doStdev = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<doStdev> %s\n', field_off, dec2bin(doStdev,8) ) ;
                                    field_off = ftell(f)-poff ;
                                    trkStat = fread(f,1,'uint8') ;
                                    fprintf(fout,'\t\t\t\t[%3d]<trkStat> %s\n', field_off, dec2bin(trkStat,8) ) ;
                                    reserved3 = fread(f,1,'uint8') ;
                                    
                                    % try to find appropriate ephemeris
                                    s_index = get_sattelite(sattelites, gnssId, svId) ;
                                    if s_index>0
                                        e_index = get_ephemeris_ontime(sattelites{s_index}.ephemerises,rcvTow) ;
                                        if e_index>0
                                            s_eph = sattelites{s_index}.ephemerises{e_index} ;
                                            s_msr = make_measurment_data() ;
                                            s_msr.gnssId = gnssId ;
                                            s_msr.svId = svId ;
                                            s_msr.msrTow = rcvTow ;
                                            s_msr.week = gpsWeek ;
                                            s_msr.leapS = leapS ;
                                            s_msr.satTow = s_eph.toc ;
                                            s_msr.s_eph = s_eph ;
                                            % gate to Denis Akos code
                                            akos_eph = make_akos_eph(s_eph) ;
                                            [satPositions, satClkCorr] = satpos(rcvTow, 1, akos_eph) ;
                                            settings.c    = 299792458 ;    % The speed of light, [m/s]
                                            s_msr.prMes = prMes ; %+ satClkCorr * settings.c ;
                                            s_msr.sat_x = satPositions(1) ;
                                            s_msr.sat_y = satPositions(2) ;
                                            s_msr.sat_z = satPositions(3) ;
                                            s_msr.sat_clk_corr = satClkCorr ;
                                            fprintf(fout, '\t\t\t\t<SATTELITE POSITION><SVID>%d <POS>(%f,%f,%f) <CLK_CORR>%f\n',...
                                                svId, s_msr.sat_x,s_msr.sat_y,s_msr.sat_z, s_msr.sat_clk_corr ) ;
                                            sat_distance = sqrt(s_msr.sat_x*s_msr.sat_x+s_msr.sat_y*s_msr.sat_y+s_msr.sat_z*s_msr.sat_z) ;
                                            fprintf(fout, '\t\t\t\t<SAT DISTANCE>%f, <PSEUDORANGE>%f\n', sat_distance, s_msr.prMes ) ;
                                            % put measurment into the queue
                                            instant_measurments{end+1} = s_msr ;
                                        end
                                    end
                                end
                                if length(instant_measurments)>=4                                    
                                    pvt_solver(fout,instant_measurments,4) ; % standart PVT solver
                                    measurments_count = measurments_count + 1 ;
                                    measurments_queue{end+1} = instant_measurments ;
                                end
                        end
                    end
                    % always seek to crc pos
                    fseek(f,crc_data_start+crc_data_length,'bof') ;
                    
                    % double check
                    field_off = ftell(f)-poff ;
                    CK_1 = fread(f,1,'uint8') ;
                    CK_2 = fread(f,1,'uint8') ;
                    fprintf(fout,'\t[%3d]<CRC> %02x%02x, %02x%02x <Ok>\n', field_off, CK_1, CK_2, CK_A, CK_B ) ;
                else
                    if msg_class==2
                        switch msg_id
                            case 19
                                fprintf(fout,'\t<RXM-SFRBX> Raw Subframe Data\n') ;
                            case 21
                                fprintf(fout,'\t<RXM-RAWX> Multi-GNSS Raw Measurement Data\n') ;                                
                        end
                    end
                    % bad crc packet
                    fprintf(fout,'\t[%3d]<CRC> %02x%02x, %02x%02x <Failed>\n', field_off, CK_1, CK_2, CK_A, CK_B ) ;
                    % bad sequence protection
                    %fseek(f,crc_data_start,'bof') ;
                end                
                fprintf(fout,'-/\n\n') ;
            end
        end
    end
    fprintf(fout, 'TOTAL PACKETS FOUND: %d\n', packet_count ) ;
    fprintf(fout, 'CRC Ok: %d packets\n', crc_ok ) ;
    fprintf(fout, 'Maximum packet size: %d bytes\n', longest_packet ) ;
    fprintf(fout, 'Maximum packet size (CRC Ok): %d bytes\n', longest_crc_ok ) ;
    fclose(f) ;
    dump_sattelites(fout,sattelites) ;
    fprintf('\n') ;
else
    fprintf(fout,'ubxdump error: Can''t open input file.\n') ;
    if fout~=1
        fprintf('ubxdump error: Can''t open input file.\n') ;
    end
end
ubxEcef = ubxEcef(:,1:ubxEcefCount) ;
ubxGeodetic = ubxGeodetic(:,1:ubxGeodeticCount) ;
save(outFile,'measurments_queue','ubxEcef','ubxEcefCount', 'ubxGeodetic', 'ubxGeodeticCount') ;
fclose(fout) ;
rmpath(gnssLib) ;

function v = bin2int(binstr)
N = length(binstr) ;
v = zeros(N,1) ;
for n=1:N
    if binstr(n)=='1'
        v(n) = 1 ;
    end
end

function v = int2bin(intvect)
N = length(intvect) ;
v = repmat('0',1,N) ;
v(logical(intvect))='1' ;

function [data_ok,D29,D30, parity_bits]=check_hamming(data_bits, D29P, D30P)
data = bin2int(data_bits) ;
data(1:24) = mod(data(1:24) + D30P,2) ;
D25 = mod(D29P+sum(data([1:3,5,6,10:14,17,18,20,23])),2) ;
D26 = mod(D30P+sum(data([2:4,6,7,11:15,18,19,21,24])),2) ;
D27 = mod(D29P+sum(data([1,3:5,7,8,12:16,19,20,22])),2) ;
D28 = mod(D30P+sum(data([2,4:6,8,9,13:17,20,21,23])),2) ;
D29 = mod(D30P+sum(data([1,3,5:7,9,10,14:18,21,22,24])),2) ;
D30 = mod(D29P+sum(data([3,5,6,8:11,13,15,19,22:24])),2) ;
parity_bits = int2bin([D25,D26,D27,D28,D29,D30]) ;

if D25==data(25) && D26==data(26) && D27==data(27) && D28==data(28) && ...
        D29==data(29) && D30==data(30)
    data_ok = 1 ;
else
    data_ok = 0 ;
end

function status = on_tlm(fout,tlm_bits)
    status = 1 ;
    % check for preamble
    if strcmpi(tlm_bits(1:8),'10001011')==1
        fprintf(fout,'\t\t\t\tTLM preamble: %s.\n', tlm_bits(1:8) ) ;                                            
    else
        fprintf(fout,'\t\t\t\t[Error]: TLM preamble missed\n') ;
        status = 0 ;
        return ;
    end
    

function [SFID] = on_how(fout,inp_bits)
    TOW = bin2dec(inp_bits(1:17)) ;
    SFID = bin2dec(inp_bits(20:22)) ;
    fprintf(fout, '\t\t\t\t<TOW> %d, %f seconds\n', TOW, TOW*6 ) ; % See Tsui, page 80
    fprintf(fout, '\t\t\t\t<SFID> %d\n', SFID ) ;

function status = on_word3(fout,inp_bits)
    status = 1 ;
    URA = bin2dec(inp_bits(13:16)) ;
    fprintf(fout, '\t\t\t\t<URA Index> %d, pseudorange precision:', URA ) ;    
    switch URA
        case 0
            fprintf(fout, '0< URA <= 2.4m\n' ) ;
        case 1
            fprintf(fout, '2.4m< URA <= 3.4m\n' ) ;
        case 2
            fprintf(fout, '3.4m< URA <= 4.85m\n' ) ;
        case 3
            fprintf(fout, '4.85m< URA <= 6.85m\n' ) ;
        case 4
            fprintf(fout, '6.85m< URA <= 9.65m\n' ) ;
        case 5
            fprintf(fout, '9.65m< URA <= 13.65m\n' ) ;
        case 6
            fprintf(fout, '13.65m< URA <= 24.0m\n' ) ;
        case 7
            fprintf(fout, '24.0m< URA <= 48.0m\n' ) ;
        case 8
            fprintf(fout, '48.0m< URA <= 96.0m\n' ) ;
        case 9
            fprintf(fout, '96.0m< URA <= 192.0m\n' ) ;
        case 10
            fprintf(fout, '192.0m< URA <= 384.0m\n' ) ;
        case 11
            fprintf(fout, '384m< URA <= 768m\n' ) ;
        case 12
            fprintf(fout, '768m< URA <= 1536m\n' ) ;
        case 13
            fprintf(fout, '1536m< URA <= 3072m\n' ) ;
        case 14
            fprintf(fout, '3072m< URA <= 6144m\n' ) ;
        case 15
            fprintf(fout, '6144m< URA <= ...\n' ) ;
    end
            
    fprintf(fout, '\t\t\t\t<SV HEALTH:> %s\n', inp_bits(17:22) ) ;
    fprintf(fout, '\t\t\t\t<bits 8,7 of IODC:> %s\n', inp_bits(23:24)) ;
    
function sattelite=on_sf1(fout,inp_bits,sattelite)
    WN = bin2dec(inp_bits(3,1:10)) ;
    URA = bin2dec(inp_bits(3,13:16)) ;
    SV_HEALTH = inp_bits(3,17:22) ;
    IODC = [inp_bits(3,23:24),inp_bits(8,1:8)] ;    % IODC
    Tgd = twosComp2dec(inp_bits(7,17:24))*2^-31 ;   % T_GD
    toc = bin2dec(inp_bits(8,9:24))*2^4 ;           % t_oc
    af2 = twosComp2dec(inp_bits(9,1:8))*2^-55 ;     % a_f2
    af1 = twosComp2dec(inp_bits(9,9:24))*2^-43 ;    % a_f1
    af0 = twosComp2dec(inp_bits(10,1:22))*2^-31 ;   % a_f0
    
    sattelite.WN = WN ;
    sattelite.URA = URA ;
    sattelite.SV_HEALTH = SV_HEALTH ;
    e_i = get_ephemeris(sattelite.ephemerises,IODC(3:end)) ;
    if e_i == -1
        e_data = make_ephemeris_data() ;
        e_data.IODC = IODC ;
        sattelite.ephemerises{end+1} = e_data ;
        e_i = length(sattelite.ephemerises) ;
    end
    sattelite.ephemerises{e_i}.IODC = IODC ; % 2 MSB bits 
    sattelite.ephemerises{e_i}.Tgd = Tgd ;
    sattelite.ephemerises{e_i}.toc = toc ;
    sattelite.ephemerises{e_i}.af2 = af2 ;
    sattelite.ephemerises{e_i}.af1 = af1 ;
    sattelite.ephemerises{e_i}.af0 = af0 ;
    sattelite.ephemerises{e_i}.sf_flags(1) = 1 ;

    fprintf(fout, '\t\t\t\tSubframe 1 parising...\n' ) ;
    fprintf(fout, '\t\t\t\t<IODC>\t%s\t%d\t /* Issue of data */\n', IODC, bin2dec(IODC) ) ;
    fprintf(fout, '\t\t\t\t<WN>\t%d\t/* Week number */\n', WN ) ;
    fprintf(fout, '\t\t\t\t<URA>\t%d\t/* User range accuracy */\n', URA ) ;
    fprintf(fout, '\t\t\t\t<SV_HEALTH>\t%s\t/* Sattelite Health */\n', SV_HEALTH ) ;
    fprintf(fout, '\t\t\t\t<Tgd>\t%f seconds \t/* Sattelite group delay differential */\n', Tgd ) ;    
    fprintf(fout, '\t\t\t\t<toc>\t%f seconds \t/* Sattelite clock correction */\n', toc ) ;    
    fprintf(fout, '\t\t\t\t<af2>\t%f seconds \t/* Sattelite clock correction */\n', af2 ) ;    
    fprintf(fout, '\t\t\t\t<af1>\t%f seconds \t/* Sattelite clock correction */\n', af1 ) ;    
    fprintf(fout, '\t\t\t\t<af0>\t%f seconds \t/* Sattelite clock correction */\n', af0 ) ;    

function sattelite=on_sf2(fout,inp_bits,sattelite)
    IODE = inp_bits(3,1:8) ;
    C_rs = twosComp2dec(inp_bits(3,9:24))*2^-5 ;                  % C_rs
    Delta_n = twosComp2dec(inp_bits(4,1:16))*(2^-43)*pi ;         % deltan 
    M_0 = twosComp2dec([inp_bits(4,17:24),inp_bits(5,1:24)])*(2^-31)*pi ; % M_0
    C_uc = twosComp2dec(inp_bits(6,1:16))*(2^-29) ;               % C_uc
    e_s = bin2dec([inp_bits(6,17:24),inp_bits(7,1:24)])*2^-33 ;   % e
    C_us = twosComp2dec(inp_bits(8,1:16))*2^-29 ;                 % C_us
    sqrt_as = bin2dec([inp_bits(8,17:24),inp_bits(9,1:24)])*2^-19 ; % sqrtA
    t_oe = bin2dec(inp_bits(10,1:16))*2^4 ;                       % t_oe

    e_i = get_ephemeris(sattelite.ephemerises,IODE) ;
    if e_i == -1
        e_data = make_ephemeris_data() ;
        e_data.IODC(3:end) = IODE ;
        sattelite.ephemerises{end+1} = e_data ;
        e_i = length(sattelite.ephemerises) ;
    end
    sattelite.ephemerises{e_i}.C_rs = C_rs ;
    sattelite.ephemerises{e_i}.Delta_n = Delta_n ;
    sattelite.ephemerises{e_i}.M_0 = M_0 ;
    sattelite.ephemerises{e_i}.C_uc = C_uc ;
    sattelite.ephemerises{e_i}.e_s = e_s ;
    sattelite.ephemerises{e_i}.C_us = C_us ;
    sattelite.ephemerises{e_i}.sqrt_as = sqrt_as ;
    sattelite.ephemerises{e_i}.t_oe = t_oe ;
    sattelite.ephemerises{e_i}.sf_flags(2) = 1 ;
    
    fprintf(fout, '\t\t\t\tSubframe 2 parising...\n' ) ;
    fprintf(fout, '\t\t\t\t<IODE>\t%s\t%d\n', IODE, bin2dec(IODE) ) ;
    fprintf(fout, '\t\t\t\t<C_rs>\t%f m \t/* Amplitude of sine harmonic correction terms to the orbit radius */\n', C_rs ) ;
    fprintf(fout, '\t\t\t\t<d_n>\t%f semicircle\t/sec /* Mean motion difference from computed value */\n', Delta_n ) ;
    fprintf(fout, '\t\t\t\t<M_0>\t%f semicircle\t/* Mean anomaly at reference time */\n', M_0 ) ;
    fprintf(fout, '\t\t\t\t<C_uc>\t%f radians\t/* Amplitude of the cosine harmonic correction term to the argument of latitude ... */\n', C_uc ) ;
    fprintf(fout, '\t\t\t\t<e_s>\t%f \t/* Orbital eccentricity */\n', e_s ) ;
    fprintf(fout, '\t\t\t\t<C_us>\t%f radians\t/* Amplitude of the sine harmonic correction term to the argument of latitude */\n', C_us ) ;
    fprintf(fout, '\t\t\t\t<Sqrt(a_s)>\t%f sqrt(m)\t/* Square root of the semimajor axis */\n', sqrt_as ) ;
    fprintf(fout, '\t\t\t\t<t_oe>\t%f seconds\t/* reference time of ephemeris */\n', t_oe ) ;

function sattelite=on_sf3(fout,inp_bits,sattelite)
    C_ic = twosComp2dec(inp_bits(3,1:16))*2^-29 ;                          % C_ic
    Om_e = twosComp2dec([inp_bits(3,17:24),inp_bits(4,1:24)])*(2^-31)*pi ; % omega_0
    C_is = twosComp2dec(inp_bits(5,1:16))*2^-29 ;                          % C_is
    I_o = twosComp2dec([inp_bits(5,17:24),inp_bits(6,1:24)])*(2^-31)*pi ;  % i_0
    C_rc = twosComp2dec(inp_bits(7,1:16))*2^-5 ;                           % C_rc
    omega = twosComp2dec([inp_bits(7,17:24),inp_bits(8,1:24)])*(2^-31)*pi ;% omega
    Omega = twosComp2dec(inp_bits(9,1:24))*(2^-43)*pi ;                    % omegaDot
    IODE = inp_bits(10,1:8) ;
    idot = twosComp2dec(inp_bits(10,9:22))*(2^-43)*pi ;                    % iDot
    
    e_i = get_ephemeris(sattelite.ephemerises,IODE) ;
    if e_i == -1
        e_data = make_ephemeris_data() ;
        e_data.IODC(3:end) = IODE ;
        sattelite.ephemerises{end+1} = e_data ;
        e_i = length(sattelite.ephemerises) ;
    end
    sattelite.ephemerises{e_i}.C_ic = C_ic ;
    sattelite.ephemerises{e_i}.Om_e = Om_e ;
    sattelite.ephemerises{e_i}.C_is = C_is ;
    sattelite.ephemerises{e_i}.I_o = I_o ;
    sattelite.ephemerises{e_i}.C_rc = C_rc ;
    sattelite.ephemerises{e_i}.omega = omega ;
    sattelite.ephemerises{e_i}.Omega = Omega ;
    sattelite.ephemerises{e_i}.idot = idot ;
    sattelite.ephemerises{e_i}.sf_flags(3) = 1 ;
    
    fprintf(fout, '\t\t\t\tSubframe 3 parising...\n' ) ;
    fprintf(fout, '\t\t\t\t<IODE>\t%s\t%d\n', IODE, bin2dec(IODE) ) ;
    fprintf(fout, '\t\t\t\t<C_ic>\t%f radians \t/* Amplitude of the cosine harmonic correction term to angle of inclination */\n', C_ic ) ;
    fprintf(fout, '\t\t\t\t<Om_e>\t%f semicircles \t/* Longitude of ascending node of orbit plane at weekly epoch */\n', Om_e ) ;
    fprintf(fout, '\t\t\t\t<C_is>\t%f radians \t/* Amplitude of sine harmonic correction term to angle of inclination */\n', C_is ) ;
    fprintf(fout, '\t\t\t\t<I_o>\t%f semicircles \t/* Inclination angle at reference time */\n', I_o ) ;
    fprintf(fout, '\t\t\t\t<C_rc>\t%f meters \t/* amplitude of the sine harmonic correction term to the orbit radius */\n', C_rc ) ;
    fprintf(fout, '\t\t\t\t<omega>\t%f semicircles \t/* Argument of perigee */\n', omega ) ;
    fprintf(fout, '\t\t\t\t<Omega>\t%f semicircles/sec \t/* Rate of right ascension */\n', Omega ) ;
    fprintf(fout, '\t\t\t\t<idot>\t%f semicircles/sec \t/* Rate of inclination angle */\n', idot ) ;

    function s_index = get_sattelite(sattelites, gnssId, svId)
        s_index = -1 ;
        for n=1:length(sattelites)
            if sattelites{n}.gnssId==gnssId && sattelites{n}.svId==svId
                s_index = n ;
            end
        end
        
    function e_index = get_ephemeris(ephemerises,IODE)
        e_index = -1 ;
        for n=1:length(ephemerises)
            if strcmpi(ephemerises{n}.IODC(3:end),IODE)==1
                e_index = n ;
            end
        end
     
    function e_index = get_ephemeris_ontime(ephemerises,rcvTow)
        e_index = -1 ;
        minDelta = 1e50 ;
        for n=1:length(ephemerises)
            if prod(ephemerises{n}.sf_flags(:))~=0 % check if all ephemeris is available
                if abs(rcvTow-ephemerises{n}.toc)<minDelta
                    e_index = n ;
                    minDelta = abs(rcvTow-ephemerises{n}.toc) ;
                end
            end
        end
        
    function s_data = make_sattelite_data()
        s_data.gnssId = -1 ;
        s_data.svId = -1 ;
        s_data.URA = -1 ;
        s_data.SV_HEALTH = 'XXXXXX' ;
        s_data.ephemerises = {} ;
        s_data.measurments = {} ;
        
    function s_eph = make_ephemeris_data()
        s_eph.sf_flags = zeros(1,3) ;
        s_eph.WN = -1 ;
        s_eph.URA = -1 ;
        s_eph.IODC = 'UUUUUUUUUU' ;
        s_eph.Tgd = 0 ;
        s_eph.toc = 0 ;
        s_eph.af2 = 0 ;
        s_eph.af1 = 0 ;
        s_eph.af0 = 0 ;
        s_eph.C_rs = 0 ;
        s_eph.Delta_n = 0 ;
        s_eph.M_0 = 0 ;
        s_eph.C_uc = 0 ;
        s_eph.e_s = 0 ;
        s_eph.C_us = 0 ;
        s_eph.sqrt_as = 0 ;
        s_eph.t_oe = 0 ;
        s_eph.C_ic = 0 ;
        s_eph.Om_e = 0 ;
        s_eph.C_is = 0 ;
        s_eph.I_o = 0 ;
        s_eph.C_rc = 0 ;
        s_eph.omega = 0 ;
        s_eph.Omega = 0 ;
        s_eph.idot = 0 ;
        
   function akos_eph = make_akos_eph(s_eph)
        akos_eph.t_oc = s_eph.toc ;
        akos_eph.a_f2 = s_eph.af2 ;
        akos_eph.a_f1 = s_eph.af1 ;
        akos_eph.a_f0 = s_eph.af0 ;
        akos_eph.T_GD = s_eph.Tgd ;
        akos_eph.sqrtA = s_eph.sqrt_as ;
        akos_eph.deltan = s_eph.Delta_n ;
        akos_eph.M_0 = s_eph.M_0 ;
        akos_eph.e = s_eph.e_s ;
        akos_eph.omega = s_eph.omega ;
        akos_eph.C_uc = s_eph.C_uc ;
        akos_eph.C_us = s_eph.C_us ;
        akos_eph.C_rc = s_eph.C_rc ;
        akos_eph.C_rs = s_eph.C_rs ;
        akos_eph.C_ic = s_eph.C_ic ;
        akos_eph.C_is = s_eph.C_is ;
        akos_eph.i_0 = s_eph.I_o ;
        akos_eph.iDot = s_eph.idot ;
        akos_eph.omega_0 = s_eph.Om_e ;
        akos_eph.omegaDot = s_eph.Omega ;
        akos_eph.t_oe = s_eph.t_oe ;
        
    function s_msr = make_measurment_data()
        s_msr.gnssId = -1 ;
        s_msr.svId = -1 ;
        s_msr.msrTow = 0 ;
        s_msr.week = 0 ;
        s_msr.leapS = 0 ;
        s_msr.prMes = 0 ;
        s_msr.satTow = 0 ;
        s_msr.sat_x = 0 ;
        s_msr.sat_y = 0 ;
        s_msr.sat_z = 0 ;
        s_msr.sat_clk_corr = 0 ;
        s_msr.s_eph = {} ;
        
    function dump_sattelites(fout,sattelites)
        fprintf(fout,'\n') ;
        fprintf(fout,'<SATTELITES LIST DUMP v1.0:>\n') ;
        fprintf(fout,'<LIST SIZE> %d\n', length(sattelites)) ;
        % count # sattelites with valid ephemeris
        validEphmSats = [] ;
        for s_i=1:length(sattelites)
            add_value = 0 ;
            for e_i=1:length(sattelites{s_i}.ephemerises)
                if sattelites{s_i}.ephemerises{e_i}.sf_flags(1)==1 && ... 
                   sattelites{s_i}.ephemerises{e_i}.sf_flags(2)==1 && ...
                   sattelites{s_i}.ephemerises{e_i}.sf_flags(3)==1               
                    add_value = 1 ;
                    break ;
                end
            end
            if add_value==1
                validEphmSats(end+1) = sattelites{s_i}.svId ;
            end
        end
        fprintf(fout, '<SAT WITH VALID EPHM> [') ;
        for nsv=1:length(validEphmSats)
            fprintf(fout, '%d, ', validEphmSats(nsv)) ;
        end
        fprintf(fout, ']\n') ;

        for s_i=1:length(sattelites)
            fprintf(fout, '\t<GPS><svId>%d\n',sattelites{s_i}.svId) ;
            fprintf(fout, '\t\t<URA>%d\n',sattelites{s_i}.URA) ;
            fprintf(fout, '\t\t<SV_HEALTH>%s\n',sattelites{s_i}.SV_HEALTH) ;
            fprintf(fout, '\t\t<EPHEMERISES LIST> %d entries\n', length(sattelites{s_i}.ephemerises) ) ;
            for e_i=1:length(sattelites{s_i}.ephemerises)
                if sattelites{s_i}.ephemerises{e_i}.sf_flags(1)==1
                    fprintf(fout, '\t\t\t<IODC> %s %d\n', sattelites{s_i}.ephemerises{e_i}.IODC, bin2dec(sattelites{s_i}.ephemerises{e_i}.IODC) ) ;
                elseif (sattelites{s_i}.ephemerises{e_i}.sf_flags(2)+sattelites{s_i}.ephemerises{e_i}.sf_flags(3))>0
                    fprintf(fout, '\t\t\t<IODC> %s %d\n', sattelites{s_i}.ephemerises{e_i}.IODC, bin2dec(sattelites{s_i}.ephemerises{e_i}.IODC(3:end)) ) ;
                else
                     fprintf(fout, '\t\t\t<IODC> %s\n', sattelites{s_i}.ephemerises{e_i}.IODC ) ;
               end
                    
                fprintf(fout, '\t\t\t\t<SF1_STATUS> %d\n', sattelites{s_i}.ephemerises{e_i}.sf_flags(1) ) ;
                fprintf(fout, '\t\t\t\t\t<Tgd>\t%f seconds \t/* Sattelite group delay differential */\n', sattelites{s_i}.ephemerises{e_i}.Tgd ) ;    
                fprintf(fout, '\t\t\t\t\t<toc>\t%f seconds \t/* Sattelite clock correction */\n', sattelites{s_i}.ephemerises{e_i}.toc ) ;    
                fprintf(fout, '\t\t\t\t\t<af2>\t%f seconds \t/* Sattelite clock correction */\n', sattelites{s_i}.ephemerises{e_i}.af2 ) ;    
                fprintf(fout, '\t\t\t\t\t<af1>\t%f seconds \t/* Sattelite clock correction */\n', sattelites{s_i}.ephemerises{e_i}.af1 ) ;    
                fprintf(fout, '\t\t\t\t\t<af0>\t%f seconds \t/* Sattelite clock correction */\n', sattelites{s_i}.ephemerises{e_i}.af0 ) ;    
                fprintf(fout, '\t\t\t\t<SF2_STATUS> %d\n', sattelites{s_i}.ephemerises{e_i}.sf_flags(2) ) ;
                fprintf(fout, '\t\t\t\t\t<C_rs>\t%f m \t/* Amplitude of sine harmonic correction terms to the orbit radius */\n', sattelites{s_i}.ephemerises{e_i}.C_rs ) ;
                fprintf(fout, '\t\t\t\t\t<d_n>\t%f semicircle\t/sec /* Mean motion difference from computed value */\n', sattelites{s_i}.ephemerises{e_i}.Delta_n ) ;
                fprintf(fout, '\t\t\t\t\t<M_0>\t%f semicircle\t/* Mean anomaly at reference time */\n', sattelites{s_i}.ephemerises{e_i}.M_0 ) ;
                fprintf(fout, '\t\t\t\t\t<C_uc>\t%f radians\t/* Amplitude of the cosine harmonic correction term to the argument of latitude ... */\n', sattelites{s_i}.ephemerises{e_i}.C_uc ) ;
                fprintf(fout, '\t\t\t\t\t<e_s>\t%f \t/* Orbital eccentricity */\n', sattelites{s_i}.ephemerises{e_i}.e_s ) ;
                fprintf(fout, '\t\t\t\t\t<C_us>\t%f radians\t/* Amplitude of the sine harmonic correction term to the argument of latitude */\n', sattelites{s_i}.ephemerises{e_i}.C_us ) ;
                fprintf(fout, '\t\t\t\t\t<Sqrt(a_s)>\t%f sqrt(m)\t/* Square root of the semimajor axis */\n', sattelites{s_i}.ephemerises{e_i}.sqrt_as ) ;
                fprintf(fout, '\t\t\t\t\t<t_oe>\t%f seconds\t/* reference time of ephemeris */\n', sattelites{s_i}.ephemerises{e_i}.t_oe ) ;
                fprintf(fout, '\t\t\t\t<SF3_STATUS> %d\n', sattelites{s_i}.ephemerises{e_i}.sf_flags(3) ) ;
                fprintf(fout, '\t\t\t\t\t<C_ic>\t%f radians \t/* Amplitude of the cosine harmonic correction term to angle of inclination */\n', sattelites{s_i}.ephemerises{e_i}.C_ic ) ;
                fprintf(fout, '\t\t\t\t\t<Om_e>\t%f semicircles \t/* Longitude of ascending node of orbit plane at weekly epoch */\n', sattelites{s_i}.ephemerises{e_i}.Om_e ) ;
                fprintf(fout, '\t\t\t\t\t<C_is>\t%f radians \t/* Amplitude of sine harmonic correction term to angle of inclination */\n', sattelites{s_i}.ephemerises{e_i}.C_is ) ;
                fprintf(fout, '\t\t\t\t\t<I_o>\t%f semicircles \t/* Inclination angle at reference time */\n', sattelites{s_i}.ephemerises{e_i}.I_o ) ;
                fprintf(fout, '\t\t\t\t\t<C_rc>\t%f meters \t/* amplitude of the sine harmonic correction term to the orbit radius */\n', sattelites{s_i}.ephemerises{e_i}.C_rc ) ;
                fprintf(fout, '\t\t\t\t\t<omega>\t%f semicircles \t/* Argument of perigee */\n', sattelites{s_i}.ephemerises{e_i}.omega ) ;
                fprintf(fout, '\t\t\t\t\t<Omega>\t%f semicircles/sec \t/* Rate of right ascension */\n', sattelites{s_i}.ephemerises{e_i}.Omega ) ;
                fprintf(fout, '\t\t\t\t\t<idot>\t%f semicircles/sec \t/* Rate of inclination angle */\n', sattelites{s_i}.ephemerises{e_i}.idot ) ;
            end
        end

function m_index=get_actual_measurment(measurments_queue,excl_sat_list,gnssId)
    m_index = -1 ;
    maxTow = 0 ;
    for n=1:length(measurments_queue)
        if measurments_queue{n}.gnssId==gnssId
            skip_flag = 0 ;
            for k=1:length(excl_sat_list)
                if measurments_queue{n}.svId==excl_sat_list(k) 
                    skip_flag = 1 ;
                end
            end
            if skip_flag
                continue ;
            end
            if measurments_queue{n}.msrTow>maxTow
                m_index = n ;
                maxTow = measurments_queue{n}.msrTow ;
            end
        end
    end
    
function pvt_solver(fout,measurments_queue,numSat)
    % prepare data for Kai Borre Procedure (https://code.google.com/p/sandiaproject/source/browse/trunk/Docs/CDs/Utah%202008-09/SiGe/GNSS_SDR/geoFunctions/leastSquarePos.m)
    satpos = zeros(3,numSat) ;
    obs = zeros(numSat,1) ;
    
    sat_list = [] ;
    gnssId = 0 ;
    n = 0 ;
    fprintf(fout, '<PVT SOLVER SAT LIST>\n') ;
    while length(sat_list)<numSat
        m_index = get_actual_measurment(measurments_queue,sat_list,gnssId) ;
        if m_index<1
            break ;
        end
        fprintf(fout, '\t<svId>%2d', measurments_queue{m_index}.svId ) ;
        fprintf(fout, '\t\t<SAT POS>(%7.2f,%7.2f,%7.2f)\n', measurments_queue{m_index}.sat_x,measurments_queue{m_index}.sat_y,measurments_queue{m_index}.sat_z) ;
        sat_distance = sqrt(measurments_queue{m_index}.sat_x*measurments_queue{m_index}.sat_x+measurments_queue{m_index}.sat_y*measurments_queue{m_index}.sat_y+measurments_queue{m_index}.sat_z*measurments_queue{m_index}.sat_z) ;
        fprintf(fout, '\t\t\t\t<SAT DISTANCE>%7.2f, <PSEUDORANGE>%7.2f\n', sat_distance, measurments_queue{m_index}.prMes ) ;
        n = n + 1 ;
        satpos(1,n) = measurments_queue{m_index}.sat_x ;
        satpos(2,n) = measurments_queue{m_index}.sat_y ;
        satpos(3,n) = measurments_queue{m_index}.sat_z ;
        obs(n) = measurments_queue{m_index}.prMes + measurments_queue{n}.sat_clk_corr*299792458 ;
        sat_list(end+1) = measurments_queue{m_index}.svId ;
    end
    
    if n==numSat
        settings.c    = 299792458 ;    % The speed of light, [m/s]
        settings.useTropCorr = 1 ;     % Use troposphere correction
        [pos, el, az, dop] = leastSquarePos(satpos, obs, settings) ;
        [phi, lambda, h] = cart2geo(pos(1), pos(2), pos(3), 5 ) ;
        fprintf(fout,'<NAVIGATION SOLUTION><ECEF><%f,%f,%f>, <GEO><%f,%f> <HEIGHT><%f>\n', ...
            pos(1), pos(2), pos(3), phi, lambda, h ) ;
    end
    
    
    