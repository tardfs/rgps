% convert nlib internal ephemeris structure 
% to EASY library format (by Kai Borre)
function easy_eph = eph2easy(s_eph,svId)
easy_eph = zeros(21,1) ;
easy_eph(1)  = svId ;
easy_eph(2)  = s_eph.af2 ;
easy_eph(3)  = s_eph.M_0 ;
easy_eph(4)  = s_eph.sqrt_as ;
easy_eph(5)  = s_eph.Delta_n ;
easy_eph(6)  = s_eph.e_s ;
easy_eph(7)  = s_eph.omega ;
easy_eph(8)  = s_eph.C_uc ;
easy_eph(9)  = s_eph.C_us ;
easy_eph(10) = s_eph.C_rc ;
easy_eph(11) = s_eph.C_rs ;
easy_eph(12) = s_eph.I_o ;
easy_eph(13) = s_eph.idot ;
easy_eph(14) = s_eph.C_ic ;
easy_eph(15) = s_eph.C_is ;
easy_eph(16) = s_eph.Om_e ;
easy_eph(17) = s_eph.Omega ;
easy_eph(18) = s_eph.t_oe ; 
easy_eph(19) = s_eph.af0 ;
easy_eph(20) = s_eph.af1 ;
easy_eph(21) = s_eph.toc ;