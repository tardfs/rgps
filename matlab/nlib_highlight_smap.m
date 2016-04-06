function [hl_smap,sat_list] = nlib_highlight_smap(smap, El)
hl_smap = smap + 1 ;
sat_list = zeros(128,1) ;
satNum = 0 ;
for n=1:size(smap,1)
    if nnz(smap(n,:)==1)==size(smap,2)
        satNum = satNum + 1 ;
        sat_list(satNum) = n ;
        hl_smap(n,:) = 3 ;
    end
end
hl_smap((El<20) & (smap~=0)) = 4 ;

sat_list = sat_list(1:satNum) ;