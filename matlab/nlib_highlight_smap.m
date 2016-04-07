function [hl_smap,sat_list] = nlib_highlight_smap(settings, smap, El)
hl_smap = smap + 1 ;
sat_list = zeros(128,1) ;
satNum = 0 ;
for n=1:size(smap,1)
    if nnz(smap(n,:)==1)==size(smap,2)
        hl_smap(n,:) = 3 ;

        if (nnz(El(n,:)>settings.minSatElevation)==size(El,2))
            satNum = satNum + 1 ;
            sat_list(satNum) = n ;
        end

    end
end
hl_smap((El<settings.minSatElevation) & (smap~=0)) = 4 ;

sat_list = sat_list(1:satNum) ;