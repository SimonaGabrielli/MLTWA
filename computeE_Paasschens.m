function Erow = computeE_Paasschens(r_i,tw,v,B0,Le,dirac,heav)
    r_isq = r_i^2;
    vtw = v*tw;
    vtw_sq = vtw.^2; 
    gs = 4*pi*v*r_isq;
    expvtwLe = exp(-vtw.*Le); 
    fourpivtw = 4*pi*vtw; 
    threeLeB0 = 3*Le*B0; 
    oneongs = 1./gs; 
    
    bnum1 = (vtw_sq-r_isq)./(vtw_sq); 
    bnum1 = single(bnum1); 
    bnum = bnum1.^(single(1/8));   
    
    % bnum = bnum1.^(1/8);
    %bden =(((fourpivtw)./(threeLeB0)).*((fourpivtw)./(threeLeB0)).*((fourpivtw)./(threeLeB0))).^(1/2);
    
    bden =sqrt(((fourpivtw)./(threeLeB0)).*((fourpivtw)./(threeLeB0)).*((fourpivtw)./(threeLeB0)));
    
    xx =vtw *Le*B0.*(bnum.^6);
    % xx = vtw *Le*B0.*(bnum.*bnum.*bnum.*bnum.*bnum.*bnum);

    % Erow  = (1./gs).*(exp(-vtw.*Le)).*dirac + ((bnum./bden).*(exp(-vtw.*Le)).*((exp(xx).*sqrt(1+(2.026./(xx+eps)))).*heav));
    
    %Erow  = (1./gs).*(expvtwLe).*dirac + ((bnum./bden).*(expvtwLe).*((exp(xx).*sqrt(1+(2.026./(xx+eps)))).*heav));
    
    %a = (1./gs).*(expvtwLe).*dirac;
    
    a = oneongs.*(expvtwLe);
    b = (bnum./bden).*(expvtwLe); 
    argumentsqrt = 1+(2.026./(xx+eps)) ; 
    %c = exp(xx).*sqrt(1+(2.026./(xx+eps))).*heav; 
    c = exp(xx).*sqrt(argumentsqrt); 
    Erow = a.*single(dirac) + b.*c.*heav; 

    %Erow  = (1./gs).*(expvtwLe).*dirac + ((bnum./bden).*(expvtwLe).*((exp(xx).*sqrt(1+(2.026./(xx+eps)))).*heav));

end