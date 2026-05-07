function E_mean = computeEnergy3(x, r, twId, p)
   v = p.v;
    Fs = p.Fs;
    deltaT = p.deltaT;
    t_coda = p.t_coda; 
    
    delta = 0.5; 
    lls = numel(r);
    B0 = x(1); 
    Le = x(2); 
    E_mean = zeros(lls,1); 
    T = r./v;

    for i = 1:lls
        
        if twId == 1
        tw = (T(i)            :1/Fs    :T(i)+deltaT);
        
        elseif twId == 2
        tw = (T(i)+2*deltaT-deltaT : 1/Fs : T(i)+2*deltaT);
        
        else
        tw = (T(i)+3*deltaT-deltaT : 1/Fs : T(i)+3*deltaT);
        
        end
        
        cond_dirac = (tw >= (r(i)/v))  & (tw<=((r(i)/v)+delta));
        dirac = cond_dirac.*(1/delta);

        cond_heav = (tw > (r(i)/v) + delta);
        heav = cond_heav.*1;
        
        E = computeE_Paasschens(r(i), tw, v, B0, Le, dirac, heav);

        tw_coda = (t_coda:1/Fs:t_coda+deltaT); 

        cond_dirac_c = (tw_coda >= (r(i)/v))  & (tw_coda<=((r(i)/v)+delta));
        dirac_c = cond_dirac_c.*(1/delta);
        cond_heav_c = (tw_coda > (r(i)/v) + delta);
        heav_c = cond_heav_c.*1;
        
        E_coda = computeE_Paasschens(r(i), tw_coda, v, B0, Le, dirac_c, heav_c);

        E_norm = 4*pi*(r(i).^2) .*(E./E_coda); 
        E_mean(i) = log10(mean(E_norm)); 

    end
end
