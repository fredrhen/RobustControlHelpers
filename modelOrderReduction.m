function [LOES] = modelOrderReduction(HOS,options)
    arguments    
        HOS
        options.freqRange (1,:) = [0.1 10];
        options.orderLOES (1,1) {mustBeNumeric} = 2
        options.plot = 0
    end
    R = reducespec(HOS,"balanced");
    R.Options.Goal = "absolute"; 
    R.Options.FreqIntervals = options.freqRange;

    LOES = getrom(R,Order=2,Method="truncate");

    if options.plot
        figure; bodemag(LOES,HOS); legend("LOES","HOS","Location","best"); xlim(options.freqRange);
        figure; stepplot(LOES,HOS); legend("LOES","HOS","Location","best"); 

        s = tf('s');
        mag_upper_bound = (3.16*s^2 + 31.61*s + 22.79) / (s^2 + 27.14*s + 1.84);
        mag_lower_bound = (0.0955*s^2 + 9.92*s + 2.15) / (s^2 + 11.6*s + 4.96);

        dif = LOES / HOS;

        figure; bodemag(dif, 'r', mag_lower_bound, 'r--', mag_upper_bound, 'r--');
        legend("LOES", "Bounds"); 
        grid on; 
        
        phase_upper_bound = (68.89*s^2+1100.12*s-275.22*exp(-0.0059*s)) / (s^2+39.94*s+9.99);
        phase_lower_bound = (475.32*s^2+184100*s+29456.1*exp(-0.0072*s)) / (s^2+11.66*s+0.0389);
        figure;  h = bodeplot(dif,'b',phase_upper_bound,'b--',phase_lower_bound,'b--');
        legend("LOES", "Bounds"); setoptions(h,'MagVisible','off'); grid on; 
    end
end


