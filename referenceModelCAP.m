function reference_model = referenceModelCAP(CAP, velocity)
    omega_sp = 5;
    zeta_sp = 1;
    
    one_t_theta_2 = omega_sp^2 / ((velocity/9.81) * CAP);

    s = tf('s');

    reference_model = (s + one_t_theta_2) / (s^2 + 2*zeta_sp*omega_sp*s + omega_sp^2);
    reference_model = reference_model /dcgain(reference_model);

end