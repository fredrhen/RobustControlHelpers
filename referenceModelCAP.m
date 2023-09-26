function reference_model = referenceModelCAP(CAP, velocity)
    omega_sp = 3;
    zeta_sp = 0.8;
    
    one_t_theta_2 = omega_sp^2 / ((velocity/9.81) * CAP);

    s = tf('s');

    reference_model = (s + one_t_theta_2) / (s^2 + 2*zeta_sp*omega_sp*s + omega_sp^2);

end