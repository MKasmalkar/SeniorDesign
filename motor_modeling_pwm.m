clear
close all
clc

% https://www.mathworks.com/help/sps/ref/pmsm.html
% https://www.mathworks.com/help/mcb/ref/mtpacontrolreference.html
% https://www.mathworks.com/help/mcb/gs/implement-motor-speed-control-by-using-field-oriented-control-foc.html

%% Simulation

% Setup
dT = 5e-6;
Tstop = 10;
tt = 0:dT:Tstop;
NN = length(tt);

ftmr = 200e3;
Ttmr = 1/ftmr;

fsw = 20e3;
Tsw = 1/fsw;

sw_timer = 0;
sw_max = Tsw/dT;

Vsupply = 50; % Assume bipolar supply, e.g. +/- 50V = 100V DC link

% Motor parameters
p = 2;          % # of poles
N = p/2;        % # of pole pairs
Rs = 0.001;      % winding resistance
Ls = 0.0002;
Lm = 0.00002;
Ms = 0.00002;
J = 0.001;
psi_m = 0.18;

% Variables
theta_m = NaN(1, NN);
theta_e = NaN(1, NN);

w_m = NaN(1, NN);
w_dot_m = NaN(1, NN);

T = NaN(1, NN);

L = NaN(3, 3, NN);

v_3ph = NaN(3, NN);
v_pwm = NaN(3, NN);
duty_cycles = NaN(3, NN);

i_3ph = NaN(3, NN);
psi_ph = NaN(3, NN);
psi_phm = NaN(3, NN);

Pin = NaN(1, NN);

i_dq0 = NaN(3, NN);

Pout = NaN(1, NN);

vref_dq0 = zeros(3, NN);

res_emf = NaN(3, NN);

v_g = NaN(1, NN);

% Initialization

theta_m(1) = 0;
theta_e(1) = 0;
w_m(1) = 0;
w_dot_m(1) = 0;
T(1) = 0;

v_3ph(:, 1) = [0 0 0].';
v_pwm(:, 1) = [0 0 0].';
duty_cycles(:, 1) = [0 0 0].';

i_3ph(:, 1) = [0 0 0].';

psi_am_0 = psi_m*cos(theta_e(1));
psi_bm_0 = psi_m*cos(theta_e(1) - 2*pi/3);
psi_cm_0 = psi_m*cos(theta_e(1) + 2*pi/3);
psi_ph(:, 1) = [psi_am_0 psi_bm_0 psi_cm_0].';
psi_phm(:, 1) = [psi_am_0 psi_bm_0 psi_cm_0].';

res_emf(:, 1) = [0 0 0].';

v_g(1) = 0;

Laa_0 = Ls + Lm*cos(2*theta_e(1));
Lbb_0 = Ls + Lm*cos(2*(theta_e(1)-2*pi/3));
Lcc_0 = Ls + Lm*cos(2*(theta_e(1)+2*pi/3));
Lab_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6));
Lba_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6));
Lbc_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6-2*pi/3));
Lcb_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6-2*pi/3));
Lca_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6+2*pi/3));
Lac_0 = -Ms - Lm*cos(2*(theta_e(1)+pi/6+2*pi/3));

L(:, :, 1) = [Laa_0 Lab_0 Lac_0
              Lba_0 Lbb_0 Lbc_0
              Lca_0 Lcb_0 Lcc_0];

% controls

wref = NaN(1, NN);
for n = 1:NN
    if n*dT < 1
        wref(n) = 10*2*pi;
    elseif n*dT < 6
        wref(n) = 20*2*pi;
    elseif n*dT < 8
        wref(n) = 10*2*pi;
    else
        wref(n) = 20*2*pi;
    end
end

kp_wref = 1.0;
ki_wref = 1.0;
wref_error_int = 0;

kp_iq = 1.0;
ki_iq = 1.0;
iq_error_int = 0;

kp_id = 1.0;
ki_id = 1.0;
id_error_int = 0;

kp_iz = 1.0;
ki_iz = 1.0;
iz_error_int = 0;

T_L = 2.5*ones(1, NN);  % Load torque
for n = 1:NN
    if n*dT < 2.5
        T_L(n) = 1;
    elseif n*dT < 5
        T_L(n) = 3;
    elseif n*dT < 7.5
        T_L(n) = 5;
    else
        T_L(n) = 7;
    end
end

for n = 1:NN-1
    % Get variable values at beginning of this time instant

    P = 2/3 * ...
        [cos(theta_e(n))  cos(theta_e(n)-2*pi/3)   cos(theta_e(n)+2*pi/3)
         -sin(theta_e(n)) -sin(theta_e(n)-2*pi/3)  -sin(theta_e(n)+2*pi/3)
               0.5                0.5                        0.5];

    i_dq0(:, n) = P * i_3ph(:, n);
    
    i_d = i_dq0(1, n);
    i_q = i_dq0(2, n);
    i_z = i_dq0(3, n);

    L_d = Ls + Ms + 3/2*Lm;
    L_q = Ls + Ms - 3/2*Lm;

    T(n) = 3/2 * N * (i_q * (i_q * L_d + psi_m) - i_d*i_q*L_q);

    Pout(n) = T(n) * w_m(n);

    Pin(n) = v_pwm(1,n)*i_3ph(1,n) + v_pwm(2,n)*i_3ph(2,n) + v_pwm(3,n)*i_3ph(3,n);

    w_dot_m(n) = (T(n) - T_L(n)) / J;

    % "Simulate" one time step passing
    w_m(n+1) = w_m(n) + w_dot_m(n) * dT;
    theta_m(n+1) = theta_m(n) + w_m(n)*dT;

    % Re-evaluate everything at the end of the time instant
    theta_e(n+1) = theta_m(n+1) * N;
    
    Laa = Ls + Lm*cos(2*theta_e(n+1));
    Lbb = Ls + Lm*cos(2*(theta_e(n+1)-2*pi/3));
    Lcc = Ls + Lm*cos(2*(theta_e(n+1)+2*pi/3));
    Lab = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6));
    Lba = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6));
    Lbc = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6-2*pi/3));
    Lcb = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6-2*pi/3));
    Lca = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6+2*pi/3));
    Lac = -Ms - Lm*cos(2*(theta_e(n+1)+pi/6+2*pi/3));

    L(:, :, n+1) = [Laa Lab Lac
                    Lba Lbb Lbc
                    Lca Lcb Lcc];

    psi_am = psi_m*cos(theta_e(n+1));
    psi_bm = psi_m*cos(theta_e(n+1) - 2*pi/3);
    psi_cm = psi_m*cos(theta_e(n+1) + 2*pi/3);
    psi_phm(:, n+1) = [psi_am psi_bm psi_cm].';

    % FOC
    if mod(n*dT, Ttmr) == 0
        wref_error = wref(n+1) - w_m(n+1);
        wref_error_int = wref_error_int + wref_error*dT;
        wref_pi_out = kp_wref*wref_error + ki_wref*wref_error_int;
    
        idref = 0;
        id_error = idref - i_d;
        id_error_int = id_error_int + id_error*dT;
        id_pi_out = kp_id*id_error + ki_id*id_error_int;
    
        iq_ref = wref_pi_out;
        iq_error = iq_ref - i_q;
        iq_error_int = iq_error_int + iq_error*dT;
        iq_pi_out = kp_iq*iq_error + ki_iq*iq_error_int;
        
        iz_ref = 0;
        iz_error = iz_ref - i_z;
        iz_error_int = iz_error_int + iz_error*dT;
        iz_pi_out = kp_iz*iz_error + ki_iz*iz_error_int;
    
        vref_dq0(1, n+1) = id_pi_out;
        vref_dq0(2, n+1) = iq_pi_out;
        vref_dq0(3, n+1) = 0;
    
        v_3ph(:, n+1) = P^(-1) * [vref_dq0(1, n+1) 
                                  vref_dq0(2, n+1)
                                  vref_dq0(3, n+1)];
    
        duty_cycles(:, n+1) = max(min(v_3ph(:, n+1) / (Vsupply*2) + 1/2, ...
                                  [0.9 0.9 0.9].'), [0.1 0.1 0.1].');
    else
        vref_dq0(:, n+1) = vref_dq0(:, n);
        v_3ph(:, n+1) = v_3ph(:, n);
        duty_cycles(:, n+1) = duty_cycles(:, n);
    end

    sw_timer = sw_timer + 1;

    if sw_timer >= sw_max
        sw_timer = 0;
    end

    % 1 if in on-state
    % -1 if in off-state
    on_or_off = sign(duty_cycles(:, n+1) * sw_max - sw_timer);
    if on_or_off == 0
        on_or_off = -1;
    end

    v_pwm(:, n+1) = on_or_off .* Vsupply;
    

    % Solve simultaneous equations that govern the three-phase load
    % assuming that no current can flow out of ground
    
    % This math is derived from the equations in the PMSM MATLAB
    % article linked above with the added conditions that
    %          v_actual_ph = v_applied_ph - v_g
    % and
    %                i_a + i_b + i_c = 0
    
    % Equations are of the form Ax = b where
    % x = [i_a i_b i_c v_g]^T
    
    A = [
        Rs + Laa/dT         Lab/dT          Lac/dT       1
           Lba/dT        Rs + Lbb/dT        Lbc/dT       1
           Lca/dT           Lcb/dT       Rs + Lcc/dT     1
             1                1              1           0
    ];

    b = [v_pwm(:,n+1) + 1/dT*(psi_ph(:,n) - psi_phm(:,n+1))
                             0                              ];

    res = A^(-1)*b;
    
    i_3ph(1,n+1) = res(1);
    i_3ph(2,n+1) = res(2);
    i_3ph(3,n+1) = res(3);
    v_g(n+1)     = res(4);

    psi_ph(:, n+1) = L(:,:,n+1)*i_3ph(:, n+1) + psi_phm(:, n+1);
end

%% Plotting

start_time = dT;
stop_time = Tstop;

range = floor(start_time/dT):floor(stop_time/dT);

figure
subplot(611)
plot(tt(range), v_3ph(:, range))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('Vph')
title('Sinusoidal phase voltages')

subplot(612)
plot(tt(range), v_pwm(:, range))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('Vpwm')
title('Switch node voltages')

subplot(613)
plot(tt(range), i_3ph(:, range))
xlabel('Time')
ylabel('Iph')
legend('a', 'b', 'c')
title('Phase currents')

subplot(614)
plot(tt(range), mod(theta_m(range)*180/pi, 360))
xlabel('Time')
ylabel('Theta mech (deg)')
title('Angle')

subplot(615)
plot(tt(range), T(range))
hold on
plot(tt(range), T_L(range))
legend('Output', 'Load')
xlabel('Time')
ylabel('T')
title('Torque')

subplot(616)
plot(tt(range), w_m(range) / (2*pi))
hold on
plot(tt(range), wref(range) / (2*pi));
xlabel('Time')
ylabel('w mech (rev/sec)')
legend('Actual', 'Ref')
title('Angular speed')

figure
subplot(511)
plot(tt(range), i_dq0(:, range))
xlabel('Time')
ylabel('dq currents')
legend('d', 'q', 'z')

subplot(512)
plot(tt, Pin)
hold on
avg_elec = movmean(Pin, Tsw/dT);
plot(tt(1:length(avg_elec)), avg_elec)
plot(tt, Pout)
xlabel('Time')
ylabel('Power')
legend('Electrical', 'Avg elec', 'Mechanical')
title('Input/Output Power')

subplot(513)
plot(tt, psi_ph(:, :))
legend('psi_a', 'psi_b', 'psi_c')
xlabel('Time')
ylabel('Wb')
title('Flux')

subplot(514)
plot(tt(1:end-1), diff(psi_ph(1, :))/dT)
hold on
plot(tt(1:end-1), diff(psi_ph(2, :))/dT)
plot(tt(1:end-1), diff(psi_ph(3, :))/dT)
legend('dpsi_a', 'dpsi_b', 'dpsi_c')
xlabel('Time')
ylabel('Wb/sec')
title('dpsi/dt')

subplot(515)
plot(tt, movmean(v_pwm(1, :), Tsw/dT))
hold on
plot(tt, movmean(v_pwm(2, :), Tsw/dT))
plot(tt, movmean(v_pwm(3, :), Tsw/dT))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('Averaged PWM')
title('Averaged PWM')

figure
plot(tt(range), duty_cycles(:, range))
xlabel('Time')
ylabel('Duty cycle')

figure
subplot(511)
plot(tt(1:end-1), v_pwm(1, 1:end-1) - diff(psi_ph(1, :))/dT)
hold on
plot(tt(1:end-1), v_pwm(2, 1:end-1) - diff(psi_ph(2, :))/dT)
plot(tt(1:end-1), v_pwm(3, 1:end-1) - diff(psi_ph(3, :))/dT)
legend('a', 'b', 'c')
xlabel('Time')
ylabel('V')
title('resistor emf (ideal)')
xlim([0, 1000*dT])

subplot(512)
plot(tt, res_emf(1, :))
hold on
plot(tt, res_emf(2, :))
plot(tt, res_emf(3, :))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('V')
title('resistor emf (used)')

subplot(513)
plot(tt, v_pwm(1, :))
hold on
plot(tt, v_pwm(2, :))
plot(tt, v_pwm(3, :))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('voltage')
title('Vpwm')
xlim([0, 1000*dT])

subplot(514)
plot(tt(1:end-1), diff(psi_ph(1, :))/dT)
hold on
plot(tt(1:end-1), diff(psi_ph(2, :))/dT)
plot(tt(1:end-1), diff(psi_ph(3, :))/dT)
legend('a', 'b', 'c')
xlabel('Time')
ylabel('flux')
title('dpsi_dt')
xlim([0, 1000*dT])

subplot(515)
plot(tt, i_3ph(1, :))
hold on
plot(tt, i_3ph(2, :))
plot(tt, i_3ph(3, :))
legend('a', 'b', 'c')
xlabel('Time')
ylabel('amps')
title('current')
xlim([0, 1000*dT])

figure
plot(tt, v_g)
hold on
xlabel('Time')
ylabel('V')
title('Neutral voltage')