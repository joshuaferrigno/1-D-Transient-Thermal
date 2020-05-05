clear all
close all
clc

eps = 1e-6;
qtp = 10^7; % w/m^3
tinf = 250; % C
k = 30;
alpha = 5e-6;
theta = 1.5; % seconds
h = 1100; 

L = 10e-3; % m
inc = L/5; % dist between nodes;
X = [0 : inc : L];
X = X * 1e3; % m to mm

bi = .0733;
fo = .375;
dtheta = .3; % timestep
tspan = 0 : dtheta : 1.5;

t0(1) = init_temp(inc*0);
t1(1) = init_temp(inc*1);
t2(1) = init_temp(inc*2);
t3(1) = init_temp(inc*3);
t4(1) = init_temp(inc*4);
t5(1) = 340.91;

t = [t0(1) t1(1) t2(1) t3(1) t4(1) t5(1)];

plot(X,t, 'LineWidth',2);
hold on
for x = 2:5000
    
    tinit = [t0(x-1) t1(x-1) t2(x-1) t3(x-1) t4(x-1) t5(x-1)];

    t0(x) = .375 * ((2 * t1(x-1)) + 2.67) + .25 * t0(x-1);
    t1(x) = .375 * (t0(x-1) + t2(x-1) + 2.67) + .25 * t1(x-1);
    t2(x) = .375 * (t1(x-1) + t3(x-1) + 2.67) + .25 * t2(x-1);
    t3(x) = .375 * (t2(x-1) + t4(x-1) + 2.67) + .25 * t3(x-1);
    t4(x) = .375 * (t3(x-1) + t5(x-1) + 2.67) + .25 * t4(x-1);
    t5(x) = .75 * (t4(x-1) + 19.67) + .195 * t5(x-1);
    t = [t0(x) t1(x) t2(x) t3(x) t4(x) t5(x)];
    
    dt = t-tinit;
    maxdt = max(dt);
    
    if maxdt < eps
        fprintf('The maximum timestep until steady state is %d\nThis corresponds to a time of %f seconds\n',x,x*.3);
        fprintf('The final temperatures are as follow:\n')
        fprintf('%3.2f\n',t)
        break
    end
    if x>1 && x<7
        tL = t(1:5);
        plot(X, t, 'LineWidth',1);
    end
end

title 'Transient Node Temperature Over Time'
xlabel 'Displacement [mm]'
ylabel 'Temeperature [C]'
hleg = legend('t = 0.0','t = 0.3','t = 0.6','t = 0.9',' t = 1.2','t = 1.5');
htitle = get(hleg,'Title');
set(htitle,'String','Time of t [seconds]')

figure(2)
hold on
plot(tspan,t0(1:6), 'LineWidth',2)
plot(tspan,t1(1:6), 'LineWidth',2)
plot(tspan,t2(1:6), 'LineWidth',2)
plot(tspan,t3(1:6), 'LineWidth',2)
plot(tspan,t4(1:6), 'LineWidth',2)
plot(tspan,t5(1:6), 'LineWidth',2)
title 'Temperature vs. Time'
xlabel 'Time [s]'
ylabel 'Temeperature [C]'
hleg = legend('Node 0','Node 1','Node 2','Node 3','Node 4','Node 5');

% saveas(figure(1),'Transient_Nodes.jpg')
% saveas(figure(2),'Temp_vs_Time.jpg')

function T = init_temp(x)
    qtp = 10^7; % w/m^3
    k = 30;
    L = 10e-3; % m
    t5(1) = 340.91;

    T = -qtp*x^2/(2*k) + qtp*L^2/(2*k) + t5(1);
end