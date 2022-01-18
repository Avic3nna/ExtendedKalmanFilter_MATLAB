
% We define a circular path 
L_speed = 0.2;  % Linear speed 0.2 m/s
timestep = 0.5;  % Sensor data update rate
total_distance = 4*pi; %2m radius
step_number = ceil(total_distance/(L_speed*timestep));
A_speed = 2*pi/(step_number*timestep); % Angular speed

% Process noise variance  
Qd = 0.01*L_speed*timestep;
Qb = 0.02*A_speed*timestep;
Qk_1 = [Qd 0; 
        0 Qb];

for k = 1:2*step_number
    pathD(k) = L_speed*timestep;
    pathB(k) = A_speed*timestep;
    pathDNoise(k) = L_speed*timestep + sqrt(Qd)*randn;
    pathBNoise(k) = A_speed*timestep + sqrt(Qb)*randn;
end

% Inicialization
xini = 5;
yini = 2;
thetaini = 0;
Xrealk = [xini; yini; thetaini];
Xk = [6; 3; pi];

Pxini = 0.001;
Pyini = 0.001;
Pthetaini = 0.001;
Pk = [  Pxini 0 0;
        0 Pyini 0 ;
        0 0 Pthetaini];

% Measurements variance
R1 = 0.001;
R2 = 0.001;
R3 = 0.001;
Rk = [  R1 0 0;
        0 R2 0;
        0 0 R3];

% Landmark positions
t1x = 4;
t1y = 8;
t2x = 1;
t2y = 1;
t3x = 11;
t3y = 3;

% Algorithm
Ktotal = zeros(3);      
for l = 1:length(pathD)
    % For the simulation only
    XrealkAUX = Xrealk;
    
    % Real movement of the robot
    Xrealk(1) = XrealkAUX(1) + pathD(l)*cos(XrealkAUX(3)+(pathB(l)/2));
    Xrealk(2) = XrealkAUX(2) + pathD(l)*sin(XrealkAUX(3)+(pathB(l)/2));
    Xrealk(3) = XrealkAUX(3) + pathB(l);
    Xreal(:,l) = Xrealk;  % To keep track of the path 

    % Landmark observation
    Zk = [(atan2(t1y-Xrealk(2),t1x-Xrealk(1)) - Xrealk(3) + sqrt(R1)*randn);
          (atan2(t2y-Xrealk(2),t2x-Xrealk(1)) - Xrealk(3) + sqrt(R2)*randn);
          (atan2(t3y-Xrealk(2),t3x-Xrealk(1)) - Xrealk(3) + sqrt(R3)*randn)];
          
    % To make notation easier and more compact
    Uk = [pathDNoise(l); pathBNoise(l)];

    % Nw cycle, k-1 = k.
    Xk_1 = Xk;
    Pk_1 = Pk;
    
    % Prediction of the state
    X_k = [(Xk_1(1) + Uk(1)*cos(Xk_1(3)+(Uk(2)/2)));
           (Xk_1(2) + Uk(1)*sin(Xk_1(3)+(Uk(2)/2)));
           (Xk_1(3) + Uk(2))];
    
    %%FIXME
    Ak =    [1, 0, -Uk(1)*sin(Uk(2)/2 + Xk_1(3));
             0, 1,  Uk(1)*cos(Uk(2)/2 + Xk_1(3));
             0, 0,                     1];
    %%FIXME
    Bk = [  cos(Uk(2)/2 + Xk_1(3)), -(Uk(1)*sin(Uk(2)/2 + Xk_1(3)))/2;
            sin(Uk(2)/2 + Xk_1(3)),  (Uk(1)*cos(Uk(2)/2 + Xk_1(3)))/2;
                                 0,                         1];
    
    %Ak = jacobian(X_k, [Xk_1(1), Xk_1(2), Xk_1(3)]);  
    %Bk = jacobian(X_k, [Uk(1), Uk(2)]);  
    
    P_k = Ak*Pk_1*((Ak)') + Bk*Qk_1*((Bk)');

    % Prediction of the measurement
    Zk_ = [(atan2(t1y-X_k(2),t1x-X_k(1)) - X_k(3));
          (atan2(t2y-X_k(2),t2x-X_k(1)) - X_k(3));
          (atan2(t3y-X_k(2),t3x-X_k(1)) - X_k(3))];
    
    Hk = [  -(imag(X_k(1)) + real(X_k(2)) - 8)/((imag(X_k(1)) + real(X_k(2)) - 8)^2 + (imag(X_k(2)) - real(X_k(1)) + 4)^2),   -(imag(X_k(2)) - real(X_k(1)) + 4)/((imag(X_k(1)) + real(X_k(2)) - 8)^2 + (imag(X_k(2)) - real(X_k(1)) + 4)^2), -1;
            -(imag(X_k(1)) + real(X_k(2)) - 1)/((imag(X_k(1)) + real(X_k(2)) - 1)^2 + (imag(X_k(2)) - real(X_k(1)) + 1)^2),   -(imag(X_k(2)) - real(X_k(1)) + 1)/((imag(X_k(1)) + real(X_k(2)) - 1)^2 + (imag(X_k(2)) - real(X_k(1)) + 1)^2), -1;
            -(imag(X_k(1)) + real(X_k(2)) - 3)/((imag(X_k(1)) + real(X_k(2)) - 3)^2 + (imag(X_k(2)) - real(X_k(1)) + 11)^2), -(imag(X_k(2)) - real(X_k(1)) + 11)/((imag(X_k(1)) + real(X_k(2)) - 3)^2 + (imag(X_k(2)) - real(X_k(1)) + 11)^2), -1];
    %Hk = jacobian(ZK_, [X_k(1), X_k(2), X_k(3)]);
    
    % Comparison
    Yk = Zk-Zk_; %innovation matrix
    for r=1:3
        if Yk(r)>pi
            Yk(r) = Yk(r) - 2*pi;
        end
        if Yk(r)<(-pi)
            Yk(r) = Yk(r) + 2*pi;
        end
    end
    Sk = Hk*P_k*((Hk)') + Rk;
    Wk = P_k*((Hk)')*inv(Sk);

    % Correction
    Xk = X_k + Wk*Yk; %%FIXME
    Pk = (eye(3) - Wk*Hk)*Pk;  %%FIXME
    
    %To store the result
    Xestimated(:,l) = Xk;
end 

% Drawing
figure(1);
subplot(2,2,1);
axis([0 12 0 9])
hold on

for m = 1:length(pathD)
    plot(Xreal(1,m),Xreal(2,m),'.k', Xestimated(1,m),Xestimated(2,m),'.r');
end
legend('Real movement','Estimation')

% Landmarks
plot(t1x,t1y,'sk');
plot(t2x,t2y,'sk');
plot(t3x,t3y,'sk');
xlabel ('X (m)')
ylabel ('Y (m)')

subplot(2,2,2);
plot(Xreal(1,:));
hold on
plot(Xestimated(1,:),'g');
xlabel ('t (timesteps)')
ylabel ('X (m)')

subplot(2,2,3);
plot(Xreal(2,:));
hold on
plot(Xestimated(2,:),'g');
xlabel ('t (timesteps)')
ylabel ('Y (m)')

subplot(2,2,4);
plot(Xreal(3,:));
hold on
plot(Xestimated(3,:),'g');
xlabel ('t (timesteps)')
ylabel ('\theta (rad)')
hold off

pause
close all
clear all

