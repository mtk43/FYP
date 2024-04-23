% Version that is to be used for debugging - Now the main version
close all
%% Script part
%User Inputs
dt = 0.02;              % Timestep size
tmax = 10;              % Simulation time
Ns = 28;                % Number of spatial steps
doPlot = "plot";        % Set equal to "plot" if dynamic plot is wanted, anything else if not
SpatialMethod = "RK4";  % Choose spatial solving scheme: "RK4" or "Euler"
TimeMethod = "BDF1";    % Choose time stepping scheme: "BDF1", "BDF2", "BDFalpha"
alpha = -0.4;           % Set alpha for BDFalpha scheme
%Choose "Straight" or "Static" initial configuration. "Static" uses the
%tension at the first timestep to staticShooting for initial robot position
InterpolateMethod = "Hermite"; %Choose interpolation method: "Linear", "Hermite"

% t = 0:dt:tmax;
% Nt = length(t);
% P_t = zeros(2,Nt);

%Demand positions
% P_t(2,1:Nt/2) = linspace(-0.043948,0.0503,Nt/2);
% P_t(2,Nt/2:end) = linspace(0.0503,-0.043948,Nt/2+1);
% 
% P_t(1,1:Nt/2) = linspace(0,0.0987,Nt/2);
% P_t(1,Nt/2:end) = linspace(0.0987,0,Nt/2+1);
%P_t(2,:) = 0.05*(sin(pi()*t+pi()*3/2)+1) %Uncomment for sinusoidal input
%in y

[solution, T] = dynamicShootingInverseExact(P_t,Ns, dt, t,doPlot, SpatialMethod, TimeMethod, alpha, InterpolateMethod);


% Edits tensions to replace negative and very small values to equivalent positive ones
for i = 1:length(T)
    for ii = 1:4
        if abs(T(ii,i)) < 0.001
            T(ii,i) = 0;
        end
    end
end

T(2,:) = (T(2,:)-T(4,:));
T(4,:) = 0;

for i = 1:length(T)
    if T(1,i) >= 0
        T(1,i) = (T(1,i)-T(3,i));
        T(3,i) = 0;
    else
        T(3,i) = (T(3,i)-T(1,i));
        T(1,i) = 0;
    end
end

function [solution,T_t] = dynamicShootingInverseExact(P_t,Ns, dt, t,doPlot, SpatialMethod, TimeMethod, alpha, InterpolateMethod)
    %%Function that solves the governing equations based upon cosserat theory
    %%of rods and strings using shooting method

    t_steps = length(t);
    options = optimset('Display','off','Algorithm','levenberg-marquardt'); %Turns off fsolves display
    
    %Set Constants
    L = 0.09;                 % Backbone Length (m)
    Area = 1.32e-7;        % Backbone Cross Sectional Area (m^2)
    TotalMass = 6.7e-3;         % Total Continuum Robot Arm Mass (kg)
    E = 75e9;                % Modulus of Elasticity of Backbone (MPa) - Have currently put in Pa value
    Gs = 28.8e9;              % Shear Modulus of Backbone (MPa) - Pa
    g = -9.81;                % Acceleration due to Gravity (m/s^2)
    u_star = [0;0;0];         % Undeformed backbone linear speed configuration
    v_star = [0;0;1];         % Undeformed backbone angular speed configuration
    fe = [0;0;TotalMass/L*g]; % Set Distributed Load on Backbone (Weight)

    %%% changing TotalMass/L*g to rho*Area*g give cantilever rod solution

    le = [0;0;0]; %Set Distributed Moment on Backbone
    Ixx = 1.39e-15; %Second moment of area 2.0106e-14
    Iyy = Ixx;
    Izz = 2.77e-15; %Polar moment of inertia 4.0212e-14
    Kse = diag([Gs*Area,Gs*Area,E*Area]); %Stiffness matrix shear and extension
    Kbt = diag([E*Ixx,E*Iyy,Gs*Izz]); %Stiffness matrix= bending and torsion
    ri = [0 8e-3 0; -8e-3 0 0; 0 -8e-3 0; 8e-3 0 0];
    ri = ri.';
    rho = 6450; %Density
    J = diag([Ixx,Iyy,Izz]);
    Bse = zeros(3); %Strain and extension damping
    Bbt = zeros(3); %Bending and torsional damping
    ds = L/(Ns-1);
    xint = linspace(0,L,Ns);
    xint_interpolate = xint(1:end-1)+ds/2;

    % Set boundary conditions
    p0 = [0;0;0];
    h0 = [0.8256;0;0.8256;0];
    q0 = [0;0;0];
    w0 = [0;0;0];
%     h0 = [0;0;0;0];

    %Initialises first step as straight, or from given static position

    T = zeros(4,t_steps);
    guess = fsolve(@staticShooting,zeros(6,1),options);

    %Guess for n,m, and T
    guess = [y(8:13,1);zeros(4,1)];
    vinit = zeros(3,Ns);
    uinit = zeros(3,Ns);
    T_t = zeros(4,t_steps);
    T_t(:,1) = guess(7:10);

    %Calculates v and u for entire solution
    for i = 1:Ns
        h = y(4:7,i);
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
        vinit(:,i) = Kse\R.'*y(8:10,i) + v_star;
        uinit(:,i) = Kbt\R.'*y(11:13,i) + u_star;
    end

    %Sets initial conditions
    z = [vinit;uinit];
    y = y(1:13,:);
    y = [y ; zeros(6,Ns)];


    y_prev_1 = y(14:19,:);
    z_prev_1 = z;
    y_prev_2 = y(14:19,:);
    z_prev_2 = z;

    switch TimeMethod
        case "BDF1"
            c0 = 1/dt;
            c1 = -1/dt;
            c2 = 0;
            d1 = 0;
        case "BDF2"
            c0 = 1.5/dt;
            c1 = -2/dt;
            c2 = 0.5/dt;
            d1 = 0;
        case "BDFalpha"
            c0 = (1.5 + alpha)/(dt*(1 + alpha));
            c1 = -2/dt;
            c2 = (0.5 + alpha)/(dt*(1 + alpha));
            d1 = alpha/(1 + alpha);
        otherwise
            error("Invalid time stepping method")
    end

    Kse_plus_c0_Bse_inv = (Kse+c0*Bse)^-1;
    Kbt_plus_c0_Bbt_inv = (Kbt+c0*Bbt)^-1;
    Kse_vstar = Kse*v_star;

    solution(:,:,1) = y(1:3,:);
    k = 1;
    visualise(doPlot,k,dt) %Creates plot at current timestep

    yh = c1*y_prev_1+c2*y_prev_2;
    zh = c1*z_prev_1+c2*z_prev_2;
    yt = zeros(6,Ns);
    zt = zeros(6,Ns);

    yh_int = zeros(6,Ns-1);
    zh_int = zeros(6,Ns-1);

    for i = 1:6
        yh_int(i,:) = interpolate(yh(i,:));
        zh_int(i,:) = interpolate(zh(i,:));
    end


    %Loops through time steps setting previous solution for n and m as
    %the next guess
    for k = 2 : t_steps
        p = P_t(:,k);

        %Shooting method solver call
        [guess, ~, exitflag] = fsolve(@dynamicShootingResidual, guess,options);
        if exitflag == 0
            error("Unstable Timestep")
        end
        y_prev_2 = y_prev_1;
        z_prev_2 = z_prev_1;

        y_prev_1 = y(14:19,:);
        z_prev_1 = z;

        yh = c1*y_prev_1+c2*y_prev_2+d1*yt;
        zh = c1*z_prev_1+c2*z_prev_2+d1*zt;

        %tic
        for l = 1:6
            yh_int(l,:) = interpolate(yh(l,:));
            zh_int(l,:) = interpolate(zh(l,:));
        end
        %timenew(k) = toc;

        solution(1:3,:,k) = y(1:3,:);
        T_t(:,k) = guess(7:10);
        visualise(doPlot,k,dt);
    end
    %solution = solution(:,:,t_steps);
    %%
    function int = interpolate(data)
        switch InterpolateMethod
            case "Linear"
                int = interp1(xint,data,xint_interpolate);
            case "Hermite"
                int = pchip(xint,data,xint_interpolate);
            otherwise
                error("Invalid interpolation method")
        end
    end
    function res = staticShooting(guess)
        %Int
        dpi = zeros(3,4); ni = zeros(3,4); mi = zeros(3,4);

        %Guesses for n and m
        n0 = guess(1:3);
        m0 = guess(4:6);

        y0 = [p0
            h0
            n0
            m0];

        %Solve governing equations
        [~,y] = ode45(@staticODE,xint,y0);
        y = y.';

        %Extract boundary values at tip
        nb = y(8:10, Ns);
        mb = y(11:13, Ns);

        %Calculate R from quaternion
        h = y(4:7, Ns);
        h1 = h(1);
        h2 = h(2);
        h3 = h(3);
        h4 = h(4);
        Rb = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];

        vb = Kse\Rb.'*nb + v_star;
        ub = Kbt\Rb.'*mb + u_star;

        %Calculates actual boundary condition for n and m, equations (3.41)
        %and (3.42) of Williams (2020)
        for i = 1:4
            dpi(:,i) = Rb*(skew(ub)*ri(:,i)+vb);

            ni(:,i) = T(i).*dpi(:,i)./norm(dpi(:,i));

            mi(:,i) = T(i)*skew(Rb*ri(:,i))*dpi(:,i)./norm(dpi(:,i));
        end

        %Error between actual BC and calculated BC. Should be equal and
        %opposite
        nb_error = nb + [sum(ni(1,:));sum(ni(2,:));sum(ni(3,:))];
        mb_error = mb + [sum(mi(1,:));sum(mi(2,:));sum(mi(3,:))];

        %Error vector
        res = [nb_error; mb_error];

    end
    %%
    function [dgds] = staticODE(s,r)

        %Initialise
        A = 0; B = 0; G = 0; H = 0; a = 0; b = 0; R = zeros(3,3); dp = zeros(3,4);

        %Collect values from state variables
        h = r(4:7);
        n = r(8:10);
        m = r(11:13);

        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];

        %Calculate v and u from n and m, Eq. 7
        v = Kse\R.'*n + v_star;
        u = Kbt\R.'*m + u_star;

        %Loop to calculate the variables for tendon load and moment
        for i = 1:4

            dp(:,i) = skew(u)*ri(:,i) + v;
            Ai = -T(i).*skew(dp(:,i))^2./norm(dp(:,i))^3;
            Bi = skew(ri(:,i))*Ai;
            A = A + Ai;
            B = B + Bi;
            G = G - Ai*skew(ri(:,i));
            H = H - Bi*skew(ri(:,i));
            a = a + Ai*(skew(u)*dp(:,i));
            b = b + skew(ri(:,i))*Ai*(skew(u)*dp(:,i));

        end

        %Eq. 13
        c = - skew(u)*Kbt*(u - u_star) - skew(v)*Kse*(v-v_star) - R.'*le - b;
        d = - skew(u)*Kse*(v - v_star) - R.'*fe - a;

        %Eq. 12
        dvdu = [Kse+A, G; B, Kbt + H]\[d; c];

        dv = [dvdu(1); dvdu(2); dvdu(3)];
        du = [dvdu(4); dvdu(5); dvdu(6)];

        %Distributed tendon load and moment, Eq. 10
        ft = R*(a + A*dv + G*du);
        lt = R*(b + B*dv + H*du);

        %Calculates state variables, Eq. 17
        dpds = R*v;
        dnds = -fe-ft;
        dmds = -lt-le - cross(dpds,n);
        dhds = [ 0  , -u(1), -u(2), -u(3);
            u(1),   0  ,  u(3), -u(2);
            u(2), -u(3),   0  ,  u(1);
            u(3),  u(2), -u(1),   0  ] * h/2;

        dgds = [dpds;dhds;dnds;dmds];


    end

    %%
    function E = dynamicShootingResidual(guess)
        %Reaction force and moment are guessed
        n0 = guess(1:3);
        m0 = guess(4:6);
        T = guess(7:10);

        y(:,1) = [p0; h0; n0; m0; q0; w0];

        %Fourth-Order Runge-Kutta Integration
        for j = 1 : Ns-1
            switch SpatialMethod
                case "RK4"
                    yj = y(:,j);
                    yhj_int = yh_int(:,j);
                    [k1,z(:,j),yt(:,j),zt(:,j)]=dynamicODE(yj,yh(:,j),zh(:,j));
                    [k2,~]=dynamicODE(yj+k1*ds/2,yhj_int,zh_int(:,j));
                    [k3,~]=dynamicODE(yj+k2*ds/2,yhj_int,zh_int(:,j));
                    [k4,~]=dynamicODE(yj+k3*ds,yh(:,j),zh(:,j+1));
                    y(:,j+1) = yj + ds*(k1 + 2*(k2+k3) + k4)/6;
                case "Euler"
                    yj = y(:,j);
                    [k1,z(:,j),yt(:,j),zt(:,j)]=dynamicODE(yj,yh(:,j),zh(:,j));
                    y(:,j+1) = yj + ds*k1;
                otherwise
                    error("Invalid spatial integration scheme")
            end
        end

        %Cantilever boundary conditions
        pb = y(2:3, Ns);
        nb = y(8:10,Ns);
        mb = y(11:13,Ns);

        %Calculates rotation vector R from quaternion at tip
        h = y(4:7,Ns);
        h1 = h(1);
        h2 = h(2);
        h3 = h(3);
        h4 = h(4);
        Rb = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];

        vb = Kse\Rb.'*nb + v_star;
        ub = Kbt\Rb.'*mb + u_star;

        %Loop for boundary condition of tendon termination
        for i = 1:4
            dpi(:,i) = Rb*(cross(ub,ri(:,i))+vb);

            ni(:,i) = T(i).*dpi(:,i)./norm(dpi(:,i));

            mi(:,i) = T(i)*skew(Rb*ri(:,i))*dpi(:,i)./norm(dpi(:,i));
        end

        %Error between actual BC and calculated BC
        nb_error = -nb - [sum(ni(1,:));sum(ni(2,:));sum(ni(3,:))];
        mb_error = -mb - [sum(mi(1,:));sum(mi(2,:));sum(mi(3,:))];
        pb_error = p - pb;

        %Sets residual
        E = [nb_error; mb_error; pb_error];
    end

    %%

    function [ys,z,yt,zt] = dynamicODE(y,yh,zh)
        %Unpack y and z
        h = y(4:7);
        n = y(8:10);
        m = y(11:13);
        q = y(14:16);
        w = y(17:19);
        vh = zh(1:3);
        uh = zh(4:6);

        %Quaternion to Rotation
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];

        v=Kse_plus_c0_Bse_inv*(R'*n+Kse_vstar-Bse*vh);
        u=Kbt_plus_c0_Bbt_inv*(R'*m-Bbt*uh);
        z=[v;u];

        % Intialise
        A = 0; B = 0; G = 0;
        H = 0; a = 0; b = 0;
        dp = zeros(3,4);

        %Calculate and set time Derivatives
        yt = c0*y(14:19) + yh;
        zt = c0*z + zh;

        vt = zt(1:3);
        ut = zt(4:6);
        qt = yt(1:3);
        wt = yt(4:6);

        for i = 1:4
            dp(:,i) = cross(u,ri(:,i)) + v;
            Ai = -T(i).*skew(dp(:,i))*skew(dp(:,i))./norm(dp(:,i))^3;
            Bi = skew(ri(:,i))*Ai;
            A = A + Ai;
            B = B + Bi;
            G = G - Ai*skew(ri(:,i));
            H = H - Bi*skew(ri(:,i));
            a = a + Ai*cross(u,dp(:,i));
            b = b + cross(ri(:,i),Ai*cross(u,dp(:,i)));

        end

        %Calculate constants c and d
        c = rho*J*wt - cross(u,Kbt*(u - u_star))...
            - cross(v,Kse*(v-v_star)) - R.'*le - b + cross(w,rho*J*w);
        d = rho*Area*qt - cross(u,Kse*(v - v_star))...
            - R.'*fe - a + rho*Area*cross(w,q);

        dvdu = [Kse+A, G; B, Kbt + H]\[d; c];

        dvds = dvdu(1:3);
        duds = dvdu(4:6);

        %Calculates total forces and moments acting along backbone
        ft = R*(a + A*dvds + G*duds);
        lt = R*(b + B*dvds + H*duds);

        %Calculates state vector dpds and dRds
        dpds = R*v;

        %Solves ODEs
        dqds = vt - cross(u,q) + cross(w,v);
        dwds = ut - cross(u,w);
        dnds = rho*Area*R*(cross(w,q) + qt) -fe-ft;
        dmds = R*(cross(w,rho*J*w) + rho*J*wt) -lt-le - cross(dpds,n);

        dhds = [ 0  , -u(1), -u(2), -u(3);
            u(1),   0  ,  u(3), -u(2);
            u(2), -u(3),   0  ,  u(1);
            u(3),  u(2), -u(1),   0  ] * h/2;

        ys = [dpds;dhds;dnds;dmds;dqds;dwds];

    end
    %%
    function visualise(doPlot,k,dt)
        %Plots dynamic robot position in 2D
        if doPlot == "plot"
            plot(y(1,:),y(3,:)); %title ( 'Cantilever Rod');

            xlabel('z (m)'); ylabel('y (m)');
            axis([0 1.1*L -0.8*L 0.8*L]);
            grid on; daspect([1 1 1]);

            legend(string((k-1)*dt))
            drawnow;

        end
    end
    %%
    function X = skew(x)
        %Function to calculate the skew of a vector
        X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];

    end
end
