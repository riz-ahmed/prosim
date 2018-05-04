function [R,P,Q,I,M0,Re0st] = RKM45_RPQI(kappa, Rgas, mu0, T0, ps, p0, h0, r0, ra)

% mach's number
M0 = sqrt((2 / (kappa-1)) * ((ps/p0)^((kappa-1)/kappa) - 1));

% reduced reynolds number
% rho0 = p0 / (Rgas*T0);
% u0bar = M0 / (sqrt(kappa*(p0/rho0)));
% Re0 = (rho0 * u0bar * r0) / mu0;    % not sure about mu0 value exactly
% Re0st = Re0 * (h0/r0)^2;
Re0st = (p0*M0*kappa^0.5*h0^2)/((Rgas*T0)^0.5*mu0*r0);

% read m and n
m_I = dlmread('I_m_2010_0bis1.txt');
n_I = dlmread('I_n_2010_0bis1.txt');

% initial conditions
I0 = 1;
Q0 = 1;
R0 = 1;
P0 = 1;

% initializing m and n
m = pchip(m_I(:,2),m_I(:,1),I0);
n = pchip(n_I(:,2),n_I(:,1),I0);

% step-size
N = 1e5;
Ra = ra/r0;
s = (Ra - R0)/N;

% Solution matrix
SM = zeros(N+1,4);
SM(1,:) = [R0 P0 Q0 I0];

% initializing ODE's
R = R0; P = P0; Q = Q0; I = I0; H = 1;

% ODE's
dQ = @(R,Q,P,H,Re0st,m)(m / (H^2*Re0st*P));
dP = @(R,Q,P,H,Re0st,kappa,M0,n) ((-kappa*n*M0^2*Q) / (H^2*Re0st));

% Runge-Kutta integration
for i = 2:N+1
    
    k1 = s * dQ(R,Q,P,H,Re0st,m);
    l1 = s * dP(R,Q,P,H,Re0st,kappa,M0,n);
    
    k2 = s * dQ(R+0.5*s,Q+0.5*k1,P+0.5*l1,H,Re0st,m); % change is H is here constant w.r.t independent varable R
    l2 = s * dP(R+0.5*s,Q+0.5*k1,P+0.5*l1,H,Re0st,kappa,M0,n);
    
    k3 = s * dQ(R+0.5*s,Q+0.5*k2,P+0.5*l2,H,Re0st,m);
    l3 = s * dP(R+0.5*s,Q+0.5*k2,P+0.5*l2,H,Re0st,kappa,M0,n);
    
    k4 = s * dQ(R+s,Q+k3,P+l3,H,Re0st,m);
    l4 = s * dP(R+s,Q+k3,P+l3,H,Re0st,kappa,M0,n);
    
    Q = Q + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 * (1/6)*k4;
    P = P + (1/6)*l1 + (1/3)*l2 + (1/3)*l3 * (1/6)*l4;
    I = 1 / (R*H*P*Q);
    R = R+s;
    
    SM(i,:) = [R P Q I];
end

% fval
R = SM(:,1);
P = SM(:,2);
Q = SM(:,3);
I = SM(:,4);

R(R==0) = [];
P(P==0) = [];
Q(Q==0) = [];
I(I==0) = [];
end
