function [R,P,Q,I,I_err,M0,Re0st] = RKM45_RPQI(kappa, Rgas, mu0, T0, ps, p0, h0, r0, ra)

% mach's number
M0 = sqrt((2 / (kappa-1)) * ((ps/p0)^((kappa-1)/kappa) - 1));

% reduced reynolds number
rho0 = p0 / (Rgas*T0);
u0bar = M0 / (sqrt(kappa*(p0/rho0)));
Re0 = (rho0 * u0bar * r0) / mu0;
Re0st = Re0 * (h0/r0)^2;

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
N = 1e4;
Ra = ra/r0;
s = (Ra - R0)/N;

% Solution matrix
SM = zeros(4,N+1);
SM(1,:) = [R0 P0 Q0 I0];

% ODE's
rho = p0 / (Rgas*T0);
RHO = rho/rho0;
dQ = @(R,Q,P,H,Re0st,m)(m / (H^2*Re0st*RHO));
end
