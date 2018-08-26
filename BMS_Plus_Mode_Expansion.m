function PHI_exp = ...
    BMS_Plus_Mode_Expansion(PHI,dof_sets_BMS,kappa,R,T_BMS,K_BMS,M_BMS)
             
% Dimitri Krattiger
% 4-4-2017
%
% Description
% ===========
% This function takes the Residual Enhanced Bloch mode synthesis (BMS)
% reduced model mode shapes and expands them to full model size using the
% reduction transformation.
%
% inputs
% ======
% PHI           = Mode Shape in BMS Coordinates 
%                 (n_dof_BMS x n_modes x n_kappa)
% 
% dof_sets_BMS  = DOF set structure for BMS reduction
%
% kappa         = wave vector(s) at which the mode shapes in "PHI" are
%                 obtained
%
% R             = Lattice vectors: [r1,r2,r3]
%
% T_BMS         = Residual Enhanced BMS transformation
%
% T_BMS         = Residual Enhanced BMS transformation
%
% Citation
% ========
% The algorithms contained in this code are described in the following
% references. Please cite them appropriately when using or modifying this 
% code.
%
% [1]   D. Krattiger and M. I. Hussein, Generalized Bloch mode synthesis 
%       for accelerated calculation of elastic band structures, Journal 
%       of Computational Physics, vol. 357, pp. 183?205, Mar. 2018.
%
% [2]   D. Krattiger and M. I. Hussein, Bloch mode synthesis: Ultrafast 
%       methodology for elastic band-structure calculations, Physical 
%       Review E, vol. 90, no. 6, Dec. 2014.
%
% [3]   D. Krattiger, Fast Band-Structure Computation for Phononic and 
%       Electronic Waves in Crystals, PhD Thesis, University of Colorado 
%       at Boulder, 2017.

%% Use the residual-enhanced transformation to expand back to full model size
% ======================================================================= %

% useful dimensions
[~,n_modes,n_kap] = size(PHI);
n_dof_full = size(T_BMS.w0,1);

% periodicity transformation
T_per = Periodic_Boundary_Conditions(dof_sets_BMS);

% preallocate expanded mode shape array
PHI_exp = zeros(n_dof_full,n_modes,n_kap);

T_BMS.w0 = full(T_BMS.w0);
T_BMS.w2 = full(T_BMS.w2);

for j = 1:n_kap

    % phase modification at boundaries due to wavevector
    kvec = kappa(:,j);
    lam  = exp(-1i*R'*kvec);
    lam(end+1:3) = 0;

    % form periodicity transformation
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);

    % Apply Periodicity Transformation to those matrices that are needed
    % for in mode expansion
    Kp.w0 = (T_per_k'*K_BMS.w0*T_per_k);
    Mp.w0 = (T_per_k'*M_BMS.w0*T_per_k);
    
    PHI_exp(:,:,j) = T_BMS.w0*(T_per_k*PHI(:,:,j))+T_BMS.w2*(T_per_k*(Mp.w0\(Kp.w0*PHI(:,:,j))));

end
