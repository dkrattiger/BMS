
function [omega, Phi] = ...
    RBME(K_free, M_free, X, R, kappaPts, varargin)
     
% THIS FUNCTION NEEDS TO BE TESTED (recently updated to simplify input and
% code itself, but it probably doesn't run yet)


% Reduced Bloch Mode Expansion dispersion computation code
% Dimitri Krattiger
% 1-27-2014
% inputs
% ======
% K_free = free stiffness matrix
% 
% M_free = free mass matrix
% 
% n_exp = number of expansion points
% 
% m_exp = number of expansion modes
% 
% kappaPts = These are the "vertices" in wave vector space that will be
% connected in order to compute the band-structure. Usually these are
% high symmetry points.
%
% X  = node coordinates
% R  = lattice vectors
%
% Citation
% ========
% The algorithms contained in this code are described in the following
% reference. Please cite it appropriately when using or modifying this 
% code.
%
% [1]   M. I. Hussein, "Reduced Bloch mode expansion for periodic media 
%       band structure calculations," Proceedings of the Royal Society 
%       A: Mathematical, Physical and Engineering Sciences, vol. 465, 
%       no. 2109, pp. 2825?2848, Sep. 2009.
%

%% Check what inputs are given and set rest to default
% ======================================================================= %

% set_default values for options
defaults.n_curves = 10;
defaults.m_exp = [];
defaults.n_exp = 2^1+1;
defaults.n_kap_seg = 2^4+1;

if nargin>=6
    options = varargin{1};
end

% fill in unspecified options fields with the default values
options = setstructfields(defaults,options);

% If m_exp is empty set it to n_curves.
if isempty(options.m_exp)
    options.m_exp = options.n_curves;
end

%% Basic Setup

% number of full model DOFs
n_dof = size(K_free,1);

% number of high symmetry points
n_segs = size(kappaPts,2)-1;

% number of dimensions and number of directions of periodicity
[n_dim,n_per] = size(R);

% unpackage options variables (for more readable code)
n_exp = options.n_exp;
m_exp = options.m_exp;
n_kap_seg = options.n_kap_seg;

%% Use geometry to create indices for interior and boundary DOF sets
% ======================================================================= %

% use node coordinates and lattice vectors to get indices of node sets 
% that are on boundaries. Output is a structure containing fields for each
% boundary set
node_sets = find_node_sets(X,R);

% convert node index sets to DOF index sets. Output is a structure
% containing same fields as "node_sets" structure.
dof_sets = node2dof(node_sets,n_dpn);


%% loop through segments of brillouin zone and compute RBME transformation
% ======================================================================= %

% loop through segments
for i = 1:n_segs
    
    % output loop information
    disp(['BZ Segment ', num2str(i), ' full model calculations'])
    
    % compute k-point vectors for expansion points in current BZ segment
    kappaExp{i} = kappaPts(:, i)*ones(1, n_exp) + ...
        (kappaPts(:, i+1) - kappaPts(:,i))*linspace(0, 1, n_exp);
    
    % compute full model band-structure solution at the expansion points
    solveoptions.n_curves = options.m_exp;
    if i==1
        [omegaExp{i}, PhiExp{i}, ~] = ...
            dispersion_solver_w_k(kappaExp{i}, K_free, M_free, dof_sets, R, solveoptions);
    else
        % skip first point because it was calculated in previous segment
        [omegaExp{i}, PhiExp{i}, ~] = ...
            dispersion_solver_w_k(kappaExp{i}(:,2:end), K_free, M_free, dof_sets, R, solveoptions);

        % add in skipped solution from previous segment
        omegaExp{i} = [omegaExp{i-1}(end, :), omegaExp{i}];
        PhiExp{i} = cat(3, PhiExp{i-1}(:, :, end), PhiExp{i-1});
    end

    % concatenate eigenvectors (Bloch modes) from each k-point into a
    % single transformation matrix
    T_RBME{i} = reshape(PhiExp{i},[size(PhiExp{i},1), n_exp*m_exp]);
    [T_RBME{i},~] = qr(T_RBME{i},0);
    
end


%% loop through segments of brillouin zone and compute dispersion
% ======================================================================= %
omega = [];
Phi = [];

for i = 1:n_segs
    
    % output loop information
    disp(['BZ Segment ', num2str(i), ' RBME reduced model calculations'])
    
    % compute k-point vectors for expansion points in current BZ segment
    kappa{i} = kappaPts(:, i)*ones(1, n_kap_seg) + ...
        (kappaPts(:, i+1) - kappaPts(:,i))*linspace(0, 1, n_kap_seg);
    
    % compute full model band-structure solution at the expansion points
    solveoptions.n_curves = options.n_curves;
    solveoptions.RBME_Transformation = T_RBME{i};
    if i==1
        [omegaTemp, PhiTemp, ~] = ...
            dispersion_solver_w_k(kappa{i}, K_free, M_free, dof_sets, R, solveoptions);
    else
        % skip first point because it was calculated in previous segment
        [omegaTemp, PhiTemp, ~] = ...
            dispersion_solver_w_k(kappa{i}(:,2:end), K_free, M_free, dof_sets, R, solveoptions);

    end
    
    % expand Phi to full size
    PhiTemp = T_RBME{i}*PhiTemp;
    
    % add in new solutions
    omega = [omega,omegaTemp];
    Phi = cat(3, Phi, PhiTemp);
end