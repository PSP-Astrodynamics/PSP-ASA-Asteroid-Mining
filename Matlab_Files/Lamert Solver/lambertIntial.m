function [a_trans, c, s, r1, r2] = lambertIntial(r1_vec,r2_vec,mu_S)
% unit vectors
r1 = norm(r1_vec);
r2 = norm(r2_vec);

c = norm(r2_vec - r1_vec);

s = 1/2 * (r1 + r2 + c);

IP_E = 365.25*24*3600; % s
IP_trans = (IP_E : IP_E : 100*IP_E);

a_trans = (mu_S*(IP_trans./(2.*pi)).^2).^(1/3);
end