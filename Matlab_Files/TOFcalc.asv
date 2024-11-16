a_E = 149597898;
mu_sun = 132712440017.99;
N = 150;
TOF_hoh = 8.6176e+07;
period_E = 31558205; % sec

r1_vec = [124175791.663549, 75089262.231372, 32549514.672663]; %Intial position vector of sataliete km
r2_vec = [-545855777.337076, -386372825.992004, -400476802.487688]; %Sataliet final position km

r1 = norm(r1_vec);
r2 = norm(r2_vec);

c = norm(r2_vec - r1_vec);

s = 1/2 * (r1 + r2 + c);

a_trans = N^(2/3)*(a_E);

alpha = 2*asin(sqrt(s/(2*a_trans)));

beta = 2*asin(sqrt((s-c)/(2*a_trans)));

TOF = ((a_trans^(3/2))*((alpha-beta)-(sin(alpha)-sin(beta))))/(sqrt(mu_sun))