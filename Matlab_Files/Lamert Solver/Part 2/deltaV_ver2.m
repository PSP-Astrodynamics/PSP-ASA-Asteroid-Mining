function [RA, RB, VA, VB, delVAmag, delVBmag, delVA, delVB] = deltaV_ver2(r1_vec, r2_vec, r1, e_solutions, p_solutions, VE)
% values are in km or km/s
RA = {};
RB = {};

% calculating position vectors for orbits
for index = 1:100
    [x1,y1,z1] = e_solutions(index,1), r1_vec, r2_vec, r1, p_solutions(index,1);
    RA{index} = [x1;y1;z1];
end

for index = 101:200
    [x2,y2,z2] = e_solutions(index,1), r1_vec, r2_vec, r1, p_solutions(index,1);
    RB{index} = [x2;y2;z2];
end

VA = {};
VB = {};

% calculate the velocity at the point t = 0
for index = 1:100
    R1A = RA{index};
    VA{index} = (R1A(:,2) - R1A(:,1))./(365.25*24*3600*index/size(R1A,2));
end

for index = 101:200
    R1B = RB{index};
    VB{index} = (R1B(:,2) - R1B(:,1))./(365.25*24*3600*index/size(R1B,2));
end

delVAmag = {};
delVBmag = {};
delVA = {};
delVB = {};

% calculate the deltaV to get into each orbit at t = 0
for index = 1:100
    VAtrans = VA{index}';
    delVAtrans = VAtrans - VE;
    delVA{index} = delVAtrans;
    delVAmag{index} = norm(delVAtrans);

end

for index = 101:200
    VBtrans = VB{index}';
    delVBtrans = VBtrans - VE;
    delVB{index} = delVBtrans;
    delVBmag{index} = norm(delVBtrans);

end

end