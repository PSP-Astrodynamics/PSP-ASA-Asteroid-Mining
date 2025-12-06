function [RA, RB, VA, VB, delVAmag, delVBmag, delVA, delVB] = deltaV(r1_vec, r2_vec, r1, e_solutions, p_solutions, VE)
RA = {};
RB = {};

for index = 1:100
    [x1,y1,z1] = plotCycler(e_solutions(index,1), r1_vec, r2_vec, r1, p_solutions(index,1), '#4DBEEE');
    [x2,y2,z2] = plotCycler(e_solutions(index,2), r1_vec, r2_vec, r1, p_solutions(index,2), '#1b10c2');
    RA{index} = [x1;y1;z1];
    RB{index} = [x2;y2;z2];
end

VA = {};
VB = {};

for index = 1:100
    R1A = RA{index};
    R1B = RB{index};
    VA{index} = (R1A(:,2) - R1A(:,1))./(365.25*24*3600*index/size(R1A,2));
    VB{index} = (R1B(:,2) - R1B(:,1))./(365.25*24*3600*index/size(R1B,2));
end

delVAmag = {};
delVBmag = {};
delVA = {};
delVB = {};

for index = 1:100
    VAtrans = VA{index}';
    VBtrans = VB{index}';
    delVAtrans = VAtrans - VE;
    delVBtrans = VBtrans - VE;
    delVA{index} = delVAtrans;
    delVB{index} = delVBtrans;
    delVAmag{index} = norm(delVAtrans);
    delVBmag{index} = norm(delVBtrans);

end
end