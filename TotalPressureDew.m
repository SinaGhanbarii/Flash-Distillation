function output = TotalPressureDew(A1,B1,C1,A2,B2,C2,T,P,y)
output = P*((y/(10^(A1 - B1/(T+C1))))+((1-y)/10^(A2-B2/(T+C2)))) - 1;
end