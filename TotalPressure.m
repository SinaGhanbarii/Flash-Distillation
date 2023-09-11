function output = TotalPressure(A1,B1,C1,A2,B2,C2,T,P,x)
output =  (1-x) * 10^(A1-B1/(T+C1)) + (x) * 10^(A2-B2/(T+C2)) - P;
end