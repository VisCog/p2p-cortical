pdList =  0.005:.005:.1;
for p =1:length(pdList)
    pd = pdList(p);  
    p_nd = 1;
    for np = 1:800
   p_nd= p_nd.*(1-pd);
   saved_p_nd(p, np) = p_nd;
    end
end
plot(1:800, 1-saved_p_nd')