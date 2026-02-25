% plot and compare

%figure(21);

if 0    
    c_load('diE3p1234')
    figure(21); h = irf_plot(diE3p1234);
    title('diE3p1234')
end

if 1
    c_load('wE3p34')
    c_load('wE3p32')
    c_eval('c_load(''Atwo?'')',sc);
    dwE3p24 = irf_add(1,wE3p34,-1,wE3p32);
    figure(22); h = irf_plot({wE3p32,wE3p34,dwE3p24,diE3sp24,Atwo3});
    ylabel(h(1),'wE3p32');
    ylabel(h(2),'wE3p34');
    ylabel(h(3),'dwE3p24');
    ylabel(h(4),'diE3sp24');
    ylabel(h(5),'Atwo');
end

if 0
    c_load('diEs3p32')
    c_load('diEs3p34')
    diE3sp24 = irf_add(1,diEs3p34,-1,diEs3p32);
    diE3sp24_abs = [diE3sp24(:,1) sqrt(diE3sp24(:,2).^2+diE3sp24(:,3).^2)];
    figure(23); h = irf_plot({diEs3p32,diEs3p34,diE3sp24,diE3sp24_abs});
    ylabel(h(1),'diEs3p32');
    ylabel(h(2),'diEs3p34');
    ylabel(h(3),'diEs3p24');
    ylabel(h(4),'diEs3p24\_abs');
    %title(h(1),'diEs3p32 diEs3p34')
end