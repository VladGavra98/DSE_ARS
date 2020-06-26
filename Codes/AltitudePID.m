Kp = 16;
u  = 0.75;
Ki = 0.8;
Kd = 1.5;
m = 15;
g = 9.81;
%Kt = m*g/u;
a = 12.38;

% lvl  = 1*ones(1,100);
% lvl2 = 2*ones(1,100);
% stepup = 1:0.01:2;
% href = [0:0.01:1 lvl stepup lvl2];

t = 0:0.01:4.01;
%[ b, i ] = min( abs( t-0.314465 ) );
href = ones(1,402);
s = 0.5*a*t.^2;
[ b, i ] = min( abs( s-1 ) );

href(1:i) = s(1:i);
% % % size(href(150:1+i))
% % % size(href(1:i))
href(150:149+i) = href(1:i)+1;
href(149+i:end) = 2;

h = tf([Kp Ki],[u/g Kd Kp Ki]);

[y, t, x] = lsim(h, href, t);
overshoot = ((max(y)-2))*100;
idx       = find(y==max(y));
osx       = t(idx-20);

%[ d, idx2 ] = min( abs( y-1 ) );
%[ q, idx3 ] = min( abs( href-1 ) );
% idx3      = find(href==0.5);
frst      = y(1:100);
[ d, idx2 ] = min( abs( frst-1 ) );
%idx2      = find(y(1:50)==1);

lag       = t(idx2) - t(i);
disp(frst(end))
lagx      = t(idx2);
figure(1)
plot(t, href, 'b');
hold on
plot(t, y, 'r');
txt = ['max. overshoot:  ', num2str(overshoot) '% \rightarrow'];
text(osx,max(y),txt,'HorizontalAlignment','right')
txt = ['\leftarrow '];
text(lagx,1,txt,'HorizontalAlignment','left')
txt = ['time lag:  ', num2str(lag) '[s]' ];
text(lagx,0.9,txt,'HorizontalAlignment','left')
hold off
title('Altitude Controller Simulation');
xlabel('Time (s)');
ylabel('Altitude (m)');
ylim([0 3]);
legend('Desired Height','UAS Height');


