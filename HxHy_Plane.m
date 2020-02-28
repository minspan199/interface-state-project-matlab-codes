figure
t = -pi:pi/100:pi;
ta = 0.02;tb = 0.008;r1 = ta - tb;r2 = ta + tb;
x = ta + tb*cos(t);
y = tb*sin(t);

plot(x,y)
axis([-0.02 0.04 -0.03 0.03])
set(gcf, 'Position', [00, 00, 350, 300])
box on
hold on
plot(ta,0,'b.')
text(ta,0,'(ta,0)','color','b');
plot(0,0,'b*')

x = r1*cos(t);
y = r1*sin(t);
plot(x,y,'b--')
x = r2*cos(t);
y = r2*sin(t);
plot(x,y,'b--')
