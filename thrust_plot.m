clear; clear; close
weit  = 0.62;
distL = 0.05;
Xcg   = 0.01;

data0 = importdata('thrust.dat');
iter0 = data0(:,1);
thrust = data0(:,2);

ww = zeros(size(iter0));
for i = 1: size(iter0)
    ww(i) = weit;
end

figure(2)
plot(iter0, thrust,'b-x', 'LineWidth',2)
hold on
plot(iter0, ww/2.)
% grid on
% xlabel('Iteration')
% ylabel('Trim equation residuals')
legend('Location','southeast')
set(gca,'fontsize',14)

% F = getframe;
% figure(3)
% imshow(F.cdata)