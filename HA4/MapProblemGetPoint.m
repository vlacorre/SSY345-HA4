
clear all
%This file draw a map of the village, and allows us to manually
%draw the trajectory of the vehicle.

figure(1)
clf
hold on
plot([1+j 1+9*j 5+9*j])
plot([7+9*j 11+9*j 11+j 7+j]);plot([5+j 1+j])
plot([2+5.2*j 2+8.3*j 4+8.3*j 4+5.2*j 2+5.2*j])%House 1
plot([2+3.7*j 2+4.4*j 4+4.4*j 4+3.7*j 2+3.7*j])%House 2
plot([2+2*j 2+3.2*j 4+3.2*j 4+2*j 2+2*j])%House 3
plot([5+j 5+2.2*j 7+2.2*j 7+j])%House 4
plot([5+2.8*j 5+5.5*j 7+5.5*j 7+2.8*j 5+2.8*j])%House 5
plot([5+6.2*j 5+9*j]);plot([7+9*j 7+6.2*j 5+6.2*j])%House 6
plot([8+4.6*j 8+8.4*j 10+8.4*j 10+4.6*j 8+4.6*j])%House 7
plot([8+2.4*j 8+4*j 10+4*j 10+2.4*j 8+2.4*j])%House 8
plot([8+1.7*j 8+1.8*j 10+1.8*j 10+1.7*j 8+1.7*j])%House 9

axis([0.8 11.2 0.8 9.2])
title('A map of the village','FontSize',20)

disp('Start clicking in the village to create a trajectory!')
disp('Press "Return" to finish.')

[X,Y]=ginput;
plot([X+Y*j],'-*')

Xk = [X';Y']
save Xk
