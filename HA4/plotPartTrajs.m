
function plotPartTrajs(k, Xk, Xkmin1, Wk, j)
%PLOTPARTTRAJS Summary of this function goes here
%   Plots lines between ith sample of Xk and j(i)th sample of Xk-1. When
%   repeated during a particle filter execution, this will produce particle
%   trajectories illustration over time.
%
%   This function is intended to be passed as a function handle into your
%   particle filter function.
%
% Inputs:
%   k           time instance index
%   Xk          [n x N] N particles of dimension n to approximate p(x_k).
%   Xkmin1      [n x N] N particles of dimension n to approximate p(x_k-1).
%   Wk          [1 x N] Corresponding weights.
%   j           Index vector such that Xk(:,i) = Xkmin1(:,j(i))

    if (size(Xk,2) <= 100) % At most 50 particles may be plotted
        for i = 1:size(Xk,2) % loop through all particles
            x_axis = [k-1 k];
            A = Xkmin1(1,j(i));
            B = Xk(1,i);
            plot(x_axis, [A B]);
            hold on
        end
        title(['Particle trajectories up to time k=', num2str(k)]);
        pause(0.01);
    else
        disp('Too many particles to plot! Taking the 100 highest weights');
        [~, I] = sort(Wk, 'descend');
        for i = 1:100
            x_axis = [k-1 k];
            A = Xkmin1(1,j(I(i)));
            B = Xk(1,I(i));
            plot(x_axis, [A B]);
            hold on
        end
        title(['Particle trajectories up to time k=', num2str(k)]);
        pause(0.01);
    end
end


