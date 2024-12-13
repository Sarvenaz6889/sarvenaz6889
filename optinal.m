[X, Y] = meshgrid(-5:0.1:5);

% Define the function to plot 
Z = X.^2 + Y.^2 - 10*(cos(2*pi*X) + cos(2*pi*Y)) + 20;
figure;
surf(X, Y, Z);
colormap('jet');
shading interp;
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Surface Plot of Rastrigin Function');
view(-37.5, 30);
axis([-5 5 -5 5 0 100]);
text(-4, 5, 80, 'f(x,y) = x^2 + y^2 - 10(cos(2\pix) + cos(2\piy)) + 20', 'FontSize', 10);