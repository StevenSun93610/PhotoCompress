function [] = compress()
A=imread('WallOfWindows.jpg','jpg');
% disp(image(A));
A_red = double(A(:, :, 1));
A_green = double(A(:, :, 2));
A_blue = double(A(:, :, 3));

r_values = [10, 20, 30, 100];
tol = 1000;

tic;

for i = 1:length(r_values)
    % compute the SVD of each color component matrix using block power svd
    r = r_values(i);
    [S_red, ~, U_red, V_red] = block_power_svd(A_red, r, tol);
    [S_green, ~, U_green, V_green] = block_power_svd(A_green, r, tol);
    [S_blue, ~, U_blue, V_blue] = block_power_svd(A_blue, r, tol);

    % obtain the rank-r approximation of each color component matrix
    A_red_r = U_red * S_red * V_red';
    A_green_r = U_green * S_green * V_green';
    A_blue_r = U_blue * S_blue * V_blue';

    % recombine the rank-r approximations into a single uint8 format matrix
    A_r = uint8(cat(3, A_red_r, A_green_r, A_blue_r));

    % display the resulting image using the image function
    figure;
    image(A_r);
    title(['Rank-' num2str(r_values(i)) ' Approximation block power svd version with tolerence of ' num2str(tol)]);
end

elapsed_time = toc;
fprintf('It took %e seconds with tolerance of %d\n', elapsed_time, tol);

end







