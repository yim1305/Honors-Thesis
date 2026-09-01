clear; clc; close all;

configs = {
    'Rooftop', [0, sin(pi/4), cos(pi/4); 0, -sin(pi/4), cos(pi/4); 0, sin(pi/4), cos(pi/4); 0, -sin(pi/4), cos(pi/4)];
    'Box', [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0];
    'Pyramid', [sin(54.74*pi/180), 0, cos(54.74*pi/180); 0, sin(54.74*pi/180), cos(54.74*pi/180); -sin(54.74*pi/180), 0, cos(54.74*pi/180); 0, -sin(54.74*pi/180), cos(54.74*pi/180)]
};

% Create Spherical Grid
[az, el] = meshgrid(linspace(0, 2*pi, 80), linspace(-pi/2, pi/2, 40)); 
ux = cos(el) .* cos(az);
uy = cos(el) .* sin(az);
uz = sin(el);

for c = 1:size(configs, 1)
    configName = configs{c, 1};
    g = configs{c, 2};

    % Setup Figure 
    fig = figure('Name', configName, 'Color', 'w', 'Position', [50 + (c-1)*50, 100, 800, 400]);
    
    % Left Subplot Setup
    subplot(1, 2, 1);
    hold on; grid on; view(3); axis equal;
    xlabel('h_x'); ylabel('h_y'); zlabel('h_z');
    title([configName, ' - External Envelope']);
    set(gca, 'FontSize', 10, 'GridAlpha', 0.15);
    
    % Right Subplot Setup
    subplot(1, 2, 2);
    hold on; grid on; view(3); axis equal;
    xlabel('h_x'); ylabel('h_y'); zlabel('h_z');
    title([configName, ' - Internal Singularities']);
    set(gca, 'FontSize', 10, 'GridAlpha', 0.15);
    
    % Initialize arrays 
    hx_all = []; hy_all = []; hz_all = [];
    
    % 4 Independent CMGs 
    h_base = zeros(size(ux,1), size(ux,2), 3, 4);
    for i = 1:4
        gx = g(i,1); gy = g(i,2); gz = g(i,3);
        udotg = ux.*gx + uy.*gy + uz.*gz;
        px = ux - udotg .* gx;
        py = uy - udotg .* gy;
        pz = uz - udotg .* gz;
        norm_p = sqrt(px.^2 + py.^2 + pz.^2);
        norm_p(norm_p == 0) = 1e-6;
        h_base(:,:,1,i) = px ./ norm_p;
        h_base(:,:,2,i) = py ./ norm_p;
        h_base(:,:,3,i) = pz ./ norm_p;
    end

    epsilons = dec2bin(0:15) - '0';
    epsilons = epsilons * 2 - 1;

    for k = 1:16
        eps = epsilons(k, :);
        hx = eps(1)*h_base(:,:,1,1) + eps(2)*h_base(:,:,1,2) + eps(3)*h_base(:,:,1,3) + eps(4)*h_base(:,:,1,4);
        hy = eps(1)*h_base(:,:,2,1) + eps(2)*h_base(:,:,2,2) + eps(3)*h_base(:,:,2,3) + eps(4)*h_base(:,:,2,4);
        hz = eps(1)*h_base(:,:,3,1) + eps(2)*h_base(:,:,3,2) + eps(3)*h_base(:,:,3,3) + eps(4)*h_base(:,:,3,4);

        % Collect points for the external envelope
        hx_all = [hx_all; hx(:)];
        hy_all = [hy_all; hy(:)];
        hz_all = [hz_all; hz(:)];

        % Plot ALL surfaces as internal structures
        subplot(1, 2, 2);
        surf(hx, hy, hz, hz, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end

    % Generate External Envelope
    subplot(1, 2, 1);
    K = convhull(hx_all, hy_all, hz_all);
    
    % Plot the boundary mesh 
    trisurf(K, hx_all, hy_all, hz_all, hz_all, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
    
    for sp = 1:2
        subplot(1, 2, sp);
        colormap(gca, 'turbo');
        camlight('headlight');
        camlight('right');
        lighting gouraud;
        material shiny;
        view(45, 30);
    end
end