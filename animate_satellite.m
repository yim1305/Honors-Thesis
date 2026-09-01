function animate_satellite(q_store, r_ECI_store, tspan, skip, sat)

figure('Color','k','Position',[100 100 900 700]);
ax = axes;
set(ax,'Color','k','XColor','w','YColor','w','ZColor','w');
hold on; grid on; axis equal;
xlim([-0.7 0.7]); ylim([-0.7 0.7]); zlim([-0.7 0.7]);
xlabel('X [m]','Color','w','FontSize',11);
ylabel('Y [m]','Color','w','FontSize',11);
zlabel('Z [m]','Color','w','FontSize',11);
title('12U CubeSat Attitude Animation','Color','w','FontSize',13);

%% Inertial frame
L = 0.55;
quiver3(0,0,0,L,0,0,'w--','LineWidth',1,'MaxHeadSize',0.3,'AutoScale','off');
quiver3(0,0,0,0,L,0,'w--','LineWidth',1,'MaxHeadSize',0.3,'AutoScale','off');
quiver3(0,0,0,0,0,L,'w--','LineWidth',1,'MaxHeadSize',0.3,'AutoScale','off');
text(L+0.02,0,0,'X_I','Color','w','FontSize',9);
text(0,L+0.02,0,'Y_I','Color','w','FontSize',9);
text(0,0,L+0.02,'Z_I','Color','w','FontSize',9);

%% Satellite body corners — defined relative to GEOMETRIC CENTER
Lx = sat.Lx; Ly = sat.Ly; Lz = sat.Lz;
corners_gc = [Lx/2*[-1  1  1 -1 -1  1  1 -1];
              Ly/2*[-1 -1  1  1 -1 -1  1  1];
              Lz/2*[-1 -1 -1 -1  1  1  1  1]];  % 3x8

%% Shift corners relative to CoM for correct rotation
CoM = sat.CoM;
corners = corners_gc - CoM;  % corners relative to CoM

faces = [1 2 3 4; 5 6 7 8;
         1 2 6 5; 3 4 8 7;
         2 3 7 6; 1 4 8 5];

face_colors = [0.5 0.5 0.8;
               1.0 0.4 0.0;  
               0.6 0.6 0.6;
               0.6 0.6 0.6;
               0.7 0.7 0.7;
               0.7 0.7 0.7];

hbody = patch('Vertices', corners', 'Faces', faces, ...
              'FaceVertexCData', face_colors, ...
              'FaceColor', 'flat', ...
              'FaceAlpha', 0.85, ...
              'EdgeColor', 'w', ...
              'LineWidth', 1.2);

%% CoM marker
hCoM = plot3(0, 0, 0, 'y+', 'MarkerSize', 12, 'LineWidth', 2.5);
text(0.02, 0, 0.02, 'CoM', 'Color', 'y', 'FontSize', 9, 'Tag', 'com_label');
h_com_label = findobj('Tag', 'com_label');

%% Geometric center marker
gc_offset = -CoM;  % geometric center relative to CoM
hGC = plot3(gc_offset(1), gc_offset(2), gc_offset(3), 'w+', 'MarkerSize', 8, 'LineWidth', 1.5);
hGC_label = text(gc_offset(1)+0.02, gc_offset(2), gc_offset(3)+0.02,'GC', 'Color', 'w', 'FontSize', 8);

%% Body axes (at CoM = origin)
hx = quiver3(0,0,0,L,0,0,'r-','LineWidth',2.5,'MaxHeadSize',0.3,'AutoScale','off');
hy = quiver3(0,0,0,0,L,0,'g-','LineWidth',2.5,'MaxHeadSize',0.3,'AutoScale','off');
hz = quiver3(0,0,0,0,0,L,'b-','LineWidth',2.5,'MaxHeadSize',0.3,'AutoScale','off');
htx = text(L+0.02,0,0,'X_B','Color','r','FontSize',9,'FontWeight','bold');
hty = text(0,L+0.02,0,'Y_B','Color','g','FontSize',9,'FontWeight','bold');
htz = text(0,0,L+0.02,'Z_B','Color','b','FontSize',9,'FontWeight','bold');

%% Nadir arrow
h_nadir     = quiver3(0,0,0,0,0,-L,'c-','LineWidth',3,'MaxHeadSize',0.4,'AutoScale','off');
h_nadir_txt = text(0,0,-L-0.05,'Earth','Color','c','FontSize',10,'FontWeight','bold');

%% CoM offset display
com_str = sprintf('CoM offset: [%.1f, %.1f, %.1f] mm', CoM*1000);
text(-0.65, -0.65, 0.35, com_str, 'Color','y','FontSize',9);

%% Status text
h_align = text(-0.65,-0.65, 0.65,'','Color','w','FontSize',10,'FontWeight','bold');
htime   = text(-0.65,-0.65, 0.55,'','Color','w','FontSize',11);

%% Animation loop
indices = 1:skip:size(q_store,2);
for k = indices
    Ex=q_store(1,k); Ey=q_store(2,k);
    Ez=q_store(3,k); n =q_store(4,k);

    R = [(1-2*(Ey^2+Ez^2)),  2*(Ex*Ey-n*Ez),   2*(Ex*Ez+n*Ey);
          2*(Ex*Ey+n*Ez),   (1-2*(Ex^2+Ez^2)),  2*(Ey*Ez-n*Ex);
          2*(Ex*Ez-n*Ey),    2*(Ey*Ez+n*Ex),   (1-2*(Ex^2+Ey^2))];

    %% Rotate body
    set(hbody, 'Vertices', (R*corners)');

    %% Rotate geometric center marker
    gc_rotated = R * gc_offset;
    set(hGC, 'XData', gc_rotated(1), 'YData', gc_rotated(2), 'ZData', gc_rotated(3));
    set(hGC_label, 'Position', gc_rotated + [0.02;0;0.02]);
    
    %% Rotate body axes
    xb=R(:,1); yb=R(:,2); zb=R(:,3);
    set(hx,'UData',L*xb(1),'VData',L*xb(2),'WData',L*xb(3));
    set(hy,'UData',L*yb(1),'VData',L*yb(2),'WData',L*yb(3));
    set(hz,'UData',L*zb(1),'VData',L*zb(2),'WData',L*zb(3));
    set(htx,'Position',(L+0.02)*xb);
    set(hty,'Position',(L+0.02)*yb);
    set(htz,'Position',(L+0.02)*zb);

    %% Nadir arrow
    r_ECI  = r_ECI_store(:,k);
    nadir_I = -r_ECI / norm(r_ECI);
    set(h_nadir, 'UData',L*nadir_I(1),'VData',L*nadir_I(2),'WData',L*nadir_I(3));
    set(h_nadir_txt, 'Position',(L+0.05)*nadir_I);

    %% Alignment check — Z_body toward nadir
    zb_nadir_angle = rad2deg(acos(max(-1,min(1,dot(zb, nadir_I)))));
    if zb_nadir_angle < 10
        set(h_align,'Color',[0 1 0],'String',sprintf('Z→Earth: %.1f° ✓',zb_nadir_angle));
    elseif zb_nadir_angle < 30
        set(h_align,'Color',[1 1 0],'String',sprintf('Z→Earth: %.1f° ~',zb_nadir_angle));
    else
        set(h_align,'Color',[1 0.3 0.3],'String',sprintf('Z→Earth: %.1f° ✗',zb_nadir_angle));
    end

    set(htime,'String',sprintf('t = %.1f s', tspan(k)));
    drawnow limitrate;
end
end