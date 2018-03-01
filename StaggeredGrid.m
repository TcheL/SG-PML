%% acoustic wave staggered grid

figure;
    axis equal;
    hold on;

    plot(0.5*ones(1,4),0:3,'r^','Markersize',25,'MarkerFacecolor','r');
    plot(0:3,2.5*ones(1,4),'gs','Markersize',25,'MarkerFacecolor','g');
    plot(0:3,3*ones(1,4),'ko','Markersize',25,'MarkerFacecolor','k');
    legend('v_x, 1/\rho','v_z, 1/\rho','P, \rho v_P^2','Location','EastOutside');
    plot(1.5*ones(1,4),0:3,'r^','Markersize',25,'MarkerFacecolor','r');
    plot(2.5*ones(1,4),0:3,'r^','Markersize',25,'MarkerFacecolor','r');
    plot(0:3,1.5*ones(1,4),'gs','Markersize',25,'MarkerFacecolor','g');
    plot(0:3,0.5*ones(1,4),'gs','Markersize',25,'MarkerFacecolor','g');
    plot(0:3,2*ones(1,4),'ko','Markersize',25,'MarkerFacecolor','k');
    plot(0:3,1*ones(1,4),'ko','Markersize',25,'MarkerFacecolor','k');
    plot(0:3,0*ones(1,4),'ko','Markersize',25,'MarkerFacecolor','k');

    plot([1 1],[0 3]);
    plot([2 2],[0 3]);
    plot([0 3],[1 1]);
    plot([0 3],[2 2]);
    plot([0 3],[0 0]);
    plot([0 3],[3 3]);
    plot([0 0],[0 3]);
    plot([3 3],[0 3]);

    text(-0.1,3.3,'i-1','Fontsize',25);
    text(1,3.3,'i','Fontsize',25);
    text(1.9,3.3,'i+1','Fontsize',25);
    text(-0.4,3,'j-1','Fontsize',25);
    text(-0.25,2,'j','Fontsize',25);
    text(-0.4,1,'j+1','Fontsize',25);
    text(1.5,3.3,'i''','Fontsize',25);
    text(-0.25,1.5,'j''','Fontsize',25);

    hold off;

%% elastic wave stagger grid

figure;
    axis equal;
    hold on;

    plot(0:3,3*ones(1,4),'rs','Markersize',25,'MarkerEdgecolor','r','LineWidth',5);
    plot(0.5:2.5,2.5*ones(1,3),'rs','Markersize',25,'MarkerFacecolor','r','LineWidth',5);
    plot(0.5*ones(1,4),0:3,'bo','Markersize',25,'MarkerEdgecolor','b','LineWidth',5);
    plot(0:3,2.5*ones(1,4),'bo','Markersize',25,'MarkerFacecolor','b','LineWidth',5);
    legend('v_x, 1/\rho','v_z, 1/\rho','\tau_{xx}, \tau_{zz}, (\lambda+2\mu), \lambda','\tau_{xz}, \mu','Location','EastOutside');
    plot(0:3,2*ones(1,4),'rs','Markersize',25,'MarkerEdgecolor','r','LineWidth',5);
    plot(0:3,1*ones(1,4),'rs','Markersize',25,'MarkerEdgecolor','r','LineWidth',5);
    plot(0:3,0*ones(1,4),'rs','Markersize',25,'MarkerEdgecolor','r','LineWidth',5);
    plot(0.5:2.5,1.5*ones(1,3),'rs','Markersize',25,'MarkerFacecolor','r','LineWidth',5);
    plot(0.5:2.5,0.5*ones(1,3),'rs','Markersize',25,'MarkerFacecolor','r','LineWidth',5);
    plot(1.5*ones(1,4),0:3,'bo','Markersize',25,'MarkerEdgecolor','b','LineWidth',5);
    plot(2.5*ones(1,4),0:3,'bo','Markersize',25,'MarkerEdgecolor','b','LineWidth',5);
    plot(0:3,1.5*ones(1,4),'bo','Markersize',25,'MarkerFacecolor','b','LineWidth',5);
    plot(0:3,0.5*ones(1,4),'bo','Markersize',25,'MarkerFacecolor','b','LineWidth',5);

    plot([1 1],[0 3]);
    plot([2 2],[0 3]);
    plot([0 3],[1 1]);
    plot([0 3],[2 2]);
    plot([0 3],[0 0]);
    plot([0 3],[3 3]);
    plot([0 0],[0 3]);
    plot([3 3],[0 3]);

    text(-0.1,3.3,'i-1','Fontsize',25);
    text(1,3.3,'i','Fontsize',25);
    text(1.9,3.3,'i+1','Fontsize',25);
    text(-0.4,3,'j-1','Fontsize',25);
    text(-0.25,2,'j','Fontsize',25);
    text(-0.4,1,'j+1','Fontsize',25);
    text(1.5,3.3,'i''','Fontsize',25);
    text(-0.25,1.5,'j''','Fontsize',25);

    hold off;
