% wave plotting program

clear all;
if(0)
    % for TDFDAWFD2DSG.f90 and TDFDEWFD2DSG.f90
    for i = 1:5:500
        file = sprintf('./data/Snapshot_%4.4d.dat', i);
        data = importdata(file);

        pcolor(data);
        axis image;
        axis ij;
        shading interp;
        xlabel('x');
        ylabel('z');
        colormap jet;
        colorbar;
        drawnow;
%         WriteGif('plot.gif', gcf, 0.1);
    end
else
    % for TDFDEWFD2DSG.cpp and TDFDEWFD2DSG.cu
    for i = 0:1:9
        idx = 0 + 100*i;
        file = sprintf('./data/P%05d.bin', idx);
        fid = fopen(file, 'rb');
        dim = fread(fid, 2, 'int');
        data = fread(fid, [dim(1), dim(2)], 'float');
        fclose(fid);

        pcolor(data');
        axis image;
        axis ij;
        shading interp;
        xlabel('x');
        ylabel('z');
        colormap jet;
        colorbar;
        drawnow;
%         WriteGif('plot.gif', gcf, 0.1);
    end
end


