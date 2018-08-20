% wave plotting program for TDFDEWFD2DSG.cpp and TDFDEWFD2DSG.cu

for i = 0:1:9
  idx = 0 + 100*i;
  file = sprintf('P%05d.bin', idx);
  fid = fopen(file, 'rb');
  dim = fread(fid, 2, 'int');
  data = fread(fid, [dim(1), dim(2)], 'float');
  fclose(fid);

  pcolor(data');
  axis image;
  shading interp;
  xlabel('x');
  ylabel('z');
  colorbar;
  drawnow;
end


