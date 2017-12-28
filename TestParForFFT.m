function TestParForFFT
clc;

N = 256;
disp(1)
f = tic;
ti = zeros(30,1);
for i=1:30
    t= tic;
    a = rand(N,N,N);
    y = fftn(a);
    ti(i) = toc(t);
end
toc(f);
fprintf('No Par No GPU average %.3f\n\n',mean(ti))



disp(2)
f=tic;
ti = zeros(30,1);
for i=1:30
    t= tic;
    a = rand(N,N,N,'gpuArray');
    y = gather(fftn(a));
    clear a
    ti(i) = toc(t);
end
toc(f);
fprintf('GPU average %.3f\n\n',mean(ti))




poolobj = gcp('nocreate'); % If pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(16);
end

disp(3)
f=tic;
ti = zeros(30,1);
parfor i=1:30
    t= tic;
    a = rand(N,N,N);
    y = fftn(a);
    ti(i) = toc(t);
end
toc(f);
fprintf('Par 16 average %.3f\n\n',mean(ti))

delete(poolobj);




poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(8);
end
disp(4)
f=tic;
ti = zeros(30,1);
parfor i=1:30
    t= tic;
    a = rand(N,N,N);
    y = fftn(a);
    ti(i) = toc(t);
end
toc(f);
fprintf('Par 8 average %.3f\n\n',mean(ti))
delete(poolobj);




poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(4);
end

disp(5)
f=tic;
ti = zeros(30,1);
parfor i=1:30
    t= tic;
    a = rand(N,N,N);
    y = fftn(a);
    ti(i) = toc(t);
end
toc(f);
fprintf('Par 4 average %.3f\n\n',mean(ti))
delete(poolobj);