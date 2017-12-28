function TestGPU

clc;

N= 512;
tic;
a = rand(N,N,N,'gpuArray');
toc;
tic;
b = gather(a);
for i=1:N
    c= b(:,:,i);
end
toc;
tic;
for i=1:N
    c = gather(a(:,:,i));
end
toc;

b = a;
tic;
b = (abs(ifftshift(fftshift(fftn((a))),3)).^2)/(N^3);
toc;