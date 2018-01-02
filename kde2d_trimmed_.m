% discrete cosine transform of initial data
a= dct2d(initial_data);
% now compute the optimal bandwidth^2
I=(0:n-1).^2; A2=a.^2;

t_star=fzero( @(t)(t-evolve(t)),[0,0.1]);

p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);


%%%HACK by bst (4/2011) for convenience
if exist('bw','var') && ~isempty(bw)
    t_x = (bw(1)/scaling(1))^2;
    t_y = (bw(2)/scaling(2))^2;
end
%%%END HACK

% smooth the discrete cosine transform of initial data using t_star
a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a;
% now apply the inverse discrete cosine transform
if nargout>1
    numel(a_t);
    prod(scaling);
    density=idct2d(a_t)*(numel(a_t)/prod(scaling));
    [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
end
bandwidth=sqrt([t_x,t_y]).*scaling;
