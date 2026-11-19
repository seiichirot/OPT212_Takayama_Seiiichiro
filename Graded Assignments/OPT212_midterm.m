%% 
% $$g(n_{eff})=\frac{(\beta^2-\alpha^2)sin(\beta d)-2\beta acos(\beta d)}{\beta^2+\alpha^2}=0$$
% 
% where $\beta=\sqrt{n_{core}^2-n_{eff}^2$, $\alpha=\sqrt{n_{eff}^2-n_{clad}^2$.
% 
% d is in units of wavelengths.

n_clad = 1.45;
n_core = 1.50;
n_eff = linspace(1.45,1.5,901)';
d = 5; %in um

beta = sqrt(n_core^2 - n_eff.^2);
alpha = sqrt(n_eff.^2 - n_clad^2);

g = ((beta.^2-alpha.^2).*sin(beta*d)-2*beta.*alpha.*cos(beta*d))./(beta.^2+alpha.^2);
figure
plot(n_eff,g)
yline(0,'--')
%% 
% Find zeros:


zeros_found = [];

for i = 1:length(g)-1
    if g(i)*g(i+1) < 0
        zeros_found = [zeros_found; n_eff(i)];
    end
end

zeros_found
%% 
% Use Newton's method to find index:

b = @(n) sqrt(n_core^2 - n.^2);
a = @(n) sqrt(n.^2 - n_clad^2);

gfun = @(n) ((b(n).^2-a(n).^2).*sin(b(n)*d) - 2*b(n).*a(n).*cos(b(n)*d))./ (b(n).^2+a(n).^2);

h = 1e-6;
tol = 1e-12; %this can be anything

for j = 1:length(zeros_found)
    x = zeros_found(j);  % Starting guess from brute force

    for i = 1:20
        % Approximate the derivative
        gp = (gfun(x+h)-gfun(x-h))/(2*h);

        % Newton's formula
        x_new = x - gfun(x)/gp;

        if abs(x_new-x) < tol
            break
        end

        x = x_new;
    end

    fprintf('Effective index = %.10f\n',x_new)
end
%% 
%