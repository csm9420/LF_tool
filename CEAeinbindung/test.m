clear

for i = 1:10
%     p_BK(i) = 60E5;
    p_BK(i) = i * 10 * 1E5;
%     ROF(i) = 1 + 0.5*i;
    ROF(i) = 6;
    d_t(i) = 28 * 1E-3;
    [rho_star(i),a_star(i)] = cea(p_BK(i),ROF(i),d_t(i));
    m_dot(i) = rho_star(i) * a_star(i) * pi/4*d_t(i)^2;
end

plot(p_BK,m_dot)
% plot(ROF,m_dot)

fid = fopen('cea.out');
C = textscan(fid,'%s');