function result = do_OPLS_DA(X,y)
T = [];
P = [];
W = [];

w_transpose = y.' * X;
w = w_transpose.';
w = w/norm(w);
% for i = 1:3
t = X * w / (w.' * w);
c_transpose = t.' * y / (t.' * t);
c = c_transpose.';
u = y * c / (c.' * c);
p_transpose = t.' * X / (t.' * t);
p = p_transpose.';
% w_ortho = p - (w.' * p / (w.' * w)) * w;
% w_ortho = w_ortho / norm(w_ortho);
% t_ortho = X * w_ortho / (w_ortho.' * w_ortho);
% p_ortho_transpose = t_ortho.' * X / (t_ortho.' * t_ortho);
% p_ortho = p_ortho_transpose.';
% E_O_PLS = X - t_ortho * p_ortho.';
% X = E_O_PLS;
% T_ortho = [T_ortho, t_ortho];
% P_ortho = [P_ortho, p_ortho];
% W_ortho = [W_ortho, w_ortho];
% end
% 
% X_ortho = T_ortho * P_ortho.';
% [coeff, score] = pca(X_ortho);
% T_pca_ortho = score;
% P_pca_ortho = coeff;
% X_ortho = T_pca_ortho * P_pca_ortho.';

result = {t,p};
end
