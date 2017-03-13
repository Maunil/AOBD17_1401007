clear;

clc;

arr= 1:5:100;

for it=1:1:length(arr)
    row = 3000;
    col = 200;
    no_col = arr(it);
    x = randn(row,col);
    [U,S,V] = svd(x,0);

    B = zeros(col+no_col,no_col);
    A = randn(row,no_col);

    sum = col+no_col;
    for i =1:1:no_col
        V = [V; zeros(1,col)];
        B(sum,i)=1;
        sum = sum-1;
    end
    
    [m_u, n_u] = size(U);
    [m_s, n_s] = size(S);
    [m_v, n_v] = size(V);
    [m_a, n_a] = size(A);
    [m_b, n_b] = size(B);

    m = m_u;    
    n = m_v;    
    r = n_u;    
    c = n_a;    

    [Q,R] = qr([U,A]);

    %Q,R

    [m_q, n_q] = size(Q);
    [m_r, n_r] = size(R);


    if (n_q == r)
        P = [];
        Ra = [];
        dimRa = 0;
    else
        P = Q(:, r+1:n_q);
        Ra = R(r+1:m_r,  r+1:n_r);              
        dimRa = m_r - r;
    end
    %M = R(1:r, r+1:r+c);

    M = U'*A;

    [Q,R] = qr([V,B]);

    [m_q, n_q] = size(Q);
    [m_r, n_r] = size(R);

    if (n_q == r)
        Q = [];
        Rb = [];
        dimRb = 0;
    else
        Q = Q(:, r+1:n_q); 
        Rb = R(r+1:m_r, r+1:n_r);
        dimRb = m_r - r;
    end

    N = V'*B;

    Saug = [M; Ra] * [N; Rb]';
    S = [S, zeros(r, dimRb); zeros(dimRa, r), zeros(dimRa, dimRb)] + Saug;
    [U2,S2,V2] = svd(S,0);
    U = [U,P]*U2;
    S = S2;
    V = [V,Q]*V2;


    %ERROR calculation 
    new_matrix = U*S*V';
    new_addition = A*B';
    n_orm = norm(new_matrix,'fro');

    new_matrix = new_matrix(:,end-no_col+1:end);
    new_addition = new_addition(:,  end-no_col+1:end);

    D = abs(new_matrix-new_addition);

    [r1,c1] = size(D);
    sum = 0;

    for i=1:1:r1
        for j=1:1:c1
            sum = sum + D(i,j);
        end
    end
    error(it) = sum/(r1*c1);
    error(it) = error(it)/n_orm;
end

figure;
plot(arr,error,'-*m','linewidth',1.5);hold on 


