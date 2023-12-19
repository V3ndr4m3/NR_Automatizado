clc; clear all; close all

%Determina o número de barras e o tipo de cada barra
%1: Barra PQ, 2: Barra PV, 3: Barra de referência
barras = [3 1 1];

%Valores especificados para cada barra
esp = [1 0;
    -0.3 -0.5;
    -0.5 -0.2];

%Separa em rótulos
for i = 1:1:length(barras)
    if barras(i) == 3
        V(i) = esp(i,1);
        theta(i) = esp(i,2);
        P(i) = 0;
        Q(i) = 0;
    end
    if barras(i) == 1
        V(i) = 1;
        theta(i) = 0;
        P(i) = esp(i,1);
        Q(i) = esp(i,2);
    end
    if barras(i) == 2
        V(i) = esp(i,2);
        theta(i) = 0;
        P(i) = esp(i,1);
        Q(i) = 0;
    end
end

%Preparação para montagem da Jacobiana e suas submatrizes
%Define os indíces das barras onde se conhece P e/ou Q
P_barras = find(barras == 1 | barras == 2);
Q_barras = find(barras == 1);

%Define os indíces das barras onde não se conhece módulo e/ou ângulo da
%tensão
V_barras = find(barras == 1);
theta_barras = find(barras == 1 | barras == 2);

%Define o tamanho da jacobiana (linhas x colunas)
tam_J = [length(P_barras)+length(Q_barras) length(theta_barras)+length(V_barras)];

%Define o tamanho da submatriz H (linhas x colunas)
tam_H = [length(P_barras) length(theta_barras)];

%Define o tamanho da submatriz N (linhasxcolunas)
tam_N = [length(P_barras) length(V_barras)];

%Define o tamanho da submatriz M (linhas x colunas)
tam_M = [length(Q_barras) length(theta_barras)];

%Define o tamanho da submatriz L (linhasxcolunas)
tam_L = [length(Q_barras) length(V_barras)];


%Define a matriz admintância
%Impedâncias e admitâncias entre barras
z12 = 1.9+5.9i;
z13 = 1.3+4.8i;      
z23 = 0.5+6i;

y12 = z12';
y13 = z13';      
y23 = z23';

%Matriz de Admitância
Y = [y12+y13 -y12 -y13;
    -y12 y23+y12 -y23;
    -y13 -y23 y13+y23];
    
%Matriz de Condutância
G = real(Y);

%Matriz de Susceptância
B = imag(Y);



%Calcula as potências das barras conhecidas
P_calc = zeros(1,length(barras));
Q_calc = zeros(1,length(barras));
for i = 1:1:length(barras)
    if barras(i)==1 || barras(i)==2
        for j=1:1:length(barras)
            P_calc(i) = P_calc(i)+V(i)*V(j)*(G(i,j)*cos(theta(i)-theta(j))...
                +B(i,j)*sin(theta(i)-theta(j)));
        end
    end
    if barras(i)==1
        for j=1:1:length(barras)
            Q_calc(i) = Q_calc(i)+V(i)*V(j)*(G(i,j)*sin(theta(i)-theta(j))...
                -B(i,j)*cos(theta(i)-theta(j)));
        end
    end
end
%Verifica os mismatches de potência
deltaP = P-P_calc;
deltaQ = Q-Q_calc;

%Determina o vetor de mismatches de potências
a = 0;
for i = 1:1:length(barras)
    if barras(i) ~= 3
        a = a+1;
        vet_pot(a) = deltaP(i);
    end
end
for i = 1:1:length(barras)
    if barras(i) == 1
        a = a+1;
        vet_pot(a) = deltaQ(i);
    end
end
vet_pot = vet_pot';

erro = max(abs(vet_pot));

%Determinação da tolerância
tol = 0.0001;
it = 0;      %iteração 0

while erro>tol
    
    %Determina as submatrizes da Jacobiana
    
    %Submatriz H
    H = zeros(tam_H(1,1),tam_H(1,2));
    for i=1:1:tam_H(1,1)
        for j=1:1:tam_H(1,2)
            ind1 = P_barras(i);
            ind2 = theta_barras(j);
            for k = 1:1:length(barras)
                if ind1 == ind2
                    if k~=ind1
                        H(i,j) = H(i,j) - V(ind1)*V(k)*(G(ind1,k)*sin(theta(ind1)-theta(k))...
                            -B(ind1,k)*cos(theta(ind1)-theta(k)));
                    end
                end
                if ind1 ~= ind2
                    H(i,j) = V(ind1)*V(ind2)*(G(ind1,ind2)*sin(theta(ind1)-theta(ind2))...
                        -B(ind1,ind2)*cos(theta(ind1)-theta(ind2)));
                end
            end
        end
    end
    
    %     %Submatriz N
    
    
    N = zeros(tam_N(1,1),tam_N(1,2));
    for i=1:1:tam_N(1,1)
        for j=1:1:tam_N(1,2)
            ind1 = P_barras(i);
            ind2 = V_barras(j);
            for k = 1:1:length(barras)
                if ind1 == ind2
                    if k~=ind1
                        N(i,j) = N(i,j) + V(k)*(G(ind1,k)*cos(theta(ind1)-theta(k))...
                            +B(ind1,k)*sin(theta(ind1)-theta(k)));
                    end
                    if k == ind1
                        N(i,j) = N(i,j)+2*G(k,k)*V(k);
                    end
                end
                if ind1 ~= ind2
                    N(i,j) = V(ind1)*(G(ind1,ind2)*cos(theta(ind1)-theta(ind2))...
                        +B(ind1,ind2)*sin(theta(ind1)-theta(ind2)));
                end
            end
        end
    end
    
    %     %Submatriz M
    M = zeros(tam_M(1,1),tam_M(1,2));
    for i=1:1:tam_M(1,1)
        for j=1:1:tam_M(1,2)
            ind1 = Q_barras(i);
            ind2 = theta_barras(j);
            for k = 1:1:length(barras)
                if ind1 == ind2
                    if k~=ind1
                        M(i,j) = M(i,j) + V(ind1)*V(k)*(G(ind1,k)*cos(theta(ind1)-theta(k))...
                            +B(ind1,k)*sin(theta(ind1)-theta(k)));
                    end
                end
                if ind1 ~= ind2
                    M(i,j) = -V(ind1)*V(k)*(G(ind1,ind2)*cos(theta(ind1)-theta(ind2))...
                        +B(ind1,ind2)*sin(theta(ind1)-theta(ind2)));
                end
            end
        end
    end
    
    %     %Submatriz L
    
    L = zeros(tam_L(1,1),tam_L(1,2));
    for i=1:1:tam_L(1,1)
        for j=1:1:tam_L(1,2)
            ind1 = Q_barras(i);
            ind2 = V_barras(j);
            for k = 1:1:length(barras)
                if ind1 == ind2
                    if k~=ind1
                        L(i,j) = L(i,j) + V(k)*(G(ind1,k)*sin(theta(ind1)-theta(k))...
                            -B(ind1,k)*cos(theta(ind1)-theta(k)));
                    end
                    if k == ind1
                        L(i,j) = L(i,j)-2*B(k,k)*V(k);
                    end
                end
                if ind1 ~= ind2
                    L(i,j) = V(ind1)*(G(ind1,ind2)*sin(theta(ind1)-theta(ind2))...
                        -B(ind1,ind2)*cos(theta(ind1)-theta(ind2)));
                end
            end
        end
    end
    
    %     %Determina a Jacobiana
    J = [H N; M L];
    
    %Calcula o inverso da Jacobiana
    J_inv = inv(J);
    
    
    %Calcula os valores de variação de módulo e ângulo das tensões não
    %especificadas
    vet_tensao = J_inv*vet_pot;
    
    %Atualiza as tensões
    theta(theta_barras) = theta(theta_barras)+vet_tensao(1:length(theta_barras))';
    V(V_barras) = V(V_barras)+vet_tensao(length(theta_barras)+1:end)';
    
    
    
    %atualiza iterações feitas
    it = it+1;
    
    %Calcula as potências das barras conhecidas
    P_calc = zeros(1,length(barras));
    Q_calc = zeros(1,length(barras));
    for i = 1:1:length(barras)
        if barras(i)==1 || barras(i)==2
            for j=1:1:length(barras)
                P_calc(i) = P_calc(i)+V(i)*V(j)*(G(i,j)*cos(theta(i)-theta(j))...
                    +B(i,j)*sin(theta(i)-theta(j)));
            end
        end
        if barras(i)==1
            for j=1:1:length(barras)
                Q_calc(i) = Q_calc(i)+V(i)*V(j)*(G(i,j)*sin(theta(i)-theta(j))...
                    -B(i,j)*cos(theta(i)-theta(j)));
            end
        end
    end
    
    %Verifica os mismatches de potência
    deltaP = P-P_calc;
    deltaQ = Q-Q_calc;
    
    %Determina o vetor de mismatches de potências
    a = 0;
    for i = 1:1:length(barras)
        if barras(i) ~= 3
            a = a+1;
            vet_pot(a) = deltaP(i);
        end
    end
    for i = 1:1:length(barras)
        if barras(i) == 1
            a = a+1;
            vet_pot(a) = deltaQ(i);
        end
    end
    
    erro = max(abs(vet_pot));
end


%Cálculo das Potências desconhecidas
P_calc = zeros(1,length(barras));
Q_calc = zeros(1,length(barras));
for i = 1:1:length(barras)
    for j=1:1:length(barras)
        P_calc(i) = P_calc(i)+V(i)*V(j)*(G(i,j)*cos(theta(i)-theta(j))...
            +B(i,j)*sin(theta(i)-theta(j)));
        Q_calc(i) = Q_calc(i)+V(i)*V(j)*(G(i,j)*sin(theta(i)-theta(j))...
            -B(i,j)*cos(theta(i)-theta(j)));
    end    
end

for i = 1:1:length(barras)
    if barras(i) == 1
        barra(1,i) = "PQ";
    end
    if barras(i) == 2
        barra(1,i) = "PV";
    end
    if barras(i) == 3
        barra(1,i) = "Vθ";
    end
end
numero_barras = 1:1:length(barras);
valores_finais = [numero_barras' barra' V' theta' P_calc' Q_calc'];
disp("    Valores Finais Calculados de Tensão e Potência nas Barras");
disp("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||");
disp('   Barra   Tipo       V           θ               P              Q');
disp(valores_finais);
disp("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||");