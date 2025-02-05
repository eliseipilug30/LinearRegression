%incarcare date
clear all
close all
clc

load('proj_fit_17.mat') 
U_id = id.X;
Y_id = id.Y;
U_val = val.X;
Y_val = val.Y;
grad_maxim = 40;
Nid = id.dims(1);
Nval = val.dims(1);
MSE_val = zeros(1,grad_maxim);
MSE_id = zeros(1,grad_maxim);
figure
mesh(U_id{1},U_id{2},Y_id);
title("Semnalul original"); xlabel("x1"); ylabel("x2"); zlabel("y");

% Antrenare aproximator de diverse grade
for m = 1:grad_maxim

aux = calcul_regresori(U_id{1}(1),U_id{2}(1),m);
dim_theta = length(aux);

phi = zeros(Nid^2,dim_theta);
index_phi = 1;
for i = 1:Nid
    for j = 1:Nid   
        phi(index_phi,:) = calcul_regresori(U_id{1}(i),U_id{2}(j),m);
        index_phi = index_phi + 1;
    end
end

Y_id_coloana = coloana(Y_id,Nid);

theta = phi \ Y_id_coloana;
%theta = pinv(phi) * Y_id_coloana;
Y_aprox_id = phi * theta;
Y_id_aprox_matrix = matrice(Y_aprox_id,Nid);


index_phi_val = 1;
phi_val = zeros(Nval^2,dim_theta);
for i = 1:Nval
    for j = 1:Nval   
        phi_val(index_phi_val,:) = calcul_regresori(U_val{1}(i),U_val{2}(j),m);
        index_phi_val = index_phi_val + 1;
    end
end

Y_aprox_val = phi_val * theta;
Y_aprox_val_matrix = matrice(Y_aprox_val,Nval);

MSE_id(m) = sum((Y_id-Y_id_aprox_matrix).^2, 'all')/Nid^2;
MSE_val(m) = sum((Y_val-Y_aprox_val_matrix).^2, 'all')/Nval^2;

end

%Calculare MSE minim
index_min_val = find(min(MSE_val)==MSE_val);
index_min_id = find(min(MSE_id)==MSE_id);
MSE_min_id = MSE_id(index_min_id);
MSE_min_val = MSE_val(index_min_val);
figure
plot(MSE_val); title(["MSE minim la validare= ",num2str(MSE_min_val),"la grad ",index_min_val]);
xlabel("puterea"); ylabel("MSE");
figure
plot(MSE_id); title(["MSE minim la identificare= ",num2str(MSE_min_id),"la grad ",index_min_id]);
xlabel("puterea"); ylabel("MSE");
m = index_min_val;


%Recalculare model validare&identificare pt. grad minim
aux = calcul_regresori(U_id{1}(1),U_id{2}(1),m);
dim_theta = length(aux);

phi = zeros(Nid^2,dim_theta);
index_phi = 1;
for i = 1:Nid
    for j = 1:Nid   
        phi(index_phi,:) = calcul_regresori(U_id{1}(i),U_id{2}(j),m);
        index_phi = index_phi + 1;
    end
end

Y_id_coloana = coloana(Y_id,Nid);
theta = phi \ Y_id_coloana;
Y_aprox_id = phi * theta;
Y_id_aprox_matrix = matrice(Y_aprox_id,Nid);

figure
mesh(U_id{1},U_id{2},Y_id_aprox_matrix);
title("Identificare model"); xlabel("x1"); ylabel("x2"); zlabel("y");

index_phi_val = 1;
phi_val = zeros(Nval^2,dim_theta);
for i = 1:Nval
    for j = 1:Nval   
        phi_val(index_phi_val,:) = calcul_regresori(U_val{1}(i),U_val{2}(j),m);
        index_phi_val = index_phi_val + 1;
    end
end

Y_aprox_val = phi_val * theta;
Y_aprox_val_matrix = matrice(Y_aprox_val,Nval);
figure
mesh(U_val{1},U_val{2},Y_aprox_val_matrix);
title("Validare model"); xlabel("x1"); ylabel("x2"); zlabel("y");

%Functie pentru calculul regresorilor
function regresori = calcul_regresori(x1,x2,power)

   index = 1;         
   for i = 0:power
      for j = 0:power
            if i+j<=power
            regresori(index) = x1^(i)*x2^(j);
            index = index+1;
            end
      end
   end
end

%Functie pentru transformarea matricelor de NxN in matrice coloana de N^2x1
function Y_coloana = coloana(Y,N)

    Y_coloana = zeros(N^2,1);
    for i = 1:N
        Y_coloana(1+(i-1)*N:i*N,1) = Y(1:N,i);
    end

end

%Functie pentru transformarea matricelor coloana de N^2x1 in matrice de NxN
function Y_matrice = matrice(Y,N)

    Y_matrice = zeros(N,N);
    for i = 1:N
        Y_matrice(1:N,i) = Y(1+(i-1)*N:i*N,1);
    end

end