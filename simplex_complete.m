%------------------------------------

% NOMBRES:
%   - JESSENIA PIZA LONDOÑO.
%   - PAULA LORENA LOPEZ ROMERO.

%------------------------------------
function[X] = simplex_complete(minomax, A, b, c, verbose)
    %INICIO FASE 1
    if(minomax == 'max')
        c = c.*(-1);
    end
    sz = size(A);
    c_B = ones(1, sz(1,1));
    c_N = zeros(1, sz(1,2));
    B = eye(sz(1,1));
    A_new = [A B];
    szA = size(A_new);
    I_N = 1:sz(1,2);
    I_B = sz(1,2) + 1:szA(1,2);
    c_N_bar = c_N - c_B*B^(-1)*A;
    art = [];
    if (verbose)
        disp('Comienza la fase I')
        disp ('-------------------------------------------')
    end
    [X, I_B, I_N] = simplex_aux(I_N, I_B, c_N_bar, B, A, c_N, c_B, b, verbose);
    for i=sz(1,2) + 1:szA(1,2)
        art = [art X(i,:)];
    end
    if (max(art) ~= 0)
        for i=sz(1,2) + 1:szA(1,2)
            for q = 1:sz(1,2)
                if(A_new(:,i) == A(:,q))
                    Xbef = X(i,:);
                    X(i,:) = X(q,:);
                    X(q,:) = Xbef;
                end
            end
        end
    end
    art = [];
    for i=sz(1,2) + 1:szA(1,2)
        art = [art X(i,:)];
    end
    B = [];
    N = [];
    c_B = [];
    c_N = [];
    if (max(art) ~= 0)
        X = [];
        disp('Problema no tiene solución factible');
        return;
    else
        if (verbose)
            disp('Solución básica factible:')
            disp(X)
            disp('Termina la fase I')
            disp ('-------------------------------------------')
            disp ('-------------------------------------------')
            disp('Comienza la fase II')
            disp ('-------------------------------------------')
        end
        while (max(I_N) > sz(1,2))
            idx = I_N == max(I_N);
            I_N(:, idx) = [];
        end
        for q=1:sz(1,1)
            if (I_B(:,q) > sz(1,2))
                I_Bbef = I_B(:,q);
                I_B(:,q) = I_N(:,q);
                I_N(:,q) = I_Bbef;
            end
        end
        while (max(I_N) > sz(1,2))
            idx = I_N == max(I_N);
            I_N(:, idx) = [];
        end
        %INICIO FASE 2
        for q = I_B
           B = [B A(:, q)];
           c_B = [c_B c(:, q)];
        end
        for q = I_N
          N = [N A(:, q)];
          c_N = [c_N c(:, q)]; 
        end
        c_N_bar = c_N - c_B*B^(-1)*N;
        X = simplex_aux(I_N, I_B, c_N_bar, B, N, c_N, c_B, b, verbose);
%         disp('Solución Óptima:')
%         disp(X)
%         disp('Termina la fase II')
%         disp ('-------------------------------------------')
    end
end

function [X, I_B, I_N] = simplex_aux(I_N, I_B, c_N_bar, B, N, c_N, c_B, b, verbose)
    szB = size(B);
    rango = szB(1,1);
    b_bar = B^(-1)*b;
    z_0bef = 0;
    for i=1:rango
            X(I_B(:,i), :) = b_bar(i,:);
    end
    for i=I_N
        X(i,:) = 0;
    end
    j = 1;
    while (min(c_N_bar) < 0)
        c_N_bar = c_N - c_B*B^(-1)*N;
        b_bar = B^(-1)*b;
        z_0 = c_B*b_bar;
        z = z_0 + z_0bef;
        for i=1:rango
            X(I_B(:,i), :) = b_bar(i,:);
        end
        for i=I_N
            X(i,:) = 0;
        end
        if (min(c_N_bar) >= 0) 
            return;
        end
        c_in = find(c_N_bar == min(c_N_bar));
        c_in = c_in(:, 1);
        a_in = N (:, c_in); %Candidata entrar
        Y_k = B^(-1)*a_in;
        if(max(Y_k) <= 0)
            X = [];
            disp('El problema no tiene óptimo finito');
            break
        else
            X_k = [];
            for i = 1:length(Y_k)
                if (Y_k(i,:) > 0)
                   X_k = [X_k b_bar(i,:)/Y_k(i,:)];
                else
                    X_k = [X_k NaN];
                end
            end
            c_out = find(X_k == min(X_k));
            if(length(c_out)> 1)
                c_out = c_out(:, 1);
            end
            a_out = B(:, c_out); %Candidata salir
            I_Bbef = I_B(:, c_out);
            I_Bant = I_B;
            cand_ent = I_N(:, c_in);
            I_B(:, c_out) = I_N(:, c_in);
            I_N(:, c_in) = I_Bbef;
            B(:, c_out) = a_in;
            N (:, c_in) = a_out;
            c_Nbef = c_N(:, c_in);
            c_Bbef = c_B(:, c_out);
            c_N(:, c_in) = c_Bbef;
            c_B(:, c_out) = c_Nbef;
            z_0bef = z_0;
            if (verbose)
                ite = ['Iteración ', num2str(j), ':'];
                disp(ite);
                disp('La solución hasta ahora es: ')
                disp(X)
                B_ant = ['Indices de la base anterior: ', num2str(I_Bant)];
                disp(B_ant)
                dir = ['La dirección de movimiento: ', 'Entra ', num2str(cand_ent), ' y sale ', num2str(I_Bbef)];
                disp(dir)
                Base = ['Indices de las variables básicas nueva: ', num2str(I_B)];
                disp(Base)
                paso = ['El tamaño del paso: ', num2str(z)];
                disp(paso);
                disp ('-------------------------------------------')
            end
            j = j+1;
        end
    end
end