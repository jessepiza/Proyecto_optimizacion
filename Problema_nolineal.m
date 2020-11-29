%------------------------------------

% NOMBRES:
%   - JESSENIA PIZA LONDOÑO.
%   - PAULA LORENA LOPEZ ROMERO.

%------------------------------------

clear all
syms x_11 x_12 x_13 x_14 x_15 x_21 x_22 x_23 x_24 x_25 x_31 x_32 x_33 x_34 x_35 x_41 x_42 x_43 x_44 x_45 x_51 x_52 x_53 x_54 x_55;
% Agregamos la función objetivo del problema no lineal.
f = @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)... 
        2000*exp(x_11) + 1400*exp(x_12) + 1300*exp(x_13) + 1800*exp(x_14) + 1900*exp(x_15) + ...
        700*exp(x_21) + 490*exp(x_22) + 456*exp(x_23) + 660*exp(x_24)+ 668*exp(x_25) +...
        5000*exp(x_31) + 3500*exp(x_32)+ 3250*exp(x_33)+ 4700*exp(x_34)+ 4750*exp(x_35) +...
        2200*exp(x_41) + 1540*exp(x_42)+ 1430*exp(x_43) + 2068*exp(x_44)+ 2090*exp(x_45) +...
        24572*exp(x_51)+ 17200*exp(x_52)+ 11180*exp(x_53)+ 16168*exp(x_54)+ 16340*exp(x_55);

% Inicializamos x_k+1 y x_k donde x_k+1 lo empezamos como matriz de
% infinitos.
x_k = [log(3); 0; 0; 0; 0;
       log(2); 0; 0; 0; 0;
       log(4); 0; 0; 0; 0;
       log(3); 0; 0; 0; 0;
       log(4); 0; 0; 0; 0];
x_k1 = inf(25, 1);

% Copiamos las restricciones del problema
rest = {@(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -114*x_11 - 24.7*x_12 - 26.6*x_13 - 68.78*x_14 - 60.8*x_15 + 300;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -114*x_21 - 24.7*x_22 - 26.6*x_23 - 68.78*x_24 - 60.8*x_25 + 150;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -114*x_31 - 24.7*x_32 - 26.6*x_33 - 68.78*x_34 - 60.8*x_35 + 140;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -114*x_41 - 24.7*x_42 - 26.6*x_43 - 68.78*x_44 - 60.8*x_45 + 180;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -114*x_51 - 24.7*x_52 - 26.6*x_53 - 68.78*x_54 - 60.8*x_55 + 190;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -3650*x_11 - 3550*x_12 - 15000*x_13 - 21500*x_14 - 22500*x_15 + 159600;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -23400*x_21 - 14040*x_22 - 54600*x_23 - 70200*x_24 - 143520*x_25 + 34713;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -27300*x_31 - 16380*x_32 - 63700*x_33 - 81900*x_34 - 167440*x_35 + 39900;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -50500*x_41 - 26250*x_42 - 87500*x_43 - 161250*x_44 - 157500*x_45 + 99750;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)...
        -50400*x_51 - 25200*x_52 - 84000*x_53 - 125000*x_54 - 151200*x_55 + 95760;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_11 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_12 - log(4);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_13 - log(5);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_14 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_15 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_21 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_22 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_23 - log(4);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_24 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_25 - log(4);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_31 - log(4);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_32 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_33 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_34 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_35 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_41 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_42 - log(5);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_43 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_44 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_45 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_51 - log(4);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_52 - log(3);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_53 - log(5);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_54 - log(2);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)x_55 - log(5);
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_11;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_12;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_13;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_14;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_15;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_21;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_22;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_23;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_24;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_25;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_31;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_32;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_33;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_34;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_35;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_41;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_42;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_43;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_44;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_45;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_51;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_52;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_53;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_54;
        @(x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55)-x_55;};

% Creamos nuestra matriz A de coeficientes de las restricciones anteriores.
A = [-114, -24.7, -26.6, -68.78, -60.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, -114, -24.7, -26.6, -68.78, -60.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -114, -24.7, -26.6, -68.78, -60.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -114, -24.7, -26.6, -68.78, -60.8, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -114, -24.7, -26.6, -68.78, -60.8;
    -3650, -3550,-15000, -21500, -22500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, -23400, -14040, -54600, -70200, -143520, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -27300, -16380, -63700, -81900, -167440, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50500, -26250, -87500, -161250, -157500, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50400, -25200, -84000, -125000, -151200;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];

% Definimos nuestro epsilon, es decir, la tolerancia del problema, y
% agregamos un booleano tal que defina cuándo se puede calcular w con la
% matriz A y cuándo no.
tol = 0.0001;
calcular_w = true;

% Comenzamos el algoritmo del Método de gradiente proyectado.
while norm(x_k1 - x_k) > tol
    ant_x = inf;
    if (x_k1 ~= inf)
        ant_x = x_k1;
    end
    % Gradiente de la función objetivo.
    g = [diff(f, x_11); diff(f, x_12); diff(f, x_13); diff(f, x_14); diff(f, x_15);
        diff(f, x_21); diff(f, x_22); diff(f, x_23); diff(f, x_24); diff(f, x_25);
        diff(f, x_31); diff(f, x_32); diff(f, x_33); diff(f, x_34); diff(f, x_35);
        diff(f, x_41); diff(f, x_42); diff(f, x_43); diff(f, x_44); diff(f, x_45);
        diff(f, x_51); diff(f, x_52); diff(f, x_53); diff(f, x_54); diff(f, x_55)];
    % Gradiente de la función evaluada en el punto x_k.
    g_k = subs(g, {x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55}, ...
        {x_k(1,1), x_k(2,1), x_k(3,1), x_k(4,1), x_k(5,1), x_k(6,1), x_k(7,1), x_k(8,1), x_k(9,1), x_k(10,1), x_k(11,1), x_k(12,1), x_k(13,1), x_k(14,1), x_k(15,1), x_k(16,1), x_k(17,1), x_k(18,1), x_k(19,1), x_k(20,1), x_k(21,1), x_k(22,1), x_k(23,1), x_k(24,1), x_k(25,1)});
    
    % Evaluamos las restricción en el punto x_k.
    rest_k = zeros(size(rest,1), 1);
    for i=1:size(rest,1)
        rest_k(i) = rest{i}(x_k(1,1), x_k(2,1), x_k(3,1), x_k(4,1), x_k(5,1), x_k(6,1), x_k(7,1), x_k(8,1), x_k(9,1), x_k(10,1), x_k(11,1), x_k(12,1), x_k(13,1), x_k(14,1), x_k(15,1), x_k(16,1), x_k(17,1), x_k(18,1), x_k(19,1), x_k(20,1), x_k(21,1), x_k(22,1), x_k(23,1), x_k(24,1), x_k(25,1));
    end
    
    % Observarmos las restricciones activas y calculamos el espacio de
    % trabajo.
    if calcular_w
        w= find(rest_k == 0)';
    end
    
    % Obtenemos la matriz A_k que representa las filas de A que están en el
    % espacio de trabajo.
    A_k = zeros(size(w, 2), size(A, 2));
    for i=1:length(w)
        A_k(i, :) = A(w(i), :);
    end

    % Calculamos P y d_k la dirección.
    P = eye(size(A,2)) - A_k'*(A_k*A_k')^-1*A_k;
    d_k = -P*g_k;

    % Evaluamos si la dirección es 0 o no.
    if (d_k == 0)
        % Calcular mu.
        mu = -(A_k*A_k')^-1*A_k*g_k;
        if(mu >= 0)
            x_k1 = x_k;
            break
        else
            mu_min = find(mu == min(mu));
            mu_min = mu_min(1,1);
            % Cambiamos el espacio de trabajo.
            w = setdiff(w, w(:, mu_min));
            calcular_w=false;
        end
    else
        syms a;
        aux = [];
        alpha1 = [];
        % Calculamos x_k+1 = x_k + a*d_k)
        x_k1 = x_k + a*d_k; 
        for i=1:size(rest, 1)
            aux = [aux; subs(rest(i, :), {x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55}, ...
                {x_k1(1,1), x_k1(2,1), x_k1(3,1), x_k1(4,1), x_k1(5,1), x_k1(6,1), x_k1(7,1), x_k1(8,1), x_k1(9,1), x_k1(10,1), x_k1(11,1), x_k1(12,1), x_k1(13,1), x_k1(14,1), x_k1(15,1), x_k1(16,1), x_k1(17,1), x_k1(18,1), x_k1(19,1), x_k1(20,1), x_k1(21,1), x_k1(22,1), x_k1(23,1), x_k1(24,1), x_k1(25,1)})];
        end
        idx = find(aux == 0);
        if ~isempty(idx)
            for i=idx
                aux(i,:) = [];
            end
        end
        % Calculamos alpha_1
        for i=1:size(aux,1)
            alpha1 = [alpha1;
                      solve(aux(i,:), a)];
        end
        alpha1 = min(alpha1(alpha1 >= 0));
        % Función para hallar alpha2 con método aureo.
        func = @(a) f(x_k(1,1) + a*d_k(1,1),... 
                      x_k(2,1) + a*d_k(2,1),... 
                      x_k(3,1) + a*d_k(3,1),... 
                      x_k(4,1) + a*d_k(4,1),... 
                      x_k(5,1) + a*d_k(5,1),...
                      x_k(6,1) + a*d_k(6,1),...
                      x_k(7,1) + a*d_k(7,1),...
                      x_k(8,1) + a*d_k(8,1),...
                      x_k(9,1) + a*d_k(9,1),...
                      x_k(10,1) + a*d_k(10,1),...
                      x_k(11,1) + a*d_k(11,1),...
                      x_k(12,1) + a*d_k(12,1),...
                      x_k(13,1) + a*d_k(13,1),...
                      x_k(14,1) + a*d_k(14,1),...
                      x_k(15,1) + a*d_k(15,1),...
                      x_k(16,1) + a*d_k(16,1),...
                      x_k(17,1) + a*d_k(17,1),...
                      x_k(18,1) + a*d_k(18,1),...
                      x_k(19,1) + a*d_k(19,1),...
                      x_k(20,1) + a*d_k(20,1),...
                      x_k(21,1) + a*d_k(21,1),...
                      x_k(22,1) + a*d_k(22,1),...
                      x_k(23,1) + a*d_k(23,1),...
                      x_k(24,1) + a*d_k(24,1),...
                      x_k(25,1) + a*d_k(25,1));
        alpha2 = MetodoAureo(tol, func, 0, alpha1, 0);
        % Reemplazamos y actualizamos x_k y x_k+1
        if (ant_x ~= inf)
            x_k = ant_x;
        end
        x_k1 = x_k + alpha2*d_k;
    end
end
% Reemplazo de los valores de x_k+1 en la función objetivo para hallar el
% minimo de la función objetivo.
x_k1 = round(x_k1, 4);
mini_f = subs(f, {x_11, x_12, x_13, x_14, x_15, x_21, x_22, x_23, x_24, x_25, x_31, x_32, x_33, x_34, x_35, x_41, x_42, x_43, x_44, x_45, x_51, x_52, x_53, x_54, x_55}, ...
        {x_k1(1,1), x_k1(2,1), x_k1(3,1), x_k1(4,1), x_k1(5,1), x_k1(6,1), x_k1(7,1), x_k1(8,1), x_k1(9,1), x_k1(10,1), x_k1(11,1), x_k1(12,1), x_k1(13,1), x_k1(14,1), x_k1(15,1), x_k1(16,1), x_k1(17,1), x_k1(18,1), x_k1(19,1), x_k1(20,1), x_k1(21,1), x_k1(22,1), x_k1(23,1), x_k1(24,1), x_k1(25,1)});
vuelos = round(exp(x_k1), 4)';
disp('Solución óptima para el problema:');
disp(x_k1');
disp('---------------------------------');
disp('Costos mínimos de la aerolinea: ');
disp(round(mini_f, 4));
disp('---------------------------------');
disp('Vuelos por avión i en el trajecto j');
disp(vuelos);