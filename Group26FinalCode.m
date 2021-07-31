%Group 26 Final Truss

% Clear Command Window and Workspace
clear all, clc;

% All distances in metres
P = -9.81;
E = 1.0392e10;
T = 0.0015875; %1/16 inch
W = [0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.0175; 0.02; 0.02;
    0.02; 0.02; 0.02; 0.0175;];
A = W .* T; %Area = Width times Thickness

L = [0.05; 0.05; 0.075; 0.12; 0.055; 0.0875; 0.095; 0.0807; 0.038891;  
    5*sqrt(2)/100; 0.0625; 0.0625; 0.0781025; 0.0781025; 0.038891;];%Length Array

theta = [90; 0; 0; 0; 0; 0; 0; -16.1892; -45; -45; 53.1301; -53.1391;
         39.8056; -39.8056; 45]; %theta array

%k coefficients
k = E.*A./L;

%k matrices
C = cosd(theta);
S = sind(theta);
CS = C.*S;


for vect = 1:15
k1 = k(1) * [C(1)^2, C(1)*S(1), -C(1)^2, -C(1)*S(1);
             S(1)^2, C(1)*S(1), -S(1)^2, -C(1)*S(1);
             -C(1)^2, -C(1)*S(1), C(1)^2, C(1)*S(1);
             -S(1)^2, -C(1)*S(1), S(1)^2, C(1)*S(1);];
end
 
k2 = k(2) * [C(2)^2, C(2)*S(2), -C(2)^2, -C(2)*S(2);
             S(2)^2, C(2)*S(2), -S(2)^2, -C(2)*S(2);
             -C(2)^2, -C(2)*S(2), C(2)^2, C(2)*S(2);
             -S(2)^2, -C(2)*S(2), S(2)^2, C(2)*S(2);];
         
k3 = k(3) * [C(3)^2, C(3)*S(3), -C(3)^2, -C(3)*S(3);
             S(3)^2, C(3)*S(3), -S(3)^2, -C(3)*S(3);
             -C(3)^2, -C(3)*S(3), C(3)^2, C(3)*S(3);
             -S(3)^2, -C(3)*S(3), S(3)^2, C(3)*S(3);];
         
k4 = k(4) * [C(4)^2, C(4)*S(4), -C(4)^2, -C(4)*S(4);
             S(4)^2, C(4)*S(4), -S(4)^2, -C(4)*S(4);
             -C(4)^2, -C(4)*S(4), C(4)^2, C(4)*S(4);
             -S(4)^2, -C(4)*S(4), S(4)^2, C(4)*S(4);];
         
k5 = k(5) * [C(5)^2, C(5)*S(5), -C(5)^2, -C(5)*S(5);
             S(5)^2, C(5)*S(5), -S(5)^2, -C(5)*S(5);
             -C(5)^2, -C(5)*S(5), C(5)^2, C(5)*S(5);
             -S(5)^2, -C(5)*S(5), S(5)^2, C(5)*S(5);];
         
k6 = k(6) * [C(6)^2, C(6)*S(6), -C(6)^2, -C(6)*S(6);
             S(6)^2, C(6)*S(6), -S(6)^2, -C(6)*S(6);
             -C(6)^2, -C(6)*S(6), C(6)^2, C(6)*S(6);
             -S(6)^2, -C(6)*S(6), S(6)^2, C(6)*S(6);];
         
k7 = k(7) * [C(7)^2, C(7)*S(7), -C(7)^2, -C(7)*S(7);
             S(7)^2, C(7)*S(7), -S(7)^2, -C(7)*S(7);
             -C(7)^2, -C(7)*S(7), C(7)^2, C(7)*S(7);
             -S(7)^2, -C(7)*S(7), S(7)^2, C(7)*S(7);];
     
k8 = k(8) * [C(8)^2, C(8)*S(8), -C(8)^2, -C(8)*S(8);
             S(8)^2, C(8)*S(8), -S(8)^2, -C(8)*S(8);
             -C(8)^2, -C(8)*S(8), C(8)^2, C(8)*S(8);
             -S(8)^2, -C(8)*S(8), S(8)^2, C(8)*S(8);];
         
k9 = k(9) * [C(9)^2, C(9)*S(9), -C(9)^2, -C(9)*S(9);
             S(9)^2, C(9)*S(9), -S(9)^2, -C(9)*S(9);
             -C(9)^2, -C(9)*S(9), C(9)^2, C(9)*S(9);
             -S(9)^2, -C(9)*S(9), S(9)^2, C(9)*S(9);];
         
k10 = k(10) * [C(10)^2, C(10)*S(10), -C(10)^2, -C(10)*S(10);
              S(10)^2, C(10)*S(10), -S(10)^2, -C(10)*S(10);
              -C(10)^2, -C(10)*S(10), C(10)^2, C(10)*S(10);
              -S(10)^2, -C(10)*S(10), S(10)^2, C(10)*S(10);];
         
k11 = k(11) * [C(11)^2, C(11)*S(11), -C(11)^2, -C(11)*S(11);
              S(11)^2, C(11)*S(11), -S(11)^2, -C(11)*S(11);
              -C(11)^2, -C(11)*S(11), C(11)^2, C(11)*S(11);
              -S(11)^2, -C(11)*S(11), S(11)^2, C(11)*S(11);];
         
k12 = k(12) * [C(12)^2, C(12)*S(12), -C(12)^2, -C(12)*S(12);
              S(12)^2, C(12)*S(12), -S(12)^2, -C(12)*S(12);
              -C(12)^2, -C(12)*S(12), C(12)^2, C(12)*S(12);
              -S(12)^2, -C(12)*S(12), S(12)^2, C(12)*S(12);];
         
k13 = k(13) * [C(13)^2, C(13)*S(13), -C(13)^2, -C(13)*S(13);
              S(13)^2, C(13)*S(13), -S(13)^2, -C(13)*S(13);
              -C(13)^2, -C(13)*S(13), C(13)^2, C(13)*S(13);
              -S(13)^2, -C(13)*S(13), S(13)^2, C(13)*S(13);];
         
k14 = k(14) * [C(14)^2, C(14)*S(14), -C(14)^2, -C(14)*S(14);
              S(14)^2, C(14)*S(14), -S(14)^2, -C(14)*S(14);
              -C(14)^2, -C(14)*S(14), C(14)^2, C(14)*S(14);
              -S(14)^2, -C(14)*S(14), S(14)^2, C(14)*S(14);];
         
k15 = k(15) * [C(15)^2, C(15)*S(15), -C(15)^2, -C(15)*S(15);
              S(15)^2, C(15)*S(15), -S(15)^2, -C(15)*S(15);
              -C(15)^2, -C(15)*S(15), C(15)^2, C(15)*S(15);
              -S(15)^2, -C(15)*S(15), S(15)^2, C(15)*S(15);];
          
%Assemble the total k matrix
K = zeros(18);

%Node 1
K(1:2,1:2) = k1(1:2,1:2) + k6(1:2,1:2) + k10(1:2,1:2);
K(1:2,3:4) = k1(1:2,3:4);
K(3:4,1:2) = k1(3:4,1:2);
K(1:2,5:6) = k10(1:2,3:4);
K(5:6,1:2) = k10(1:2,3:4);
K(13:14,1:2) = k6(1:2,3:4);
K(1:2,13:14) = k6(1:2,3:4);

%Node 2
K(3:4,3:4) = k1(3:4,3:4) + k2(1:2,1:2);
K(3:4,5:6) = k2(1:2,3:4);
K(5:6,3:4) = k2(3:4,1:2);

%Node 3
K(5:6,5:6) = k2(3:4,3:4) + k3(1:2,1:2) + k10(1:2,1:2) + k11(1:2,1:2);
K(5:6,7:8) = k3(1:2,3:4);
K(7:8,5:6) = k3(3:4,1:2);
K(13:14,5:6) = k11(3:4,1:2);
K(5:6,13:14) = k11(3:4,1:2);

%Node 4
K(7:8,7:8) = k3(3:4,3:4) + k4(1:2,1:2) + k12(1:2,1:2) + k13(1:2,1:2);
K(7:8,9:10) = k4(1:2,3:4);
K(9:10,7:8) = k4(3:4,1:2);
K(7:8,13:14) = k12(1:2,3:4);
K(13:14,7:8) = k12(1:2,3:4);
K(7:8,15:16) = k13(1:2,3:4);
K(15:16,7:8) = k13(1:2,3:4);

%Node 5
K(9:10,9:10) = k4(3:4,3:4) + k5(1:2,1:2) + k14(1:2,1:2) + k15(1:2,1:2);
K(9:10,11:12) = k5(1:2,3:4); 
K(11:12,9:10) = k5(1:2,3:4);
K(9:10,15:16) = k14(1:2,3:4);
K(15:16,9:10) = k14(1:2,3:4);
K(17:18,9:10) = k15(1:2,3:4);
K(9:10,17:18) = k15(1:2,3:4);

%Node 6
K(11:12,11:12) = k5(3:4,3:4) + k9(3:4,3:4);
K(11:12,17:18) = k9(1:2,3:4);
K(17:18,11:12) = k9(1:2,3:4);

%Node 7
K(13:14,13:14) = k6(3:4,3:4) + k7(1:2,1:2) + k11(1:2,1:2) + k12(1:2,1:2);
K(13:14,15:16) = k7(1:2,3:4);
K(15:16,13:14) = k7(1:2,3:4);

%Node 8
K(15:16,15:16) = k7(3:4,3:4) + k8(1:2,1:2) + k13(3:4,3:4) + k14(3:4,3:4);
K(15:16,17:18) = k8(1:2,3:4);
K(17:18,15:16) = k8(1:2,3:4);

%Node 9
K(17:18,17:18) = k5(3:4,3:4) + k8(1:2,1:2) + k14(3:4,3:4) + k15(3:4,3:4);

% Solve for unknown displacements (u2y, u3-9)
u = zeros(15);
u = [K(4,4:18); K(5,4:18); K(6,4:18); K(7,4:18); K(8,4:18); K(9,4:18);
     K(10,4:18); K(11,4:18); K(12,4:18); K(13,4:18); K(14,4:18); K(15,4:18); 
     K(16,4:18); K(17,4:18); K(18,4:18);];
u2 = u\[0; 0; 0; 0; 0; 0; 0; 0; P; 0; 0; 0; 0; 0; 0;]; 

U = [0; 0; 0; u2(1); u2(2); u2(3); u2(4); u2(5); u2(6); u2(7); u2(8); 
    u2(9); u2(10); u2(11); u2(12); u2(13); u2(14); u2(15);];

F = K * U;

%member forces
fe1 = k(1) * ((U(3)-U(1))*C(1) + (U(4)-U(2))*S(1));
fe2 = k(2) * ((U(5)-U(3))*C(2) + (U(6)-U(4))*S(2));
fe3 = k(3) * ((U(7)-U(5))*C(3) + (U(8)-U(6))*S(3));
fe4 = k(4) * ((U(9)-U(7))*C(4) + (U(10)-U(8))*S(4));
fe5 = k(5) * ((U(11)-U(9))*C(5) + (U(12)-U(10))*S(5));
fe6 = k(6) * ((U(13)-U(1))*C(6) + (U(14)-U(2))*S(6));
fe7 = k(7) * ((U(15)-U(13))*C(7) + (U(16)-U(14))*S(7));
fe8 = k(8) * ((U(17)-U(15))*C(8) + (U(18)-U(16))*S(8));
fe9 = k(9) * ((U(11)-U(17))*C(9) + (U(12)-U(18))*S(9));
fe10 = k(10) * ((U(5)-U(1))*C(10) + (U(6)-U(2))*S(10));
fe11 = k(11) * ((U(13)-U(5))*C(11) + (U(14)-U(6))*S(11));
fe12 = k(12) * ((U(13)-U(7))*C(12) + (U(14)-U(8))*S(12));
fe13 = k(13) * ((U(15)-U(7))*C(13) + (U(16)-U(8))*S(13));
fe14 = k(14) * ((U(15)-U(9))*C(14) + (U(16)-U(10))*S(14));
fe15 = k(15) * ((U(17)-U(9))*C(15) + (U(18)-U(10))*S(15));

Fe = [fe1;fe2;fe3;fe4;fe5;fe6;fe7;fe8;fe9;fe10;fe11;fe12;fe13;fe14;fe15;];

% Format and print results
fprintf('Nodal Displacements:\n');
j = 1;
for i = 1:2:length(U)
    fprintf('u%ix = %.3f \n', j, U(i));
    fprintf('u%iy = %.3f \n', j, U(i+1));
    j = j+1;
end

fprintf('\nNodal Forces:\n');
j = 1;
for i = 1:2:length(F)
    fprintf('f%ix = %.3f \n', j, F(i));
    fprintf('f%iy = %.3f \n', j, F(i+1));
    j = j+1;
end

fprintf('\nMember Forces:\n');
for i = 1:length(Fe)
    fprintf('f(%i) = %.3f \n', i, Fe(i));
end

%1. Pin Failure due to Shear

maxShearStrengthForDowels = 54e6;
PJoint = [59.686, 58.875, 58.875, 59.686, 59.0785];
Area =  pi*(0.0047625/2)^2;
shearAllowedForPins =  (PJoint./Area)./1e6;

shearAllowedForPinsAssert = (shearAllowedForPins < maxShearStrengthForDowels);
% Format and print test1 results
fprintf('Test 1 Results:\n');  
for i = 1:length(shearAllowedForPins)               
    if shearAllowedForPinsAssert(i) ~= true
        fprintf('Pin %i failed due to shear (Test 1 Failed). \n', i);
    else
         fprintf('Pin %i passed Pin Failure due to shear (Test 1 Passed). \n', i);
    end
end

%2. Plate rupture due to tension (MPa)

maxBasswoodStress = 65.163e6;

PlateRuptureTensionStress = (Fe/(T.*(W-0.0047625)))/(1e6);
PlateRuptureTensionsAssert = (PlateRuptureTensionStress < maxBasswoodStress);

for i = 1:length(PlateRuptureTensionStress)
    if PlateRuptureTensionsAssert(i) ~= true
        fprintf('Link %i failed the link rupture due to tension test \n', i);
    else
        fprintf('Link %i passed the link rupture due to tension test \n', i);
    end
end

%3. Bearing Failure
BearingFailureStress = (abs(PJoint)./(T.*W))/(1e6);
BearingFailureStressAssert = (BearingFailureStress < maxBasswoodStress);

for i=1:length(BearingFailureStress)
    if BearingFailureStressAssert(i) ~= true
        fprintf('Pin %i failed and caused bearing failure \n', i);
    else
        fprintf('Pin %i passed and did not cause bearing failure \n', i);
    end
end


%4. Plate/link tear out due to shear
maxBasswoodShear = 4.2e6;
LinkTearoutShear = (abs(Fe)/(2*W*T))/(1e6);
LinkTearoutShearAssert = LinkTearoutShear < maxBasswoodShear;

for i=1:length(LinkTearoutShear)
    if LinkTearoutShearAssert(i) ~= true
        fprintf('Link %i failed the link tearout test \n', i);
    else
        fprintf('Link %i passed the link tearout test \n', i);
    end
end

%5. Plate/Link Buckling due to Compression

%Area Moment of Inertia
I = (1/12)*W*(T^3);

CriticalLoad = (((pi^2)*E*I)./(L.^2));

for i = 1:length(Fe)
   if(Fe > 0)
     if abs(CriticalLoad(i)) > abs(Fe(i))
        fprintf('Link %i passed the link buckling due to compression test \n', i);
     else 
        fprintf('Link %i failed the link buckling due to compression test \n', i);       
     end
   end
end

%Mass Calculations
density_dowels = 0.747;
density_basswood = 0.385;
radius_dowel = 0.0047625;
number_pins = 9;
area_links = 0;
area_ends = pi.*(W).^2;
for i = 1:length(L) %loop through all members and get total area
    area_links = area_links + L(i) * W(i) + area_ends(i) - 2 * pi *(radius_dowel)^2;
end
area_links = area_links * 10000; %cm^2
vol_links = area_links * T * 100; %cm^3
mass_links = vol_links * density_basswood; %grams

vol_dowels = pi*(radius_dowel)^2 * number_pins * 4 * T * 1e6; %cm^3
mass_dowels = vol_dowels * density_dowels; %grams

mass_truss = mass_links + mass_dowels; %grams
PV = P/mass_truss;
