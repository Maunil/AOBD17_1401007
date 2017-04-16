%LDAs Updation 

clear;
clc;

load fisheriris

meas = meas(1:100,3:4)';
species = species(1:100)';

xc1 = mean(meas(:,1:49)');
xc2 = mean(meas(:,51:99)');

x_bar = (xc1 + xc2)./2;

sw1 = zeros(2,2);
sw2 = zeros(2,2);
sw = zeros(2,2);

no_points = 49;

for i=1:1:no_points
     sw1 = sw1 +(meas(:,i)-xc1')*(meas(:,i)-xc1')';   
end

for j=51:1:no_points+50
    sw2 = sw2 +(meas(:,j)-xc2')*(meas(:,j)-xc2')';
end

sw = sw1+sw2;

sb1 = no_points*(xc1-x_bar)*(xc1-x_bar)'; 
sb2 =  no_points*(xc2-x_bar)*(xc2-x_bar)';
sb = sb1 + sb2;

lda= inv(sw)*sb;

[v d] = eig(lda);

new_data = v'*meas;


%%To check inner and outer distance varies 

% new_mean1 = mean(new_data(:,1:50)');
% new_mean2 = mean(new_data(:,51:100)');
% 
% %new_mean = mean(new_data(:,1:100)');
% new_mean = (new_mean1 + new_mean2)./2;
% 
% %Inner Distance 
% disp(sum((new_data(:,1)-new_mean1').^2));
% disp(sum((meas(:,1)-xc1').^2));
% 
% 
% %Outer Distance 
% disp(sum((new_mean2-new_mean1).^2));
% disp(sum((xc2-xc1).^2));
% 


%Online version 
%New point added
new_point = meas(:,50);

%Labeling check 
dist1=(sum(new_point'-xc1))^2;
dist2=(sum(new_point'-xc2))^2;

if(dist1<dist2)
    class = 1;
else
    class = 2;
end


%Updating sw and sb
if(class==1)
    new_sw = sw2+sw1+(no_points/(no_points+1)).*((new_point- xc1')*(new_point-xc1')');
    new_x_bar = (2*no_points.*x_bar+ new_point')/(2*no_points+1);
    new_xc1 = (no_points.*xc1+ new_point')/(no_points+1);
    updated_sb1 = (no_points+1).*(new_xc1-new_x_bar)*(new_xc1-new_x_bar)';
    new_sb = updated_sb1+sb;    
end

%Computed New Ldas
lda_updated= inv(new_sw)*new_sb;

[v_updated d_updated] = eig(lda_updated);

































