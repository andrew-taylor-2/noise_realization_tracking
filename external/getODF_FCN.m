function [coeff varargout] = getODF_FCN(dt,kt,alpha)
%Calculates the kurtosis dODF from the diffusion tensor, dt, the kurtosis
%tensor, kt, and the radial weighting power, alpha.  
%
%   COEFF is kurtosis dODF functino coefficients
%   FK is kurtosis dODF function in spherical form
%   FK2 is kurtosis dODF function in cartesian form
%
%Non-linear optmization is necessary with non-linear optimization because there's 
%only 2 parameters to vary and it is confied to the spherical surface. Brute-force 
%optimization works better with FK2 since there's not a million trigonometric 
%functions  to evaluate over and over and over again
%
%These coefficients are defined in Glenn et al. Optimization of white matter fiber 
%tractography with diffusional kurtosis imaging. NMR in Biomed. 2015.
%
%Author: Russell Glenn
%Medical University of South Carolina

D = [dt(1) dt(4) dt(5); dt(4) dt(2) dt(6); dt(5) dt(6) dt(3)]; %Diffusion Tensor

W = zeros(3,3,3,3); %Kurtosis Tensor
W(1,1,1,1) = kt(1);
W(2,2,2,2) = kt(2);
W(3,3,3,3) = kt(3);
W(1,1,1,2) = kt(4);  W(1,1,2,1) = W(1,1,1,2); W(1,2,1,1) = W(1,1,1,2); W(2,1,1,1) = W(1,1,1,2);
W(1,1,1,3) = kt(5);  W(1,1,3,1) = W(1,1,1,3); W(1,3,1,1) = W(1,1,1,3); W(3,1,1,1) = W(1,1,1,3);
W(1,2,2,2) = kt(6);  W(2,1,2,2) = W(1,2,2,2); W(2,2,1,2) = W(1,2,2,2); W(2,2,2,1) = W(1,2,2,2);
W(1,3,3,3) = kt(7);  W(3,1,3,3) = W(1,3,3,3); W(3,3,1,3) = W(1,3,3,3); W(3,3,3,1) = W(1,3,3,3);
W(2,2,2,3) = kt(8);  W(2,2,3,2) = W(2,2,2,3); W(2,3,2,2) = W(2,2,2,3); W(3,2,2,2) = W(2,2,2,3);
W(2,3,3,3) = kt(9);  W(3,2,3,3) = W(2,3,3,3); W(3,3,2,3) = W(2,3,3,3); W(3,3,3,2) = W(2,3,3,3);
W(1,1,2,2) = kt(10); W(1,2,1,2) = W(1,1,2,2); W(1,2,2,1) = W(1,1,2,2); W(2,1,1,2) = W(1,1,2,2); W(2,1,2,1) = W(1,1,2,2); W(2,2,1,1) = W(1,1,2,2);
W(1,1,3,3) = kt(11); W(1,3,1,3) = W(1,1,3,3); W(1,3,3,1) = W(1,1,3,3); W(3,1,1,3) = W(1,1,3,3); W(3,1,3,1) = W(1,1,3,3); W(3,3,1,1) = W(1,1,3,3);
W(2,2,3,3) = kt(12); W(2,3,2,3) = W(2,2,3,3); W(2,3,3,2) = W(2,2,3,3); W(3,2,2,3) = W(2,2,3,3); W(3,2,3,2) = W(2,2,3,3); W(3,3,2,2) = W(2,2,3,3);
W(1,1,2,3) = kt(13); W(1,1,3,2) = W(1,1,2,3); W(1,2,1,3) = W(1,1,2,3); W(1,2,3,1) = W(1,1,2,3); W(1,3,1,2) = W(1,1,2,3); W(1,3,2,1) = W(1,1,2,3); W(2,1,1,3) = W(1,1,2,3); W(2,1,3,1) = W(1,1,2,3); W(2,3,1,1) = W(1,1,2,3); W(3,1,1,2) = W(1,1,2,3); W(3,1,2,1) = W(1,1,2,3); W(3,2,1,1) = W(1,1,2,3);
W(1,2,2,3) = kt(14); W(1,2,3,2) = W(1,2,2,3); W(1,3,2,2) = W(1,2,2,3); W(2,1,2,3) = W(1,2,2,3); W(2,1,3,2) = W(1,2,2,3); W(2,2,1,3) = W(1,2,2,3); W(2,2,3,1) = W(1,2,2,3); W(2,3,1,2) = W(1,2,2,3); W(2,3,2,1) = W(1,2,2,3); W(3,1,2,2) = W(1,2,2,3); W(3,2,1,2) = W(1,2,2,3); W(3,2,2,1) = W(1,2,2,3);
W(1,2,3,3) = kt(15); W(1,3,2,3) = W(1,2,3,3); W(1,3,3,2) = W(1,2,3,3); W(2,1,3,3) = W(1,2,3,3); W(2,3,1,3) = W(1,2,3,3); W(2,3,3,1) = W(1,2,3,3); W(3,1,2,3) = W(1,2,3,3); W(3,1,3,2) = W(1,2,3,3); W(3,2,1,3) = W(1,2,3,3); W(3,2,3,1) = W(1,2,3,3); W(3,3,1,2) = W(1,2,3,3); W(3,3,2,1) = W(1,2,3,3);

Davg = trace(D)/3;
U = Davg*D^-1;

A1=0;B11=0;B12=0;B13=0;B22=0;B23=0;B33=0;
C1111=0;C1112=0;C1113=0;C1122=0;C1123=0;C1133=0;C1222=0;
C1223=0;C1233=0;C1333=0;C2222=0;C2223=0;C2233=0;C2333=0;C3333=0;

for i=1:3; for j = 1:3; for k = 1:3; for l = 1:3; %ITERATE THROUGH EVERYTHING ONCE

    %COEFFICIENTF FOR: 3UijWijklUkl

    A1 = A1 + 3*U(i,j)*W(i,j,k,l)*U(k,l); 

    %COEFFICIENTS FOR: -6(a+1)UijWijklVkl

    %B11*(sin(x(1))*cos(x(2)))^2   
    %B12*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))
    %B13*(sin(x(1))*cos(x(2)))*cos(x(1))
    %B22*(sin(x(1))*sin(x(2)))^2   
    %B23*(sin(x(1))*sin(x(2)))*cos(x(1))
    %B33*cos(x(1))^2   

    B0 = -6*(alpha+1)*U(i,j)*W(i,j,k,l); 

    B11 = B11 + B0*(U(k,1)*U(l,1));              
    B12 = B12 + B0*(U(k,1)*U(l,2)+U(k,2)*U(l,1));
    B13 = B13 + B0*(U(k,1)*U(l,3)+U(k,3)*U(l,1));
    B22 = B22 + B0*(U(k,2)*U(l,2));              
    B23 = B23 + B0*(U(k,2)*U(l,3)+U(k,3)*U(l,2));
    B33 = B33 + B0*(U(k,3)*U(l,3));  

    %COEFFICIENTS FOR: (alpha+1)(alpha+3)W(i,j,k,l)VijVkl

    %C1111*(sin(x(1))*cos(x(2)))^4          
    %C1112*(sin(x(1))*cos(x(2)))^3*(sin(x(1))*sin(x(2)))     
    %C1113*(sin(x(1))*cos(x(2)))^3*cos(x(1))     
    %C1122*(sin(x(1))*cos(x(2)))^2*(sin(x(1))*sin(x(2)))^2   
    %C1123*(sin(x(1))*cos(x(2)))^2*(sin(x(1))*sin(x(2)))*cos(x(1))
    %C1133*(sin(x(1))*cos(x(2)))^2*cos(x(1))^2   
    %C1222*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))^3     
    %C1223*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))^2*cos(x(1))
    %C1233*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))*cos(x(1))^2
    %C1333*(sin(x(1))*cos(x(2)))*cos(x(1))^3     
    %C2222*(sin(x(1))*sin(x(2)))^4          
    %C2223*(sin(x(1))*sin(x(2)))^3*cos(x(1))     
    %C2233*(sin(x(1))*sin(x(2)))^2*cos(x(1))^2   
    %C2333*(sin(x(1))*sin(x(2)))*cos(x(1))^3     
    %C3333*cos(x(1))^4 

    C0 = (alpha+1)*(alpha+3)*W(i,j,k,l);

    C1111 = C1111 + C0*(U(i,1)*U(j,1)*U(k,1)*U(l,1));                                                                                                                                                                                                                                                                                                                    
    C1112 = C1112 + C0*(U(i,1)*U(j,1)*U(k,1)*U(l,2)+U(i,1)*U(j,1)*U(k,2)*U(l,1)+U(i,1)*U(j,2)*U(k,1)*U(l,1)+U(i,2)*U(j,1)*U(k,1)*U(l,1));                                                                                                                                                                                                                                
    C1113 = C1113 + C0*(U(i,1)*U(j,1)*U(k,1)*U(l,3)+U(i,1)*U(j,1)*U(k,3)*U(l,1)+U(i,1)*U(j,3)*U(k,1)*U(l,1)+U(i,3)*U(j,1)*U(k,1)*U(l,1));                                                                                                                                                                                                                                
    C1122 = C1122 + C0*(U(i,1)*U(j,1)*U(k,2)*U(l,2)+U(i,1)*U(j,2)*U(k,1)*U(l,2)+U(i,1)*U(j,2)*U(k,2)*U(l,1)+U(i,2)*U(j,1)*U(k,1)*U(l,2)+...
        U(i,2)*U(j,1)*U(k,2)*U(l,1)+U(i,2)*U(j,2)*U(k,1)*U(l,1));                                                                                                                                                                        
    C1123 = C1123 + C0*(U(i,1)*U(j,1)*U(k,2)*U(l,3)+U(i,1)*U(j,1)*U(k,3)*U(l,2)+U(i,1)*U(j,2)*U(k,1)*U(l,3)+U(i,1)*U(j,2)*U(k,3)*U(l,1)+...
        U(i,1)*U(j,3)*U(k,1)*U(l,2)+U(i,1)*U(j,3)*U(k,2)*U(l,1)+U(i,2)*U(j,1)*U(k,1)*U(l,3)+U(i,2)*U(j,1)*U(k,3)*U(l,1)+U(i,2)*U(j,3)*U(k,1)*U(l,1)+...
        U(i,3)*U(j,1)*U(k,1)*U(l,2)+U(i,3)*U(j,1)*U(k,2)*U(l,1)+U(i,3)*U(j,2)*U(k,1)*U(l,1));
    C1133 = C1133 + C0*(U(i,1)*U(j,1)*U(k,3)*U(l,3)+U(i,1)*U(j,3)*U(k,1)*U(l,3)+U(i,1)*U(j,3)*U(k,3)*U(l,1)+U(i,3)*U(j,1)*U(k,1)*U(l,3)+...
        U(i,3)*U(j,1)*U(k,3)*U(l,1)+U(i,3)*U(j,3)*U(k,1)*U(l,1));                                                                                                                                                                        
    C1222 = C1222 + C0*(U(i,1)*U(j,2)*U(k,2)*U(l,2)+U(i,2)*U(j,1)*U(k,2)*U(l,2)+U(i,2)*U(j,2)*U(k,1)*U(l,2)+U(i,2)*U(j,2)*U(k,2)*U(l,1));                                                                                                                                                                                                                                
    C1223 = C1223 + C0*(U(i,1)*U(j,2)*U(k,2)*U(l,3)+U(i,1)*U(j,2)*U(k,3)*U(l,2)+U(i,1)*U(j,3)*U(k,2)*U(l,2)+U(i,2)*U(j,1)*U(k,2)*U(l,3)+...
        U(i,2)*U(j,1)*U(k,3)*U(l,2)+U(i,2)*U(j,2)*U(k,1)*U(l,3)+U(i,2)*U(j,2)*U(k,3)*U(l,1)+U(i,2)*U(j,3)*U(k,1)*U(l,2)+U(i,2)*U(j,3)*U(k,2)*U(l,1)+...
        U(i,3)*U(j,1)*U(k,2)*U(l,2)+U(i,3)*U(j,2)*U(k,1)*U(l,2)+U(i,3)*U(j,2)*U(k,2)*U(l,1));
    C1233 = C1233 + C0*(U(i,1)*U(j,2)*U(k,3)*U(l,3)+U(i,1)*U(j,3)*U(k,2)*U(l,3)+U(i,1)*U(j,3)*U(k,3)*U(l,2)+U(i,2)*U(j,1)*U(k,3)*U(l,3)+...
        U(i,2)*U(j,3)*U(k,1)*U(l,3)+U(i,2)*U(j,3)*U(k,3)*U(l,1)+U(i,3)*U(j,1)*U(k,2)*U(l,3)+U(i,3)*U(j,1)*U(k,3)*U(l,2)+U(i,3)*U(j,2)*U(k,1)*U(l,3)+...
        U(i,3)*U(j,2)*U(k,3)*U(l,1)+U(i,3)*U(j,3)*U(k,1)*U(l,2)+U(i,3)*U(j,3)*U(k,2)*U(l,1));
    C1333 = C1333 + C0*(U(i,1)*U(j,3)*U(k,3)*U(l,3)+U(i,3)*U(j,1)*U(k,3)*U(l,3)+U(i,3)*U(j,3)*U(k,1)*U(l,3)+U(i,3)*U(j,3)*U(k,3)*U(l,1));                                                                                                                                                                                                                                
    C2222 = C2222 + C0*(U(i,2)*U(j,2)*U(k,2)*U(l,2));                                                                                                                                                                                                                                                                                                                    
    C2223 = C2223 + C0*(U(i,2)*U(j,2)*U(k,2)*U(l,3)+U(i,2)*U(j,2)*U(k,3)*U(l,2)+U(i,2)*U(j,3)*U(k,2)*U(l,2)+U(i,3)*U(j,2)*U(k,2)*U(l,2));                                                                                                                                                                                                                                
    C2233 = C2233 + C0*(U(i,2)*U(j,2)*U(k,3)*U(l,3)+U(i,2)*U(j,3)*U(k,2)*U(l,3)+U(i,2)*U(j,3)*U(k,3)*U(l,2)+U(i,3)*U(j,2)*U(k,2)*U(l,3)+...
        U(i,3)*U(j,2)*U(k,3)*U(l,2)+U(i,3)*U(j,3)*U(k,2)*U(l,2));                                                                                                                                                                        
    C2333 = C2333 + C0*(U(i,2)*U(j,3)*U(k,3)*U(l,3)+U(i,3)*U(j,2)*U(k,3)*U(l,3)+U(i,3)*U(j,3)*U(k,2)*U(l,3)+U(i,3)*U(j,3)*U(k,3)*U(l,2));                                                                                                                                                                                                                                
    C3333 = C3333 + C0*(U(i,3)*U(j,3)*U(k,3)*U(l,3));         

end; end; end; end;

coeff = [A1;B11;B12;B13;B22;B23;B33;...
C1111;C1112;C1113;C1122;C1123;C1133;C1222;...
C1223;C1233;C1333;C2222;C2223;C2233;C2333;C3333;...
U(1,1);U(2,2);U(3,3);U(1,2);U(1,3);U(2,3);alpha];


if nargout>1; 

%GET KURTOSIS ODF FUNCTION IN SPHERICAL FORM

varargout{1} = @(x,A)-(1./((sin(x(:,1)).*cos(x(:,2))).^2.*A(23)+(sin(x(:,1)).*sin(x(:,2))).^2.*A(24)+cos(x(:,1)).^2.*A(25)+2.*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).*A(26)+...
2.*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1)).*A(27)+2.*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).*A(28))).^((A(29)+1)./2).*(1+(A(1)+...
(A(2).*(sin(x(:,1)).*cos(x(:,2))).^2+A(3).*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2)))+A(4).*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1))+A(5).*(sin(x(:,1)).*sin(x(:,2))).^2+...
A(6).*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1))+A(7).*cos(x(:,1)).^2)./((sin(x(:,1)).*cos(x(:,2))).^2.*A(23)+(sin(x(:,1)).*sin(x(:,2))).^2.*A(24)+cos(x(:,1)).^2.*A(25)+...
2.*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).*A(26)+2.*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1)).*A(27)+2.*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).*A(28))+...
(A(8).*(sin(x(:,1)).*cos(x(:,2))).^4+A(9).*(sin(x(:,1)).*cos(x(:,2))).^3.*(sin(x(:,1)).*sin(x(:,2)))+A(10).*(sin(x(:,1)).*cos(x(:,2))).^3.*cos(x(:,1))+A(11).*(sin(x(:,1)).*cos(x(:,2))).^2.*(sin(x(:,1)).*sin(x(:,2))).^2+...
A(12).*(sin(x(:,1)).*cos(x(:,2))).^2.*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1))+A(13).*(sin(x(:,1)).*cos(x(:,2))).^2.*cos(x(:,1)).^2+A(14).*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).^3+...
A(15).*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).^2.*cos(x(:,1))+A(16).*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).^2+A(17).*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1)).^3+A(18).*(sin(x(:,1)).*sin(x(:,2))).^4+...
A(19).*(sin(x(:,1)).*sin(x(:,2))).^3.*cos(x(:,1))+A(20).*(sin(x(:,1)).*sin(x(:,2))).^2.*cos(x(:,1)).^2+A(21).*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).^3+A(22).*cos(x(:,1)).^4)./((sin(x(:,1)).*cos(x(:,2))).^2.*A(23)+(sin(x(:,1)).*sin(x(:,2))).^2.*A(24)+...
cos(x(:,1)).^2.*A(25)+2.*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).*A(26)+2.*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1)).*A(27)+2.*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).*A(28)).^2)./24);


if nargout > 2

%GET KURTOSIS ODF FUNCTION IN CARTESIAN FORM

varargout{2} = @(x,A)-(1./((x(:,1)).^2.*A(23)+(x(:,2)).^2.*A(24)+x(:,3).^2.*A(25)+2.*(x(:,1)).*(x(:,2)).*A(26)+...
2.*(x(:,1)).*x(:,3).*A(27)+2.*(x(:,2)).*x(:,3).*A(28))).^((A(29)+1)./2).*(1+(A(1)+...
(A(2).*(x(:,1)).^2+A(3).*(x(:,1)).*(x(:,2))+A(4).*(x(:,1)).*x(:,3)+A(5).*(x(:,2)).^2+...
A(6).*(x(:,2)).*x(:,3)+A(7).*x(:,3).^2)./((x(:,1)).^2.*A(23)+(x(:,2)).^2.*A(24)+x(:,3).^2.*A(25)+...
2.*(x(:,1)).*(x(:,2)).*A(26)+2.*(x(:,1)).*x(:,3).*A(27)+2.*(x(:,2)).*x(:,3).*A(28))+...
(A(8).*(x(:,1)).^4+A(9).*(x(:,1)).^3.*(x(:,2))+A(10).*(x(:,1)).^3.*x(:,3)+A(11).*(x(:,1)).^2.*(x(:,2)).^2+...
A(12).*(x(:,1)).^2.*(x(:,2)).*x(:,3)+A(13).*(x(:,1)).^2.*x(:,3).^2+A(14).*(x(:,1)).*(x(:,2)).^3+...
A(15).*(x(:,1)).*(x(:,2)).^2.*x(:,3)+A(16).*(x(:,1)).*(x(:,2)).*x(:,3).^2+A(17).*(x(:,1)).*x(:,3).^3+A(18).*(x(:,2)).^4+...
A(19).*(x(:,2)).^3.*x(:,3)+A(20).*(x(:,2)).^2.*x(:,3).^2+A(21).*(x(:,2)).*x(:,3).^3+A(22).*x(:,3).^4)./((x(:,1)).^2.*A(23)+(x(:,2)).^2.*A(24)+...
x(:,3).^2.*A(25)+2.*(x(:,1)).*(x(:,2)).*A(26)+2.*(x(:,1)).*x(:,3).*A(27)+2.*(x(:,2)).*x(:,3).*A(28)).^2)./24);

if nargout > 3
    
%GET GAUSSIAN ODF FUNCTION IN SPHERICAL FORM    
    
varargout{3} = @(x,A)-(1./((sin(x(:,1)).*cos(x(:,2))).^2.*A(23)+(sin(x(:,1)).*sin(x(:,2))).^2.*A(24)+cos(x(:,1)).^2.*A(25)+2.*(sin(x(:,1)).*cos(x(:,2))).*(sin(x(:,1)).*sin(x(:,2))).*A(26)+...
2.*(sin(x(:,1)).*cos(x(:,2))).*cos(x(:,1)).*A(27)+2.*(sin(x(:,1)).*sin(x(:,2))).*cos(x(:,1)).*A(28))).^((A(29)+1)./2); 

if nargout > 4
    
%GET GAUSSIAN ODF FUNCTION IN CARTESIAN FORM 
    
varargout{4} = @(x,A)-(1./((x(:,1)).^2.*A(23)+(x(:,2)).^2.*A(24)+x(:,3).^2.*A(25)+2.*(x(:,1)).*(x(:,2)).*A(26)+...
2.*(x(:,1)).*x(:,3).*A(27)+2.*(x(:,2)).*x(:,3).*A(28))).^((A(29)+1)./2);

end
end
end
end