function [Gam_opt,Kopt, Xopt,Yopt] = SLPV(Gsym,plant_paras,design_paras)
%% LPV controller design based on the synthesis conditions presented in [Apkarian 1998, Advanced gain-scheduling techniques for uncertain systems] 
% Note
% (1) can be used for both one and two scheduling parameters. For one parameter case, need to set:
%       - Theta2Grid = 0,  d_Theta ={xx,0},
%       - Fcn_theta is only function of theta1, e.g., Fcn_theta = @(x) [1 x(1)]  
% (4) include both basic and projected conditions presented in [Apkarian 1998, Advanced gain-scheduling techniques for uncertain systems] 

% Input parameters
% Gsym: generalized plant with symbol variables
% ThetaGrid, gridded points of theta1
% Theta2Grid, gridded points of theta2
% d_Theta, cell array, bounds for derivatives of theta, d_Thetah{1} for theta1, d_Thetah{2} for theta2

% XY_PD: choose whether X or Y is constant. 
    %XY_PD=1: Y is constant;
    %XY_PD=2: X is constant;
    %XY_PD=0: both X and Y are constant;
    %XY_PD=3: both X and Y are parameter dependent
% Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
% the matrix variables, for instance in X = X0+ f1(theta)X1+f2(theta)X2,
% Fcn_theta(theta) = [1 f1(theta) f2(theta)]
% d_Fcn_theta: a function handle for the derivative of the scalar functions
% FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as a function of theta1, one matrix as a function of theta2,
%            while [1 2] means one constant matrix, two matrices as a function of theta1, no theta2. 
% IsProject, boolean, 0 for using the basic conditions (Theorem 2.1 in [Apkarian 1998], 1 for using the projected condition (Theorem 2.2 in [Apkarian 1998]
% SolverSel, integer, 1 for using LMI Lab, 2 for using CVX (not complete yet)   
% StrucUncert, boolean, 0 for normal LPV synthesis, 1 for LPV synthesis considering the structured uncertainty (Section V in [Apkarian 98])

% Output parameters:
% Gam_opt: optimal gamma value
% Kopt: structure array, contains the matrices Ahat, Bhat, Chat and Dhat  
% Xopt: optimal X
% Yopt: optimal Y
% Fcn_theta: function handle, Fcn_theta(theta) will give all the scalar functions for the matrix variables. To realize X = X0+ f1(theta)X1+f2(theta)X2, set Fcn_theta(theta) = [1 f1(theta) f2(theta)]
% d_Fcn_theta: function handle for the derivative of the scalar functions
% FthetaNum: 1 by 3 vector in the form of [1 n1 n2], where n1 (n2) denotes the number of scalar functions depending on theta1 (theta2)

Theta1Grid= plant_paras.Theta1Grid;
Theta2Grid = plant_paras.Theta2Grid;
d_Theta = plant_paras.d_Theta;
StrucUncert = plant_paras.StrucUncert;

XY_PD = design_paras.XY_PD;
Fcn_theta = design_paras.Fcn_theta;
d_Fcn_theta = design_paras.d_Fcn_theta;
FthetaNum = design_paras.FthetaNum;
IsProject = design_paras.IsProject;
SolverSel = design_paras.SolverSel;

% eigShift = 1.1; % to improve condition number 
Sv = 1;
Svinf = 1/Sv;
% if isempty(Theta2) || norm(Theta2,'inf') == 0
%     GSParaNum = 1;
% else
%     GSParaNum = 2;
% end

n =  size(Gsym.A(:,:,1),1);
nw = size(Gsym.B1(:,:,1),2);
nu = size(Gsym.B2(:,:,1),2);
nz = size(Gsym.C1(:,:,1),1);
ny = size(Gsym.C2(:,:,1),1);
%% Scaling matrix defined for considering structured uncertainty
% need to modify this for specific examples
if StrucUncert == 1
%     S = sqrtm(diag([0.5e-3,0.5e-3,1,1]));
    S = diag([0.5e-3,0.5e-3,1,1]);
    Sinv = inv(S);
else
    S = eye(nw);
    Sinv = eye(nz);
end
%% Defining LMI variable
% gamma,1  
% X_0,2
% Ahat_0,3
% Bhat_0,4
% Chat_0,5
% Dhat_0,6
% X_1,7
% Ahat_1,8
% Bhat_1,9
% Chat_1,10
% Dhat_1,11
% X_2,12
% Ahat_2,13
% Bhat_2,14
% Chat_2,15
% Dhat_2,16
% ...
% Y_0, 17
if XY_PD == 0
    d_Theta = {0,0};
    disp('Parameter variation rate is set to 0 due to use of constant Lyapunov function, or both parameter-dependent X and Y');
end
if SolverSel == 1 %LMI lab  
        setlmis([]);
        Gam = lmivar(1,[1 1]);
        %% Note the order of matrices: X = X0+ theta1*X1+ theta1^2*X2+ theta2*X3+ theta2^2*X4
        gam = lmivar(1,[1 1]);
        for Id_Ftheta=1:sum(FthetaNum)% 
            if XY_PD == 1           %Y0 is constatnt, so define X
               X(Id_Ftheta)=lmivar(1,[n 1]); 
            elseif XY_PD == 2       %X0 is constatnt
               Y(Id_Ftheta)=lmivar(1,[n 1]);
            elseif XY_PD == 3       %both X and Y are parameter-dependent
               X(Id_Ftheta)=lmivar(1,[n 1]); 
               Y(Id_Ftheta)=lmivar(1,[n 1]); 
            end
            K.Ah(Id_Ftheta)=lmivar(2,[n n]);
            K.Bh(Id_Ftheta)=lmivar(2,[n ny]);
            K.Ch(Id_Ftheta)=lmivar(2,[nu n]);
            K.Dh(Id_Ftheta)=lmivar(2,[nu ny]); 
        end            
         if XY_PD == 1 %Y0 is constant 
            Y = lmivar(1,[n 1]);  %Y = Y0,constant
        elseif XY_PD ==2 
            X = lmivar(1,[n 1]);  % X = X0,constant
        elseif XY_PD == 0 % X= X0; Y = Y0;
            X = lmivar(1,[n 1]); 
            Y = lmivar(1,[n 1]);
         end
         
        lminum = 0;       
        %% LMI conditions              
        for Id_theta1 = 1:length(Theta1Grid) %Repeat loop based on Number of gridded parameters
            for Id_theta2= 1:length(Theta2Grid) 
                theta = [Theta1Grid(Id_theta1); Theta2Grid(Id_theta2)];
                Ga = AugPltEval(Gsym,theta);
                % Get Nx and Ny as the basis for null space of ..., for
                % projected technique
                Nx = null([Ga.C2 Ga.D21]);
                Ny = null([Ga.B2' Ga.D12']);                    
                %  clear theta
                %% Function of theta and its derivative                          
%                         Ftheta = [1 theta(1,1) theta(2,1)];
%                         d_Ftheta = [0 1 1];
                % For Wu96 Example            
                % Ftheta = [1 cos(theta) sin(theta)];
                % d_Ftheta = [0 -sin(theta) cos(theta)];
                Ftheta = Fcn_theta(theta);
                d_Ftheta = d_Fcn_theta(theta);

                for d_theta1 = d_Theta{1}                           
                    for d_theta2 = d_Theta{2}
                        if IsProject == 0
                            lminum = lminum+1;                            
                            for Id_Ftheta=1:sum(FthetaNum) %Repeat based on number of Ftheta
                                % Equation (7)
                                %(1,1) Dromin*DX+XA+Bhat*C2+(*)
                                %(2,1) Ahatk'+A+B2*Dhatk*C2
                                %(2,2) -DY+A*Y+B2*Chat+(*)
                                %(3,1) (X*B1+Bhatk*D21)'
                                %(3,2) (B1+B2*Dk*D21)'
                                %(3,3) -gamma
                                %(4,1) C1+D12*Dk*C2
                                %(4,2) C1*Y+D12*Chat
                                %(4,3) D11+D12*Dk*D21
                                %(4,4) -gamma
                                if XY_PD == 1 % Y = Y0;
                                    if Id_Ftheta >1 && Id_Ftheta <= 1+FthetaNum(2) %% add the Xdot
                                        lmiterm([lminum 1 1 X(Id_Ftheta)],d_Ftheta(Id_Ftheta)*d_theta1,1);
                                    elseif Id_Ftheta > 1+FthetaNum(2) && Id_Ftheta <= sum(FthetaNum)
                                        lmiterm([lminum 1 1 X(Id_Ftheta)],d_Ftheta(Id_Ftheta)*d_theta2,1);                                       
                                    end
                                    %add the term X.
                                    lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A,'s'); %XA+A'X
                                    lmiterm([lminum 3 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B1',1);
                                elseif XY_PD == 2 % X is constant
                                    if  Id_Ftheta>1 && Id_Ftheta <= 1+FthetaNum(2) %% add the Xdot
                                        lmiterm([lminum 2 2 Y(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta1,1);
                                    elseif Id_Ftheta > 1+FthetaNum(2) && Id_Ftheta <= sum(FthetaNum)
                                        lmiterm([lminum 2 2 Y(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta2,1);                                       
                                    end
                                    lmiterm([lminum 2 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.A,1,'s');
                                    lmiterm([lminum 4 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.C1,1);
                                elseif XY_PD == 3                                        
                                    lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A,'s'); %XA+A'X
                                    lmiterm([lminum 3 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B1',1);
                                    lmiterm([lminum 2 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.A,1,'s');
                                    lmiterm([lminum 4 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.C1,1);                                        
                                end
                                lmiterm([lminum 1 1  K.Bh(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.C2,'s');
                                lmiterm([lminum 2 1 -K.Ah(Id_Ftheta)],Ftheta(Id_Ftheta),1);
                                lmiterm([lminum 2 1  K.Dh(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B2,Ga.C2);
                                lmiterm([lminum 2 2  K.Ch(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B2,1,'s');
                                lmiterm([lminum 3 1 -K.Bh(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.D21',1);
                                lmiterm([lminum 3 2 -K.Dh(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.D21',Ga.B2');
                                lmiterm([lminum 4 1  K.Dh(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.D12,Ga.C2);
                                lmiterm([lminum 4 2  K.Ch(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.D12,1);
                                lmiterm([lminum 4 3  K.Dh(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.D12,Ga.D21);
                            end

                            %add the term containing X0 or Y0
                            if XY_PD == 1 
                                lmiterm([lminum 2 2 Y],Ga.A,1,'s');
                                lmiterm([lminum 4 2 Y],Ga.C1,1);
                            elseif XY_PD == 2
                                lmiterm([lminum 1 1 X],1,Ga.A,'s');
                                lmiterm([lminum 3 1 X],Ga.B1',1);
                            elseif XY_PD == 0
                                lmiterm([lminum 1 1 X],1,Ga.A,'s'); %XA+A'X
                                lmiterm([lminum 3 1 X],Ga.B1',1);
                                lmiterm([lminum 2 2 Y],Ga.A,1,'s');
                                lmiterm([lminum 4 2 Y],Ga.C1,1);
                            end                                      
                            lmiterm([lminum 2 1 0],Ga.A);
                            lmiterm([lminum 3 2 0],Ga.B1');
                            lmiterm([lminum 3 3 gam],-Sv,1);
                            lmiterm([lminum 4 1 0],Ga.C1);
                            lmiterm([lminum 4 3 0],Ga.D11);
                            lmiterm([lminum 4 4 gam],-Svinf,1);     
                        elseif IsProject == 1
                            lminum = lminum + 1;
                            lmiterm([lminum 0 0 0],[Nx zeros(n+nw,nz);zeros(nz,n+nw-ny) eye(nz)]);% outer factor
                            for Id_Ftheta=1:sum(FthetaNum)%Repeat based on number of Ftheta                                                        
                                if XY_PD == 1 % Y = Y0;     
                                    if Id_Ftheta >1 && Id_Ftheta <= 1+FthetaNum(2) %% add the Xdot
                                        lmiterm([lminum 1 1 X(Id_Ftheta)],d_Ftheta(Id_Ftheta)*d_theta1,1);
                                    elseif Id_Ftheta > 1+FthetaNum(2) && Id_Ftheta <= sum(FthetaNum)
                                        lmiterm([lminum 1 1 X(Id_Ftheta)],d_Ftheta(Id_Ftheta)*d_theta2,1);                                       
                                    end                                        
                                    %add the term X.
                                    lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A,'s'); %XA+A'X
                                    lmiterm([lminum 2 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B1',1);                                                         
                                elseif XY_PD == 3
                                    lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A,'s'); %XA+A'X
                                    lmiterm([lminum 2 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.B1',1); 
                                end
                            end
                            if XY_PD == 2 || XY_PD == 0
                               lmiterm([lminum 1 1 X],1,Ga.A,'s');
                               lmiterm([lminum 2 1 X],Ga.B1',1);   
                            end                            
                            lmiterm([lminum 3 1 0],Ga.C1);
                            lmiterm([lminum 3 2 0],Ga.D11);
                            lmiterm([lminum 2 2 gam],-Sv*S,1); % S and Sinv are for structured uncertainty
                            lmiterm([lminum 3 3 gam],-Svinf*Sinv,1);                

                        %% Equation (13) in Apkarian 98
                            lminum = lminum + 1;
                            lmiterm([lminum 0 0 0],[Ny zeros(n+nz,nw);zeros(nw,n+nz-nu) eye(nw)]); % outer factor
                            for Id_Ftheta=1:sum(FthetaNum)%Repeat based on number of Ftheta                                                       
                                if XY_PD == 2 % Y = Y(theta), X= X0;
                                    if  Id_Ftheta>1 && Id_Ftheta <= 1+FthetaNum(2) %% add the Xdot
                                        lmiterm([lminum 1 1 Y(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta1,1);
                                    elseif Id_Ftheta > 1+FthetaNum(2) && Id_Ftheta <= sum(FthetaNum)
                                        lmiterm([lminum 1 1 Y(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta2,1);                                       
                                    end
                                    %add the term X.
                                    lmiterm([lminum 1 1 Y(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A','s'); %YA'+AY
                                    lmiterm([lminum 2 1 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.C1,1);                               
                                elseif XY_PD == 3
                                    lmiterm([lminum 1 1 Y(Id_Ftheta)],Ftheta(Id_Ftheta),Ga.A','s'); %XA+A'X
                                    lmiterm([lminum 2 1 Y(Id_Ftheta)],Ftheta(Id_Ftheta)*Ga.C1,1); 
                                end
                            end
                            if XY_PD == 1 || XY_PD == 0
                               lmiterm([lminum 1 1 Y],1,Ga.A','s');
                               lmiterm([lminum 2 1 Y],Ga.C1,1);  
                            end  
                            lmiterm([lminum 3 1 0],Ga.B1');
                            lmiterm([lminum 3 2 0],Ga.D11');
                            lmiterm([lminum 2 2 gam],-Sv*Sinv,1); % S and Sinv are for structured uncertainty
                            lmiterm([lminum 3 3 gam],-Svinf*S,1);     
                        end                            
                    end %Id_d_theta2
                end %Id_d_theta1

                %% Equation (8) [-X -I;-I -Y]<0
                lminum = lminum + 1;
                %     lmiterm([lminum 2 1 0],-Eigshift);
                lmiterm([lminum 2 1 0],-1);
                if XY_PD == 1
                    lmiterm([lminum 2 2 Y],1,-1);
                    for Id_Ftheta=1:sum(FthetaNum)
                        lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),-1); %%%%%No x0
                    end
                elseif XY_PD ==2
                    lmiterm([lminum 1 1 X],1,-1);
                    for Id_Ftheta=1:sum(FthetaNum)
                        lmiterm([lminum 2 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta),-1);
                    end
                elseif XY_PD == 3
                    for Id_Ftheta=1:sum(FthetaNum)
                        lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),-1);
                        lmiterm([lminum 2 2 Y(Id_Ftheta)],Ftheta(Id_Ftheta),-1);                        
                    end 
                elseif XY_PD == 0
                    lmiterm([lminum 1 1 X],1,-1);
                    lmiterm([lminum 2 2 Y],1,-1);
                end
                % for avoiding ill-conditioned matrices, we impose [X I; I Y] >= epsilon*I
                % % should comment this if just interested in gamma value computation 
                epsilon = 1e-6; % need to tune for each example, large number may lead to infeasiblity                    
                lmiterm([lminum 1 1 0],epsilon);  
                lmiterm([lminum 2 2 0],epsilon); 
            end %Id_theta2
        end  %Id_theta1    
        % gam < Gam
        lminum = lminum + 1;
        lmiterm([lminum 1 1 gam],1,1);
        lmiterm([lminum 1 1 Gam],-1,1);  

        % for avoiding big matrices that may cause numerical problems, impose [X(1) Y(1] <= epsilon2*I
        % should comment this if just interested in gamma value computation 
%         epsilon2 = 1e5; % need to tune for each example, small number may lead to infeasiblity
%         lminum = lminum + 1;
%         lmiterm([lminum 1 1 X(1)],1,1);
%         lmiterm([lminum 2 2 Y(1)],1,1);
%         lmiterm([lminum 1 1 0],-epsilon2);
%         lmiterm([lminum 2 2 0],-epsilon2);  

        %% Get and slove LMIs
        lminum
        lmisys = getlmis;        
        

        nvar = decnbr(lmisys); %number of decision variable.
        c = zeros(nvar,1);
        c(1)=1; 

        options(1)= 1e-4;
        options(2)= 500; %Number of Iteration
        options(3) = 2e5; % for avoiding large numbers in the decision variables
        options(4) =  5; % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
        options(5)= 0; % 1 not show the process
                        %% initial value
                        % if stpcont~=1
                        %     if XY_PD==1 %Y0 is constatnt
                        %         xinit = mat2dec(lmisys,gamma0,0(:,:,1),Khat0.A(:,:,1),Khat0.B(:,:,1),Khat0.C(:,:,1),Khat0.D(:,:,1),X00(:,:,2),Khat0.A(:,:,2),Khat0.B(:,:,2),Khat0.C(:,:,2),Khat0.D(:,:,2),X00(:,:,3),Khat0.A(:,:,3),Khat0.B(:,:,3),Khat0.C(:,:,3),Khat0.D(:,:,3),Y00);
                        %     else
                        %         xinit = mat2dec(lmisys,gamma0,Y00(:,:,1),Khat0.A(:,:,1),Khat0.B(:,:,1),Khat0.C(:,:,1),Khat0.D(:,:,1),Y00(:,:,2),Khat0.A(:,:,2),Khat0.B(:,:,2),Khat0.C(:,:,2),Khat0.D(:,:,2),Y00(:,:,3),Khat0.A(:,:,3),Khat0.B(:,:,3),Khat0.C(:,:,3),Khat0.D(:,:,3),X00);
                        %     end
                        %     [gammaopt,xopt] = mincx(lmisys,c,options,xinit);
                        % else
                        % end
        [Gam_opt,xopt] = mincx(lmisys,c,options);
        % [gammaopt,xopt] = mincx(lmisys,c);
        % [gammaopt,xopt] = mincx(lmisys,c,options,[Khatopt0,Yopt0,Xopt0]);

        %% retrive the decision matrices          
        for Id_Ftheta = 1:sum(FthetaNum)
            Kopt.Ah(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,K.Ah(Id_Ftheta));
            Kopt.Bh(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,K.Bh(Id_Ftheta));
            Kopt.Ch(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,K.Ch(Id_Ftheta));
            Kopt.Dh(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,K.Dh(Id_Ftheta));
            switch XY_PD
                case 0
                    if Id_Ftheta == 1
                        Xopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,X);
                        Yopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,Y);
                    else
                        Xopt(:,:,Id_Ftheta) = zeros(n,n);
                        Yopt(:,:,Id_Ftheta) = zeros(n,n);                            
                    end                            
                case 1
                    Xopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,X(Id_Ftheta));
                    if Id_Ftheta == 1
                        Yopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,Y);
                    else
                        Yopt(:,:,Id_Ftheta) = zeros(n,n);
                    end
                case 2   
                    Yopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,Y(Id_Ftheta));
                    if Id_Ftheta == 1
                        Xopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,X);
                    else
                        Xopt(:,:,Id_Ftheta) = zeros(n,n);
                    end 
                case 3
                    Xopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,X(Id_Ftheta));
                    Yopt(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,Y(Id_Ftheta));
            end
        end
%             if XY_PD == 1
%                 Yopt(:,:) = dec2mat(lmisys,xopt,Y);
%             elseif XY_PD ==2
%                 Xopt(:,:) = dec2mat(lmisys,xopt,X);
%             elseif XY_PD == 0
%                 Xopt(:,:) = dec2mat(lmisys,xopt,X);
%                 Yopt(:,:) = dec2mat(lmisys,xopt,Y);
%             end                  
elseif SOLVER == 2 % using YALMIP to formulate and solve the LMI problem. To be completed. 
  
end

