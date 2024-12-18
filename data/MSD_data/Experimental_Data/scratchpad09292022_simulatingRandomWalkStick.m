%%
%https://jb.asm.org/content/192/18/4535
%  Dgfp=3.5*10^5 %um^2/s
%  agfp =4*10^-3 ;%4nm x 2 nm
%  a =50*10^-3; %40nm
%  D = Dgfp*agfp/a;%D1a1=D2a2
  mu =0;
  t = 30*10^-3;%20ms
  %t=0.1;
%  sigma = (2*D*t)^(1/2);

%% These are the diffusion coefficients from the interpolation
sigmaLin40=3.2*10^5; %points 1 to 5
sigmaLin40later=2.7*10^5; %points 3 to 8

sigmaLin50=2*10^5; %points 1 to 5
sigmaLin50later=1.8*10^5; %points 3 to 8

sigmaLin15= 6*10^4; %points 1 to 5
sigmaLin15later= 5.2*10^4; %points 3 to 8

sigmaLin18= 1.3*10^5; %points 1 to 5
sigmaLin18later= 1.1*10^5; %points 3 to 8

sigmaLin20= 2.2*10^5; %points 1 to 5
sigmaLin20later= 1.5*10^5; %points 3 to 8

%%
 sigma = (sigmaLin50/3*t)^(1/2);
%sigmaLin20Temp = 2.75*10^5;

%sigma = (sigmaLin20Temp/3*t)^(1/2);
 %R = 280; %radius of nucleoid
 %hl=700; %half length nucleoid
 R=400%330
 hl=950%775%775
 
 thickness =100%60;%50;%R%250; %thickness of ring
 Rmin=R-thickness;
 r0 = randInterval (Rmin,R,1000,1);
 th0 = randInterval (0,2*pi(),1000,1);
 x0 = randInterval (-hl,hl,1000,1);
 y0= r0.*sin(th0);
 z0=r0.*cos(th0);
 %figure
 %plot(x0,y0,'.')
 %plot3(x0,y0,z0,'.')
 simNum=floor(randInterval(1,1000,1,1));
 counter =0;
 
 tsteps =1000;
 x=zeros(tsteps,1);
 y=zeros(tsteps,1);
 z=zeros(tsteps,1);
 disp=normrnd(mu,sigma,tsteps,3);
 counter =0;
 for iiStep = 1:tsteps
     iiStep;
     if (iiStep == 1)
        x(iiStep)=x0(simNum);
        y(iiStep)=y0(simNum);
        z(iiStep)=z0(simNum);
     else
        
        xtemp = x(iiStep-1)+disp(iiStep,1);
        ytemp = y(iiStep-1)+disp(iiStep,2);
        ztemp = z(iiStep-1)+disp(iiStep,3);
        rt= sqrt(ytemp^2+ztemp^2);
        ht=abs(xtemp);
        counter =0;
        
        rt= sqrt(ytemp^2+ztemp^2);
        tht= atan2(ytemp,ztemp);
        if(rt>R)
            y(iiStep)=R*sin(tht);
            z(iiStep)=R*cos(tht);
            counter =counter +1;
        elseif (rt<Rmin)
            y(iiStep)=Rmin*sin(tht);
            z(iiStep)=Rmin*cos(tht);
            counter =counter+1;
            
        else
            y(iiStep)=ytemp;
            z(iiStep)=ztemp;
        end
        if(abs(xtemp)>hl)
            x(iiStep)= sign(xtemp)*hl;
        else
            x(iiStep)=xtemp;
        end
        
        
        
        
        
        
%         while (rt>R|rt<Rmin|ht>hl)
%             temp=normrnd(mu,sigma,1,3);
%             xtemp= x(iiStep-1) + temp(1);
%             ytemp= y(iiStep-1) + temp(2);
%             ztemp= z(iiStep-1) + temp(3);
%             rt= sqrt(ytemp^2+ztemp^2);
%             ht=abs(xtemp);
%             counter =counter +1;
%         
%         end
%         
%         x(iiStep) =xtemp;
%         y(iiStep)= ytemp;
%         z(iiStep) = ztemp;

    end
        
 end
 
Xprime =x;%/80
Yprime =y;%/80
Zprime =z;%(z/100)*0.59

figure
plot3(Xprime,Yprime,Zprime)
%figure
%plot(Yprime,Zprime)





