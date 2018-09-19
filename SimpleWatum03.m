clc , clear ;
length_river = 1000;
dx      = 1 ;
decay        = 0 ;
q            = 1 ;
b            = 1 ;
y            = 1 ;
z            = 0 ;
s            = 0.0001;
A            = b*y ;
u            = q/A ;
vol          = A*dx;
Alpha        = 1; 
Betha        = 1 - Alpha ;
%% LDC maker
% chossing LDC
    disp('.* is right');
    disp('you can chosse your desired LDC by the List which is shown below');
    disp('pleas insert the number which is assigned to your choice...');
    disp(' (1) Fischer 1975'                );
    disp(' (2) Seo & Cheong 1998'           );
    disp(' (3) Deng 2001'                   );
    disp(' (4) Kashefipour & Falconer 2002' );
    disp(' (5) Rajeev and Dutta 2009'       );
    disp(' (6) Jhang 2015'                  );
    disp(' (7) Qiasi  2015 '                );
    disp(' (8) if you want study on a not dispersive model '                );
    
       LCD_num = input('Enter a number: ');
      filename = input('and gimme a filename to save your work ','s');

if     LCD_num == 1
        LCD0 = fischer(b,y,z,s,q) ;
elseif LCD_num == 2
        LCD0 = seo(b,y,z,s,q) ;
elseif LCD_num == 3
        LCD0 = deng(b,y,z,s,q) ;
elseif LCD_num == 4
        LCD0 = kashefipour(b,y,z,s,q) ;
elseif LCD_num == 5
        LCD0 = rajeev(b,y,z,s,q) ;
elseif LCD_num == 6
        LCD0 = jhang(b,y,z,s,q) ;
elseif LCD_num == 7
        LCD0 = qiasi(b,y,z,s,q) ;
elseif LCD_num == 8
        LCD0 = 0 ;        
end 

epr          = LCD0*A/dx  ;

%end LDC maker

total_time   = floor(length_river/u)+1 ;
dt      = ((dx ^ 2) / (( u * dx * (Alpha - Betha) + 2 * LCD0 + decay * dx^2)));%delta_X/u;%
nT           = floor(total_time/dt)+1;
nX           = floor(length_river/dx)+1;
c0           = zeros(nT,nX);
co       = u*dt/dx; 
Landa        = (LCD0*dt)/(dx^2);
pause()
%%initial Condition
c1 = c0      ;
wi = 0 .* c0 ;

xx = -5: 0.1 : 5 ;
% next line is a initial condition maker 
c0( 2 , : ) = [normpdf(xx,0,1) .* 2.5 , zeros(1,length(c0(2,:))-length(xx))];


%%Boundary      
c0(:,1) = 0       ;

%Load




% matterial for saveing the results
           

mp = 1/.746 ;
etha = decay*LCD0/(u^2);
%% inform user
if etha > 1 
   disp('Diffusion predominates')
elseif etha <=1
    disp('Advection predominates')
end
formatSpec = '  Total Time number is %8.1f \n  nT number is %8.1f \n  nX number is %8.1f \n  Courant Number is %8.1f \n  Landa is %8.5f \n  delta_T is %8.2f \n  delta_X is %8.1f \n  E is %8.3f \n  E'' is %8.3f \n  Press any key to continue';
fprintf(formatSpec,total_time,nT,nX,co,Landa,dt,dx,LCD0,epr)
pause()
%% main solver
tic
for t = 2 : 1 : nT -1
      for x = 2 : 1 : nX -1
          
          c0(t+1,x) = wi(t,x)*dt/vol - dt*(-q*Alpha-epr)*c0(t,x-1)/vol + (1-(dt/vol)*(-1*q*Betha+q*Alpha+2*epr+ decay*vol))*c0(t,x) - dt*(q*Betha-epr)*c0(t,x+1)/vol ;%page 215 Chapra
          
          %c1(t+1,x) = Landa*c1(t,x-1) + (1-2*Landa)*c1(t,x) + Landa*c1(t,x+1);
          %c1(t,x)   = ( mp / ( 2 * sqrt( pi * LCD0 * t))) * exp ( -1 * ((x - u * t )^2)/(4 * LCD0 * t) - (decay * t));    
          %c0(t+1,x)   = ((u*delta_T)/(delta_X))*c0(t,x-1) + (1-u*delta_T/delta_X)*c0(t,x); %montazeri
          %c00(t, x)   = (((-1 * U * delta_T) / (2 * delta_X)) - ( (LCD0 * delta_T)/(delta_X^2)))*c00(t+1,x-1) + (1+ ((2*LCD0*delta_T)/(delta_X^2)) + decay*delta_T)*c00(t+1,x) + (((U*delta_T)/(2*delta_X)) - ((LCD0 * delta_T)/(delta_X^2)))*c00(t+1,x+1);
          %c000(t+1,x) = c000(t,x) + LCD0 * ((c000(t,x+1)-2*c000(t,x)+c000(t,x-1))*delta_T)/(delta_X^2) ;
     
      end
end
fprintf('\n  solver Time %8.1f ',toc)

%% presentation
   %  *****  %
    %  ***  %
    
steper = 1 ;
  %% statistics
    c_peak_anal = zeros(floor(nT/steper),1);
    c_peak_num  = c_peak_anal;
    
  for i = 1 :steper: nT   
      
      c_peak_anal(i,1) = max(c1(i,:));   
      c_peak_num(i,1)  = max(c0(i,:));
      
  end 
  
  fprintf('\n  statistics Time %8.1f ',toc)
  %% plotter
  grid on
  hold on
  subplot(1,4,1);plot(c1(:,1:30:end))
  grid on
  hold on
  subplot(1,4,2);plot(c0(:,1:30:end))
  grid on
  hold on
  subplot(1,4,3);plot(c_peak_anal)
  grid on
  hold on
  subplot(1,4,4);plot(c_peak_num)
  
  fprintf('\n  plotter Time %8.1f ',toc)
 %% file generator
         
 %xlswrite(filename,c0, 'c_numeric' )
 %xlswrite(filename,c1, 'c_anal' )
 xlswrite ( filename, c_peak_anal, 'peak_anal' )
 xlswrite ( filename, c_peak_num, 'peak_num' )
 
 fprintf('\n  file generatoin Time %8.1f ',toc)