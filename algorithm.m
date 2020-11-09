#Trabalho de Transmissão de Calor - EMC5403
#Nomes: Johnattan Lima Santana/Maicon Secchi

clear
N = 10;
M = 10;
l= 1; #Comprimento da Barra
dx= l/(N - 1); #Tamanho do volume de controle; dx=dy=dz
dy=dx;
h=100; #Coeficiente de Transmissão de calor
maior = 0;
dt =1000; #chute
tempo = 0;
#-------------------

	for i=1:N
		for j=1:M
			T(i,j)=0; #Temperatura Inicial
		    Ta(i,j)=0; #Temperatura Antiga Inicial
			Tc(i,j)=0;
			Tinf=0;
		endfor
	endfor
#-------------------

#Para o Material A - COBRE (BOM CONDUTOR)
#k= 401; #[W/m.K] - Condutividade do material
#ro= 8933; #[kg/m³] - Massa específica do material 
#cp= 385; #[J/kg.K] - Calor específico a 300k do material

#Para o Material B - Tijolo de argila (BOA APLICAÇÃO)
k= 1; #[W/m.K] - Condutividade do material
ro= 2050; #[kg/m³] - Massa específica do material 
cp= 960; #[J/kg.K] - Calor específico a 300k do material

#Para o Material C - Poliestireno Expandido - Pérolas Moldadas (BOM ISOLANTE)
#k= 0.04; #[W/m.k] - Condutividade do material
#ro= 16; #[kg/m³] - Massa específica do material 
#cp= 1210; #[J/kg.K] - Calor específico a 300k do material

#Para o Material D - Aluminio (BOM CONDUTOR)
#k= 237; #[W/m.K] - Condutividade do material
#ro= 2702; #[kg/m³] - Massa específica do material 
#cp= 903; #[J/kg.K] - Calor específico a 300k do material

bi=(h*dx)/k #Bi do problema
alfa = k/(ro*cp) #Alfa
fo= alfa*dt/(dx*dx); #Fourier

#	RESOLVENDO AS TEMPERATURAS
niteracoes=500
	for p=1:niteracoes
		breakout = 0;
			while breakout == 0 #	REFINAMENTO DAS TEMPERATURA - LAÇO INTERNO
#------------------- (QUINA SUPERIOR ESQUERDA)

	T(1,1)=(2*fo*T(1+1,1) +2*fo*T(1,1+1) +Ta(1,1) +(4*fo*bi*0))/(4*fo*bi+4*fo+1);
#------------------- (BORDA ESQUERDA)
	for i=2:N-1
		T(i,1)=(fo*T(i-1,1) +2*fo*T(i,1+1) +fo*T(i+1,1) +Ta(i,1) +(2*bi*fo*0))/(1+4*fo+2*bi*fo);
	endfor
#-------------------- (QUINA INFERIOR ESQUERDA)

		T(N,1)=(2*fo*T(N-1,1) +2*fo*T(N,1+1) +Ta(N,1) +(4*fo*bi*0))/(4*fo*bi+4*fo+1);
#-------------------- (BORDA SUPERIOR)
	for j=2:(M-1)			
		g=0;
		T(1,j)=(fo*T(1,j-1) +2*fo*T(1+1,j) +fo*T(1,j+1) +(2*bi*fo*g) +Ta(1,j))/(1+4*fo+2*bi*fo);
	endfor
#--------------------- (MEIO)
	for i=2:(N-1)
		for j=2:(M-1)	
			T(i,j)=(fo*T(i,j-1) +fo*T(i+1,j) +fo*T(i,j+1) +fo*T(i-1,j) +Ta(i,j))/(1+4*fo);					
		endfor
	endfor
#---------------------- (BORDA INFERIOR)
	for j=2:(M-1)
		T(N,j)=(fo*T(N,j-1) +2*fo*T(N-1,j) +fo*T(N,j+1) +Ta(N,j) +(2*bi*fo*0))/(1+4*fo+2*bi*fo);
	endfor
#---------------------- (QUINA SUPERIOR DIREITA)

		T(1,M)=(2*fo*T(1+1,M) +2*fo*T(1,M-1) +Ta(1,M) +(4*fo*bi*0))/(4*fo*bi+4*fo+1);		
#---------------------- (BORDA DIREITA)
	for i=2:N-1
		T(i,M)=(fo*T(i-1,M) +2*fo*T(i,M-1) +fo*T(i+1,M) +Ta(i,M) +(2*bi*fo*0))/(1+4*fo+2*bi*fo);
	endfor
#----------------------- (QUINA INFERIOR DIREITA)
	
		T(N,M)=(2*fo*T(N-1,M) +2*fo*T(N,M-1) +Ta(N,M) +(4*fo*bi*0))/(4*fo*bi+4*fo+1);	
#-----------------------
	
#	CALCULO DO ERRO
maior = 0;

for i=1:N
	for j=1:M
		erro = abs(T(i,j)-Tc(i,j));
		if erro > maior
			maior = erro;
		endif
	endfor
endfor

if maior < 10^-4
	breakout = 1;
end

if maior > 10^-4
	breakout = 0;
end	

for i=1:N
	for j=1:M
		Tc(i,j) = T(i,j);
	end
end


endwhile

# TEMPERATURA ANTIGA RECEBE TEMPERATURA NOVA
for i=1:N
	for j=1:M
		Ta(i,j) = T(i,j);
	endfor
endfor

# DADOS PARA CONSTRUÇÃO DOS GRÁFICOS
tempo = tempo + dt;
p;
Vt(p)=tempo;
Vtgr(p)=T(1,5);
Vtgr1(p)=T(5,5);
Vtgr2(p)=T(1,1);
Vtgr3(p)=T(N,1);


for i=1:N
	for j=1:M
		T(i,j); #MOSTRA TEMPERATURAS PARA CADA P
	endfor
endfor

endfor

for i=1:N
	for j=1:M
		T(i,j) #MOSTRA TEMPERATURAS PARA ÚLTIMO P
	endfor
endfor

#TAXA NAS QUINAS SUPERIORES, LATERAIS E INFERIORES (ESQUERDA E DIREITA)
	g=0;
	qse =(g-T(1,1))*h*(dx/2)
	qsd =(g-T(1,M))*h*(dx/2)
	qles =(0-T(1,1))*h*(dy/2)
	qlds =(0-T(1,M))*h*(dy/2)
	qlei=(0-T(N,1))*h*(dy/2)
	qldi=(0-T(N,M))*h*(dy/2)
	qie =(0-T(N,1))*h*(dx/2)
	qid=(0-T(N,M))*h*(dx/2)
quina=qse+qsd+qles+qlds+qlei+qldi+qie+qid
#TAXA NA BORDA SUPERIOR
soma1=0;
for j=2:M-1
	g=100*sin((dx/2)*(2*j-2)*pi/l);
	q1 =(g-T(1,j))*h*(dx)
	soma1 = soma1 + q1;
endfor
q1in=soma1;
q1in

#TAXA NA BORDA ESQUERDA
soma2=0;
for i=2:N-1
	q2 =(0-T(i,1))*h*(dx)
	soma2 = soma2 + q2;
endfor
q2out=soma2;
q2out

#TAXA NA BORDA DIREITA
soma3=0;
for i=2:N-1
	q3 =(0-T(i,M))*h*(dx)
	soma3 = soma3 + q3;
endfor
q3out=soma3;
q3out

#TAXA NA BORDA INFERIOR
soma4=0;
for j=1:M-1
	q4 =(0-T(N,j))*h*(dx)
	soma4 = soma4 + q4;
endfor
q4out=soma4;
q4out

tempo
qintotal=q1in
qouttotal=quina+q2out+q3out+q4out
errop = qintotal-qouttotal;



plot(Vt,Vtgr,"-r;Ponto médio superior;", Vt,Vtgr1,"-b;Ponto médio;", Vt,Vtgr2,"-g;Quina Superior Esquerda;", Vt,Vtgr3,"-k;Quina Inferior Esquerda;")
grid on;
title('TEMPERATURA X TEMPO');
xlabel('Tempo [s]');
ylabel('Temperatura [ºC]');