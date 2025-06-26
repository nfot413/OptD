
data=importdata('F:\RAMZI_gene\fsr1.xlsx');%Ŀ�꺯��

%n_vector:���Ƶĵ缫�ĸ���the number of adjustable electrodes
%vol_range: �缫�ĵ�ѹ��Χthe range of chip voltage (V)
%vol_precision: �缫��ѹ�ľ��ȣ�С�����vol_precisionλ��the accurate decimal places of voltage source,0.01V~2
global n_vector vol_range vol_precision num 
n_vector=11;
vol_range=2;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vol_precision=1;%С���������λ��
%vol_channel=[0]%the voltage source channel number used,0~channel 1th
vol_channel=1:1:11;

W_out = zeros(1,n_vector);%parameter vector
Y_max = zeros(1, n_vector);%actual optimal output result
Y_obj = zeros(1,n_vector);%best output result

W_out_generationbest=zeros(100,n_vector);
W_out_generationbest_f=zeros(100,1);
x=0;

Y_obj = data;%Ŀ�꺯��
%---------------------------------------------------------------------------------
%Genetic algorithm optimization, using floating-point coding method
num=100;%ÿһ���ĸ�����Number of individuals per generation
length_cro=5;%Cross gene length �������ĳ��ȣ�Ҫ�������õ��㽻�滹�����������棬���㽻���ǽ���Ƭ�Σ���length_cro�йأ�������������������ϣ���length_cro�޹أ�
length_vari=2;%Variant gene length  �������ĳ��ȣ�����û�õ�
rate_cro=0.7;%Crossover probability  �������
rate_vari=0.5;%Mutation probability  �������
rate_copy=0.3;%Replication probability  ���Ƹ���
value_obj=2;%Target evaluation value  Ŀ��ƥ���
generation_max=50;%Population algebraic threshold  ������������
rate_ET=0.02;%Elite ratio  ��Ӣ��
k_rv=0;%Mutation probability adaptive factor, default 0 is no change �����������Ӧ���ӣ�����Ϊ�㣬��ʵû������Ӧ����������޹�
generation_interval=200;%Population interval for fast convergence
min_conver=1;%Minimum optimization degree of population spacing

%---------------------------------------------------------------------------------
self=species_origin1();%��һ������100�����壬ÿ��������n_vector��Ⱦɫ�壨����ͬ�ĵ�ѹ��ϣ�
%(���д��������һ����Ϊ species_origin1 �ĺ����������������ɳ�ʼ��Ⱥ�����壩��self ���������������ʼ��Ⱥ����Ϣ��������ԣ��������Ӧ�ø��𴴽��ʺ��Ż�����ĸ��壬ͨ��������ʼ������Ļ���ͼ����ʼ��Ӧ�ȣ�fitness��)
generation=1;
conver1=1000;%Initial optimization degree

max_value=zeros(1,generation_max);%Maximum fitness value of each generation
mean_value=zeros(1,generation_max);%Average fitness value of each generation
value_species=zeros(1,num);%Population fitness

%------------------�豸ͨѶ

COM_order = 4; % newsetup��4
obj2 = instrfind('Port','COM4');
if isempty(obj2)
    obj2 = SiliconExtreme(COM_order);
else
    fclose(obj2);
    obj2 = SiliconExtreme(COM_order);
end


open(obj2);
Imax=10; % mA,�������õ���󣬲�һ�����Դﵽ��ȡ��������������
Vmax=4;  % V���������õ���󣬲�һ�����Դﵽ��ȡ��������������
for i=1:n_vector
setImax(obj2,i,Imax);
setVmax(obj2,i,Vmax);
% ��һ��
setV(obj2,i,0);
end

       pause(0.5); %ֻ���²��ȣ�������ÿ�ζ����ļ��0.5s
        % Vt = getV(obj1,channel);
        % It = getI(obj1,channel);
open(obj2);

%%%%%%%64·��ѹԴͨѶ

obj1 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 01, 'Tag', '');
if isempty(obj1)
    obj1 = gpib('NI', 0, 01);
else
    fclose(obj1);
    obj1 = obj1(1);
    a=1;
end
set(obj1, 'InputBufferSize', 180009*3);
set(obj1, 'OutputBufferSize', 180009*3);
set(obj1, 'Timeout', 150);
fopen(obj1);%OSA��ͨѶ
%---------------

%Calculation of fitness of primary population
for i_species=1:num
    %Periodically update the voltage values of all used channels
for i=1:n_vector
setV(obj2,i,self(i_species,i));
end
        
    pause(0.01);
    [osa_Lamb, osa_Power] = OSA_AQ6370_Wavelength_sweep(obj1);%Acquire spectral data
    value_species(i_species) = f_value_20250423(osa_Lamb, osa_Power, Y_obj);
end
%pause(0.01);

while generation<=generation_max&&max(value_species)<value_obj&&(conver1>=min_conver)%�����̶�  �Ż�������ƥ���Խ��Խ��
    clear sort_num;clear sort_value;
    
    species_new=zeros(num,n_vector);%New species
    newspecies_num=1;% Number of existing individuals in the new population + 1
    
    if generation>1
        value_species=zeros(1,num);%Population fitness
        %Calculation of fitness of population
        for i_species=1:num
            %Periodically update the voltage values of all used channels

          for i=1:n_vector
     setV(obj2,i,self(i_species,i));
         end

            pause(0.01);
            [osa_Lamb, osa_Power] = OSA_AQ6370_Wavelength_sweep(obj1);%Acquire spectral data
            value_species(i_species) =  f_value_20250423(osa_Lamb, osa_Power, Y_obj);
        end
       
    end
      
        %pause(0.01);
    
    [sort_value,sort_num]=sort(value_species,'descend');%Fitness ranking
    max_value(generation)=sort_value(1);
    mean_value(generation)=mean(value_species);
    
    if generation>generation_interval%Convenient for rapid convergence
        conver1=max_value(generation)-max_value(generation-generation_interval);%Degree of optimization
    end%�����ò���
    
    generation=generation+1;%�µ�һ��
    rvn=rv_exp(rate_vari, k_rv, generation);%�ȼ��������� Calculate the adaptive mutation probability of the current generation population
   
    value_species1=zeros(1,num);
    for i = 1:num
        [p,q] = find(sort_num==i);
        value_species1(i) = value_species(sort_num(num+1-q));
    end
    value_one = value_species1./sum(value_species);%�ٰѵ�һ���ĸ�������ۺ�����һ�� normalization

    
    %Elite mechanism
    for i=1:round(rate_ET*num) %(�޸���)ѡ����һ�����ۺ�����õ� If it is not an integer, round off the decimal part
        species_new(newspecies_num,:) = self(sort_num(i),:);
        newspecies_num=newspecies_num+1;
    end
    
    while (newspecies_num-1)<num
        %Cross
        if rand<=rate_cro&&(newspecies_num+1<=num)%��Ϊ�������ǲ����������壬���������������newspecies_num+1<=num
            indiv1=indiv_selection(self, value_one);
            indiv2=indiv_selection(self, value_one);%�������̶�ѡ������ĸ��常��
        
            cross_flag=2;%���������棬���µĸ�����������������һ��Ȩ����϶��ɵ�
            %Crossing mode sign (1 refers to single point crossing, 2 refers to unconditional crossing and 3 refers to intermediate crossing)
            indiv_crossnew=crossover(indiv1, indiv2, cross_flag, length_cro);
        
            species_new(newspecies_num,:)=indiv_crossnew(1,:);
            newspecies_num=newspecies_num+1;
            species_new(newspecies_num,:)=indiv_crossnew(2,:);
            newspecies_num=newspecies_num+1;
        end
        %copy
        if rand<=rate_copy&&(newspecies_num<=num)
            species_new(newspecies_num,:)=indiv_selection(self, value_one);
            newspecies_num=newspecies_num+1;%���̶�ȷ�����Ƶ���һ������
        end
    end%���� ���ƺ��ٱ���

    for i = round(rate_ET*num)+1:num%(�޸���)
        %variation
        if rand<=rvn
            species_new(i,:)=varition(species_new(i,:));
        end
    end
 self=species_new;    
 disp(generation)
 %disp(max_value(generation));
end%�µ�һ���ǣ�����ƥ�����ߵĸ��壬�Լ����� ���� ����ĸ���

clear sort_num;clear sort_value;

value_species=zeros(1,num);

for i_species=1:num
for i=1:n_vector
setV(obj2,i,self(i_species,i));
end

    pause(0.01);
    [osa_Lamb, osa_Power] = OSA_AQ6370_Wavelength_sweep(obj1);
    value_species(i_species) = f_value_20250423(osa_Lamb, osa_Power, Y_obj);
end
%pause(0.01);


[sort_value,sort_num]=sort(value_species,'descend');
max_value(generation)=sort_value(1);
mean_value(generation)=mean(value_species);

W_out = self(sort_num(1),:);%The best individual
%writematrix(W_out,'GA_W/best_indiv.csv');
%disp("The optimal weight vector is:");
%disp(W_out);
%W_out;
%fprintf("The fitness of the optimal vector is:%f",sort_value(1));

for i=1:n_vector
setV(obj2,i,W_out(i));
end

pause(0.01);
[osa_Lamb, osa_Power] = OSA_AQ6370_Wavelength_sweep(obj1);

fclose(obj1);
delete(obj1);
% plot(Y_obj);title('best output result');
% figure;
%plot(Y_max);title('actual optimal output result');
plot(max_value(1:generation));title("Optimization process -- optimal fitness of each generation");
figure;
plot(mean_value(1:generation));title("Optimization process -- average fitness of each generation");



function [osa_Lamb, osa_Power] = OSA_AQ6370_Wavelength_sweep(obj1)
clc;
format long

lamstar = 1550;
lamstop = 1551;% optical source

res =0.02;
% res =0.1;
sensitivity = 'normal';
%sensitivity = 'high2';
pointnum = (lamstop-lamstar)*5/res+1;%采样点数�?
speed = '2x';

reflevel = 0;
channel = 'a';
r=1;% 1:制定波长范围�?0：默认光源范�?5
%---------------------------------------------------------------------------------

% Communicating with instrument object, obj1.
%-----------------------------------------------------------------

fprintf(obj1, [':sens:wav:star ',num2str(lamstar),'nm']);
fprintf(obj1, [':sens:wav:stop ',num2str(lamstop),'nm']);
fprintf(obj1, [':disp:trace:y1:rlev ',num2str(reflevel),'dbm']);
fprintf(obj1,[':sens:swe:points ',num2str(pointnum)]);
% fprintf(obj1,':sens:swe:points:auto on');
fprintf(obj1,[':sens:band:res ',num2str(res),'nm']);
fprintf(obj1,[':sens:sens ',sensitivity]);
fprintf(obj1,[':sens:swe:spe ',speed]);
fprintf(obj1,[':trac:stat:tr',channel,' on']);
fprintf(obj1,[':trac:attr:tr',channel,' write']);
fprintf(obj1, ':initiate:smode single');
fprintf(obj1,'*CLS');
fprintf(obj1,':init');
fprintf(obj1,':stat:oper:even?');
finish = fscanf(obj1);
while ~finish
    fprintf(obj1,':stat:oper:even?')
    fscanf(obj1);
end
fprintf(obj1,':format:data ascii');
fprintf(obj1,[':trac:y? tr',channel]);
Power = fscanf(obj1);
fprintf(obj1,[':trac:x? tr',channel]);
Lambda = fscanf(obj1);

% Preset the parameter for alignment
%-------------------------------------------------------------
if r==1
    fprintf(obj1, [':sens:wav:star ',num2str(lamstar),'nm']);
    fprintf(obj1, [':sens:wav:stop ',num2str(lamstop),'nm']);
elseif r==0
    fprintf(obj1, [':sens:wav:star ','1510','nm']);
    fprintf(obj1, [':sens:wav:stop ','1610','nm']);
end
fprintf(obj1, ':disp:trace:y1:rlev -20dbm'); %%% 这里是OSA上面显示的ref
fprintf(obj1,':sens:band:res 1nm');
fprintf(obj1,':sens:sens norm');
fprintf(obj1,':sens:swe:points:auto on');
fprintf(obj1,':sens:swe:spe 2x');
fprintf(obj1, ':initiate:smode repeat');
fprintf(obj1, ':init');


% Dealing with the data
%----------------------------------------------------------
Points = pointnum;
Power = char(Power);
Lambda = char(Lambda);
P = zeros(1,Points);
L = zeros(1,Points);
I = strfind(Power,',');
P(1) = str2double(Power(1:I(1)-1));
for Loop = 1:length(I)-1
    P(Loop+1) = str2double(Power(I(Loop)+1:I(Loop+1)-1));
end
P(end) = str2double(Power(I(end)+1:end));
I = strfind(Lambda,',');
L(1) = str2double(Lambda(1:I(1)-1));
for Loop = 1:length(I)-1
    L(Loop+1) = str2double(Lambda(I(Loop)+1:I(Loop+1)-1));
end
L(end) = str2double(Lambda(I(end)+1:end));
L=L*1e9;%1*pointnum的行向量
osa_Lamb = L;
osa_Power = P;%���δ������Ҫ�����ǽ����ַ�����ʽ�Ĺ��ʺͲ������ݣ�����ת��Ϊ��ֵ���飬�����е�λת��

% Plot the data
%------------------------------------------------------------
%figure(1);
%hold on;
%plot(L,P,'k');
%grid on;
%axis ([lamstar lamstop -80 reflevel]);
end
 
%%%%%%%%%%%%%%%%%%% ���� osa �õ�������

%%%%%%%%%%%%%%%%%%%