%Figure4.m
clear
close all
clc
here = pwd;
labelsize=.02;

f1 = figure('Position',[20 500 1000 500]);

N =70 ;

%% Load quadrant cell results
load('./quadrantcells.mat');
A = arrayfun(@(x) abs(cellArea(tri2cell(x.modelshape))),quadrantcells);
D = [quadrantcells.divisionmode];
L = [quadrantcells.l1; quadrantcells.l2; quadrantcells.l3];
[junk mode1] = min(L,[],1);
[junk mode3] = max(L,[],1);

L_sorted = sort(L,1);
delta12 = ((L_sorted(2,:)- L_sorted(1,:))./sqrt(A));

D1 = (D==mode1);
display('Quadrant cells');
display(['Number of cells: ' num2str(length(D1))]);
display(['Percentage of mode 1: ' num2str(sum(D1)/length(D1))]);
D2 = (D~=mode3 & D~=mode1);
%D3 = (D==mode3);

delta12_QC= delta12(D1|D2);
D12_QC = D1(D1|D2);


%% Load Dionaea results
load('./dionaea_results.mat');
A = arrayfun(@(x) abs(cellArea(tri2cell(x.modelshape,2*pi/3))),dionaea_results);
D = [dionaea_results.divisionmode];
mode1 = [dionaea_results.mode1];
mode2 = [dionaea_results.mode2];
mode3 = [dionaea_results.mode3];
D1 = (D==mode1);
D2 = (D==mode2);
D3 = (D==mode3);
display('Triangular cells');
display(['Number of cells: ' num2str(length(D1))]);
display(['Percentage of mode 1: ' num2str(sum(D1)/length(D1))]);

l1 = [dionaea_results.l1];
l2 = [dionaea_results.l2];
l3 = [dionaea_results.l3];

delta12 = (l2-l1)./sqrt(A);
delta13 = (l3-l1)./sqrt(A);
delta23 = (l3-l2)./sqrt(A);

% Equilibriun solution 1- solution2

delta12_TC= delta12(D1|D2);
D12_TC = D1(D1|D2);
[delta12 index] = sort([delta12_QC delta12_TC]);
D12 = [D12_QC D12_TC];
D12 =D12(index);

delta12_Dionaea = zeros(floor(length(delta12)/N),1);
ddelta12_Dionaea = zeros(size(delta12_Dionaea));
r12_Dionaea = zeros(size(delta12_Dionaea));
dr12_Dionaea = zeros(size(delta12_Dionaea));
for i =1:floor(length(delta12)/N)
    delta12_Dionaea(i) = mean(delta12((i-1)*N+1:i*N));
    ddelta12_Dionaea(i) = std(delta12((i-1)*N+1:i*N));
    r12_Dionaea(i) = mean(D12((i-1)*N+1:i*N));
end

% Equilibriun solution 1- solution3
delta13= delta13(D1|D3);
D13 = D1(D1|D3);

[delta13 index] = sort(delta13);
D13 =D13(index);
delta13_Dionaea = zeros(floor(length(delta13)/N),1);
ddelta13_Dionaea = zeros(size(delta13_Dionaea));
r13_Dionaea = zeros(size(delta13_Dionaea));
dr13_Dionaea = zeros(size(delta13_Dionaea));
for i =1:floor(length(delta13)/N)
    delta13_Dionaea(i) = mean(delta13((i-1)*N+1:i*N));
    ddelta13_Dionaea(i) = std(delta13((i-1)*N+1:i*N));
    r13_Dionaea(i) = mean(D13((i-1)*N+1:i*N));
end

% Equilibriun solution 2- solution3
delta23= delta23(D2|D3);
D23 = D2(D2|D3);

[delta23 index] = sort(delta23);
D23 =D23(index);
delta23_Dionaea = zeros(floor(length(delta23)/N),1);
ddelta23_Dionaea = zeros(size(delta23_Dionaea));
r23_Dionaea = zeros(size(delta23_Dionaea));
dr23_Dionaea = zeros(size(delta23_Dionaea));
for i =1:floor(length(delta23)/N)
    delta23_Dionaea(i) = mean(delta23((i-1)*N+1:i*N));
    ddelta23_Dionaea(i) = std(delta23((i-1)*N+1:i*N));
    r23_Dionaea(i) = mean(D23((i-1)*N+1:i*N));
end


delta12_quadrantcells = zeros(floor(length(delta12)/N),1);
ddelta12_quadrantcells = zeros(size(delta12_quadrantcells));
r12_quadrantcells = zeros(size(delta12_quadrantcells));
dr12_quadrantcells = zeros(size(delta12_quadrantcells));
for i =1:floor(length(delta12)/N)
    delta12_quadrantcells(i) = mean(delta12((i-1)*N+1:i*N));
    ddelta12_quadrantcells(i) = std(delta12((i-1)*N+1:i*N));
    r12_quadrantcells(i) = mean(D12((i-1)*N+1:i*N));
end

%% Load Coleochaete results
load('./coleochaeteresults.mat');

A = arrayfun(@(x) abs(cellArea(col2cell(x.modelshape))),coleochaeteresults);
D = [coleochaeteresults.divisionmode];

mode1 = arrayfun(@(x) x.modes(1),coleochaeteresults);
mode2 = arrayfun(@(x) x.modes(2),coleochaeteresults);
nmodes = arrayfun(@(x) length(x.modes),coleochaeteresults);
mode3(nmodes~=2) = arrayfun(@(x) x.modes(3),coleochaeteresults(nmodes~=2));
mode3(nmodes==2) = NaN;
%mode3 = arrayfun(@(x) x.modes(3),coleochaeteresults);
%mode1 = [coleochaeteresults.mode1];
%mode2 = [coleochaeteresults.mode2];
l1 = arrayfun(@(x) x.planes(1).walllength,coleochaeteresults);
l2 = arrayfun(@(x) x.planes(2).walllength,coleochaeteresults);
l3(nmodes~=2) = arrayfun(@(x) x.planes(3).walllength,coleochaeteresults(nmodes~=2));
l3(nmodes==2) = NaN;

delta12 = (l2-l1)./sqrt(A);
delta13 = (l3-l1)./sqrt(A);
delta23 = (l3-l2)./sqrt(A);
D1 = (D==mode1);
D2 = (D==mode2);
D3 = (D==mode3);
D3(isnan(D3))=0;
display('Coleochaete cells');
display(['Number of cells: ' num2str(length(D1))]);
display(['Percentage of mode 1: ' num2str(sum(D1)/length(D1))]);

delta12= delta12(D1|D2);
D12 = D1(D1|D2)./(D1(D1|D2)+D2(D1|D2));  % Because we cannot differentiate between the anticlinal divisions

[delta12 index] = sort(delta12);
D12 =D12(index);
delta12_Coleochaete = zeros(floor(length(delta12)/N),1);
ddelta12_Coleochaete = zeros(size(delta12_Coleochaete));
r12_Coleochaete = zeros(size(delta12_Coleochaete));
dr12_Coleochaete = zeros(size(delta12_Coleochaete));
for i =1:floor(length(delta12)/N)
    delta12_Coleochaete(i) = mean(delta12((i-1)*N+1:i*N));
    ddelta12_Coleochaete(i) = std(delta12((i-1)*N+1:i*N));
    r12_Coleochaete(i) = mean(D12((i-1)*N+1:i*N));
end

% Equilibriun solution 1- solution3
delta13= delta13(D1|D3);
D13 = D1(D1|D3)./(D1(D1|D3)+D3(D1|D3));

[delta13 index] = sort(delta13);
D13 =D13(index);
delta13_Coleochaete = zeros(floor(length(delta13)/N),1);
ddelta13_Coleochaete = zeros(size(delta13_Coleochaete));
r13_Coleochaete = zeros(size(delta13_Coleochaete));
dr13_Coleochaete = zeros(size(delta13_Coleochaete));
for i =1:floor(length(delta13)/N)
    delta13_Coleochaete(i) = mean(delta13((i-1)*N+1:i*N));
    ddelta13_Coleochaete(i) = std(delta13((i-1)*N+1:i*N));
    r13_Coleochaete(i) = mean(D13((i-1)*N+1:i*N));
end
% 
% N=40
% % Equilibriun solution 1- solution3
% delta23= delta23(D2|D3);
% D23 = D2(D2|D3)./(D2(D2|D3)+D3(D2|D3));
% 
% [delta23 index] = sort(delta23);
% D23 =D23(index);
% delta23_Coleochaete = zeros(floor(length(delta23)/N),1);
% ddelta23_Coleochaete = zeros(size(delta23_Coleochaete));
% r23_Coleochaete = zeros(size(delta23_Coleochaete));
% dr23_Coleochaete = zeros(size(delta23_Coleochaete));
% for i =1:floor(length(delta23)/N)
%     delta23_Coleochaete(i) = mean(delta23((i-1)*N+1:i*N));
%     ddelta23_Coleochaete(i) = std(delta23((i-1)*N+1:i*N));
%     r23_Coleochaete(i) = mean(D23((i-1)*N+1:i*N));
% end

%% Load Zinnia results

load('./Zinnia1_5_22_A_ssd_010-treated.mat');
tissue = simplifiedtissue;

ncells = length(originaltissue.c)-length(tissue.c); 
real_edges = sort(edges_list,2);

divisionmodes = NaN(1,ncells);

for i=1:ncells
    ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)';
    [a j] = intersect(ideal_edges,real_edges(i,:),'rows');
    if ~isempty(j)
        divisionmodes(i)=j;
    end
end

A =  abs(cellArea(tissue,1:ncells))';
D1 = (divisionmodes==1);
D2 = (divisionmodes==2);
D3 = (divisionmodes==3);
D_Zinnia = [sum(D1) sum(D2) sum(D3) length(D)-sum(D1)-sum(D2)-sum(D3)];
I =  arrayfun(@(x) isempty(x.walllength),planes(:,3))';
display('Zinnia cells');
display(['Number of cells: ' num2str(length(D1))]);
display(['Percentage of mode 1: ' num2str(sum(D1)/length(D1))]);

l1 = [planes(:,1).walllength];
l2 = [planes(:,2).walllength];
l3 = [planes(:,3).walllength];
delta12 = (l2-l1)./sqrt(A);
delta13 = (l3-l1(~I))./sqrt(A(~I));
delta23 = (l3-l2(~I))./sqrt(A(~I));

 
delta12= delta12(D1|D2);
D12 = D1(D1|D2);

[delta12 index] = sort(delta12);
D12 =D12(index);
delta12_Zinnia = zeros(floor(length(delta12)/N),1);
ddelta12_Zinnia = zeros(size(delta12_Zinnia));
r12_Zinnia = zeros(size(delta12_Zinnia));
dr12_Zinnia = zeros(size(delta12_Zinnia));
for i =1:floor(length(delta12)/N)
    delta12_Zinnia(i) = mean(delta12((i-1)*N+1:i*N));
    ddelta12_Zinnia(i) = std(delta12((i-1)*N+1:i*N));
    r12_Zinnia(i) = mean(D12((i-1)*N+1:i*N));
end

delta13= delta13(D1(~I)|D3(~I));
d1 = D1(~I);
D13 = d1(D1(~I)|D3(~I));

[delta13 index] = sort(delta13);
D13 =D13(index);
delta13_Zinnia = zeros(floor(length(delta13)/N),1);
ddelta13_Zinnia = zeros(size(delta13_Zinnia));
r13_Zinnia = zeros(size(delta13_Zinnia));
dr13_Zinnia = zeros(size(delta13_Zinnia));
for i =1:floor(length(delta13)/N)
    delta13_Zinnia(i) = mean(delta13((i-1)*N+1:i*N));
    ddelta13_Zinnia(i) = std(delta13((i-1)*N+1:i*N));
    r13_Zinnia(i) = mean(D13((i-1)*N+1:i*N));
end

delta23= delta23(D2(~I)|D3(~I));
d2 = D2(~I);
D23 = d2(D2(~I)|D3(~I));

N=70;
[delta23 index] = sort(delta23);
D23 =D23(index);
delta23_Zinnia = zeros(floor(length(delta23)/N),1);
ddelta23_Zinnia = zeros(size(delta23_Zinnia));
r23_Zinnia = zeros(size(delta23_Zinnia));
dr23_Zinnia = zeros(size(delta23_Zinnia));
for i =1:floor(length(delta23)/N)
    delta23_Zinnia(i) = mean(delta23((i-1)*N+1:i*N));
    ddelta23_Zinnia(i) = std(delta23((i-1)*N+1:i*N));
    r23_Zinnia(i) = mean(D23((i-1)*N+1:i*N));
end

%% Load Fern results

A=[];
divisionmodes = [];
I=[];
l1=[];
l2=[];
l3=[];
for i=1:5
    load(['./Microsorum' num2str(i) '_1-treated.mat']);
    tissue = simplifiedtissue;

    ncells = length(originaltissue.c)-length(tissue.c); 
    real_edges = sort(edges_list,2);

    %cell_list = 1:ncells;
    D = NaN(1,ncells);

    for i=1:ncells
        ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)';
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows');
        if ~isempty(j)
            D(i)=j;
        end
    end
    divisionmodes=[divisionmodes D];
    A =  [A abs(cellArea(tissue,1:ncells))'];

    I =  [I arrayfun(@(x) isempty(x.walllength),planes(:,3))'];

    l1 = [l1 [planes(:,1).walllength]];
    l2 = [l2 [planes(:,2).walllength]];
    l3 = [l3 [planes(:,3).walllength]];
end

D1 = (divisionmodes==1);
D2 = (divisionmodes==2);
D3 = (divisionmodes==3);
delta12 = (l2-l1)./sqrt(A);
delta13 = (l3-l1(~I))./sqrt(A(~I));
delta23 = (l3-l2(~I))./sqrt(A(~I));

display('Microsorum cells');
display(['Number of cells: ' num2str(length(D1))]);
display(['Percentage of mode 1: ' num2str(sum(D1)/length(D1))]);

delta12= delta12(D1|D2);
D12 = D1(D1|D2);

[delta12 index] = sort(delta12);
D12 =D12(index);
delta12_Microsorum = zeros(floor(length(delta12)/N),1);
ddelta12_Microsorum = zeros(size(delta12_Microsorum));
r12_Microsorum = zeros(size(delta12_Microsorum));
dr12_Microsorum = zeros(size(delta12_Microsorum));
for i =1:floor(length(delta12)/N)
    delta12_Microsorum(i) = mean(delta12((i-1)*N+1:i*N));
    ddelta12_Microsorum(i) = std(delta12((i-1)*N+1:i*N));
    r12_Microsorum(i) = mean(D12((i-1)*N+1:i*N));
end

delta13= delta13(D1(~I)|D3(~I));
d1 = D1(~I);
D13 = d1(D1(~I)|D3(~I));

[delta13 index] = sort(delta13);
D13 =D13(index);
delta13_Microsorum = zeros(floor(length(delta13)/N),1);
ddelta13_Microsorum = zeros(size(delta13_Microsorum));
r13_Microsorum = zeros(size(delta13_Microsorum));
dr13_Microsorum = zeros(size(delta13_Microsorum));
for i =1:floor(length(delta13)/N)
    delta13_Microsorum(i) = mean(delta13((i-1)*N+1:i*N));
    ddelta13_Microsorum(i) = std(delta13((i-1)*N+1:i*N));
    r13_Microsorum(i) = mean(D13((i-1)*N+1:i*N));
end

delta23= delta23(D2(~I)|D3(~I));
d2 = D2(~I);
D23 = d2(D2(~I)|D3(~I));

N=70;
[delta23 index] = sort(delta23);
D23 =D23(index);
delta23_Microsorum = zeros(floor(length(delta23)/N),1);
ddelta23_Microsorum = zeros(size(delta23_Microsorum));
r23_Microsorum = zeros(size(delta23_Microsorum));
dr23_Microsorum = zeros(size(delta23_Microsorum));
for i =1:floor(length(delta23)/N)
    delta23_Microsorum(i) = mean(delta23((i-1)*N+1:i*N));
    ddelta23_Microsorum(i) = std(delta23((i-1)*N+1:i*N));
    r23_Microsorum(i) = mean(D23((i-1)*N+1:i*N));
end


axes('OuterPosition',[0 0 1/2 1]);
hold on
msize =  10;
deltamax = .5;
index = delta12_Coleochaete<deltamax;
errorbar2(delta12_Coleochaete(index),ddelta12_Coleochaete(index),r12_Coleochaete(index),dr12_Coleochaete(index),...
    'sr','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize);
index = delta13_Coleochaete<deltamax;
errorbar2(delta13_Coleochaete(index),ddelta13_Coleochaete(index),r13_Coleochaete(index),dr13_Coleochaete(index),...
    'sr','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize);
errorbar2(delta12_Dionaea,ddelta12_Dionaea,r12_Dionaea,dr12_Dionaea,...
    '^b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize);
errorbar2(delta13_Dionaea,ddelta13_Dionaea,r13_Dionaea,dr13_Dionaea,...
    '^b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize);
errorbar2(delta23_Dionaea,ddelta23_Dionaea,r23_Dionaea,dr23_Dionaea,...
    '^b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize);
errorbar2(delta12_Zinnia,ddelta12_Zinnia,r12_Zinnia,dr12_Zinnia,...
    'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
index = delta13_Zinnia<deltamax;
errorbar2(delta13_Zinnia(index),ddelta13_Zinnia(index),r13_Zinnia(index),dr13_Zinnia(index),...
    'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
errorbar2(delta23_Zinnia,ddelta23_Zinnia,r23_Zinnia,dr23_Zinnia,...
    'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
errorbar2(delta12_Microsorum,ddelta12_Microsorum,r12_Microsorum,dr12_Microsorum,...
    'dg','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize);
index = delta13_Microsorum<deltamax;
errorbar2(delta13_Microsorum(index),ddelta13_Microsorum(index),r13_Microsorum(index),dr13_Microsorum(index),...
    'dg','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize);
errorbar2(delta23_Microsorum,ddelta23_Microsorum,r23_Microsorum,dr23_Microsorum,...
    'dg','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize);
%plot([.08 .08],[.4 1],'--k','Linewidth',1);
%plot([.15 .15],[.4 1],'--k','Linewidth',1);
axis([0 .5 .4 1.05])
set(gca,'FontName','Helvetica','FontSize',16)
xlabel('\delta_{ij}','FontSize',20,'FontName','Helvetica');
ylabel('n_i/(n_i+n_j)','FontSize',20,'FontName','Helvetica');

set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
box on
%% Fit
fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint',.1);
sigmfit = fittype('1/(1+exp(-alpha*x))','options',fit_options);

% Fit both series of data
% deltatot = [delta12_Coleochaete;delta12_quadrantcells;delta12_Dionaea;delta13_Dionaea;delta23_Dionaea;...
%     delta12_Zinnia;delta13_Zinnia;delta23_Zinnia;delta12_Microsorum;delta13_Microsorum;delta23_Microsorum];
% rtot = [r12_Coleochaete;r12_quadrantcells;r12_Dionaea;r13_Dionaea;r23_Dionaea;...
%     r12_Zinnia;r13_Zinnia;r23_Zinnia;r12_Microsorum;r13_Microsorum;r23_Microsorum];
delta_Coleochaete = delta12_Coleochaete;
delta_Dionaea = [delta12_quadrantcells;delta12_Dionaea;delta13_Dionaea;delta23_Dionaea];
delta_Zinnia = [delta12_Zinnia;delta13_Zinnia;delta23_Zinnia];
delta_Microsorum = [delta12_Microsorum;delta13_Microsorum;delta23_Microsorum];
r_Coleochaete = r12_Coleochaete;
r_Dionaea = [r12_quadrantcells;r12_Dionaea;r13_Dionaea;r23_Dionaea];
r_Zinnia = [r12_Zinnia;r13_Zinnia;r23_Zinnia];
r_Microsorum = [r12_Microsorum;r13_Microsorum;r23_Microsorum];
deltatot = [delta_Coleochaete;delta_Dionaea;delta_Zinnia;delta_Microsorum];
rtot = [r_Coleochaete;r_Dionaea;r_Zinnia;r_Microsorum];

%[globalfit,globalfit_info] = fit(deltatot(deltatot<=.25),rtot(rtot<=.25),sigmfit);
[fit_Coleochaete] = fit(delta_Coleochaete,r_Coleochaete,sigmfit);
[fit_Dionaea] = fit(delta_Dionaea,r_Dionaea,sigmfit);
[fit_Zinnia] = fit(delta_Zinnia,r_Zinnia,sigmfit);
[fit_Microsorum] = fit(delta_Microsorum,r_Microsorum,sigmfit);
[globalfit,globalfit_info] = fit(deltatot,rtot,sigmfit);

deltafit = 0:.01:.5;
rfit = 1./(1+exp(-globalfit.alpha*deltafit));
plot(deltafit,rfit,'-k');


alpha= globalfit.alpha
N=70;
%return
%% Load quadrant cell results
load('./quadrantcells.mat');
A = arrayfun(@(x) abs(cellArea(tri2cell(x.modelshape))),quadrantcells);
D = [quadrantcells.divisionmode];
L = [quadrantcells.l1; quadrantcells.l2; quadrantcells.l3];
[junk mode1] = min(L,[],1);

L_sorted = sort(L,1);

D1_QC = (D==mode1);
P1_QC = arrayfun(@(x) 1/sum(exp(-alpha*(L_sorted(:,x)-L_sorted(1,x))/sqrt(A(x)))),1:length(L));


%% Load Dionaea results
load('./dionaea_results.mat');
A = arrayfun(@(x) abs(cellArea(tri2cell(x.modelshape,2*pi/3))),dionaea_results);
D = [dionaea_results.divisionmode];
mode1 = [dionaea_results.mode1];
mode2 = [dionaea_results.mode2];
mode3 = [dionaea_results.mode3];
D1_TC = (D==mode1);
D2 = (D==mode2);

L = [dionaea_results.l1; dionaea_results.l2; dionaea_results.l3]';

P1_TC = arrayfun(@(x) 1/sum(exp(-alpha*(L(x,:)-L(x,1))/sqrt(A(x)))),1:length(L));
P2 = arrayfun(@(x) 1/sum(exp(-alpha*(L(x,:)-L(x,2))/sqrt(A(x)))),1:length(L));

P1 = [P1_QC P1_TC];
D1 = [D1_QC D1_TC];
[P1 index] = sort(P1);
D1 =D1(index);
P1_Dionaea = zeros(floor(length(P1)/N),1);
dP1_Dionaea = zeros(size(P1_Dionaea));
r1_Dionaea = zeros(size(P1_Dionaea));
dr1_Dionaea = zeros(size(P1_Dionaea));
for i =1:floor(length(P1)/N)
    P1_Dionaea(i) = mean(P1((i-1)*N+1:i*N));
    dP1_Dionaea(i) = std(P1((i-1)*N+1:i*N));
    r1_Dionaea(i) = mean(D1((i-1)*N+1:i*N));
end

[P2 index] = sort(P2);
D2 =D2(index);
P2_Dionaea = zeros(floor(length(P2)/N),1);
dP2_Dionaea = zeros(size(P2_Dionaea));
r2_Dionaea = zeros(size(P2_Dionaea));
dr2_Dionaea = zeros(size(P2_Dionaea));
for i =1:floor(length(P2)/N)
    P2_Dionaea(i) = mean(P2((i-1)*N+1:i*N));
    dP2_Dionaea(i) = std(P2((i-1)*N+1:i*N));
    r2_Dionaea(i) = mean(D2((i-1)*N+1:i*N));
end

%% Load Coleochaete results
load('./coleochaeteresults.mat');
A = arrayfun(@(x) abs(cellArea(col2cell(x.modelshape))),coleochaeteresults);
D = [coleochaeteresults.divisionmode];
%mode1 = [coleochaeteresults.mode1];
mode1 = arrayfun(@(x) x.modes(1),coleochaeteresults);
mode2 = arrayfun(@(x) x.modes(2),coleochaeteresults);
N1 = arrayfun(@(x) sum(x.modes == x.modes(1)),coleochaeteresults);
N2 = arrayfun(@(x) sum(x.modes == x.modes(2)),coleochaeteresults);
D1 = (D==mode1)./N1;
D2 = (D==mode2)./N2;

% P1_old = arrayfun(@(x) 1/sum(exp(-2*alpha*...
%     ([coleochaeteresults(x).planes.walllength]-coleochaeteresults(x).planes(1).walllength)/sqrt(A(x)))),...
%     1:length(coleochaeteresults));
P1 = arrayfun(@(x) 1/sum(exp(-alpha*...
    ([coleochaeteresults(x).planes.walllength]-coleochaeteresults(x).planes(1).walllength)/sqrt(A(x)))),...
    1:length(coleochaeteresults));
P2 = arrayfun(@(x) 1/sum(exp(-alpha*...
    ([coleochaeteresults(x).planes.walllength]-coleochaeteresults(x).planes(2).walllength)/sqrt(A(x)))),...
    1:length(coleochaeteresults));

[P1 index] = sort(P1);
D1 =D1(index);
P1_Coleochaete = zeros(floor(length(P1)/N),1);
dP1_Coleochaete = zeros(size(P1_Coleochaete));
r1_Coleochaete = zeros(size(P1_Coleochaete));
dr1_Coleochaete = zeros(size(P1_Coleochaete));
for i =1:floor(length(P1)/N)
    P1_Coleochaete(i) = mean(P1((i-1)*N+1:i*N));
    dP1_Coleochaete(i) = std(P1((i-1)*N+1:i*N));
    r1_Coleochaete(i) = mean(D1((i-1)*N+1:i*N));
end


[P2 index] = sort(P2);
D2 =D2(index);
P2_Coleochaete = zeros(floor(length(P2)/N),1);
dP2_Coleochaete = zeros(size(P2_Coleochaete));
r2_Coleochaete = zeros(size(P2_Coleochaete));
dr2_Coleochaete = zeros(size(P2_Coleochaete));
for i =1:floor(length(P1)/N)
    P2_Coleochaete(i) = mean(P2((i-1)*N+1:i*N));
    dP2_Coleochaete(i) = std(P2((i-1)*N+1:i*N));
    r2_Coleochaete(i) = mean(D2((i-1)*N+1:i*N));
end

%% Load Zinnia results

load('./Zinnia1_5_22_A_ssd_010-treated.mat');
tissue = simplifiedtissue;

ncells = length(originaltissue.c)-length(tissue.c); 
real_edges = sort(edges_list,2);

cell_list = 1:ncells;
divisionmodes = NaN(1,ncells);

for i=1:ncells
    ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)';
    [a j] = intersect(ideal_edges,real_edges(i,:),'rows');
    if ~isempty(j)
        divisionmodes(i)=j;
    end
end

A =  abs(cellArea(tissue,1:ncells))';
D1 = (divisionmodes==1);
D2 = (divisionmodes==2);
D3 = (divisionmodes==3);

P1 = arrayfun(@(x) 1/sum(exp(-alpha*([planes(x,:).walllength]-planes(x,1).walllength)/sqrt(A(x)))),1:length(planes));
P2 = arrayfun(@(x) 1/sum(exp(-alpha*([planes(x,:).walllength]-planes(x,2).walllength)/sqrt(A(x)))),1:length(planes));

[P1 index] = sort(P1);
D1 =D1(index);
P1_Zinnia = zeros(floor(length(P1)/N),1);
dP1_Zinnia = zeros(size(P1_Zinnia));
r1_Zinnia = zeros(size(P1_Zinnia));
dr1_Zinnia = zeros(size(P1_Zinnia));
for i =1:floor(length(P1)/N)
    P1_Zinnia(i) = mean(P1((i-1)*N+1:i*N));
    dP1_Zinnia(i) = std(P1((i-1)*N+1:i*N));
    r1_Zinnia(i) = mean(D1((i-1)*N+1:i*N));
end


[P2 index] = sort(P2);
D2 =D2(index);
P2_Zinnia = zeros(floor(length(P2)/N),1);
dP2_Zinnia = zeros(size(P2_Zinnia));
r2_Zinnia = zeros(size(P2_Zinnia));
dr2_Zinnia = zeros(size(P2_Zinnia));
for i =1:floor(length(P2)/N)
    P2_Zinnia(i) = mean(P2((i-1)*N+1:i*N));
    dP2_Zinnia(i) = std(P2((i-1)*N+1:i*N));
    r2_Zinnia(i) = mean(D2((i-1)*N+1:i*N));
end

%% Load Fern results

divisionmodes = [];
P1=[];
P2=[];
for i=1:5
    load(['./Microsorum' num2str(i) '_1-treated.mat']);
    tissue = simplifiedtissue;

    ncells = length(originaltissue.c)-length(tissue.c); 
    real_edges = sort(edges_list,2);

    %cell_list = 1:ncells;
    D = NaN(1,ncells);

    for i=1:ncells
        ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)';
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows');
        if ~isempty(j)
            D(i)=j;
        end
    end
    divisionmodes=[divisionmodes D];
    A =  abs(cellArea(tissue,1:ncells))';
    P1 = [P1 arrayfun(@(x) 1/sum(exp(-alpha*([planes(x,:).walllength]-planes(x,1).walllength)/sqrt(A(x)))),1:length(planes))];
    P2 = [P2 arrayfun(@(x) 1/sum(exp(-alpha*([planes(x,:).walllength]-planes(x,2).walllength)/sqrt(A(x)))),1:length(planes))];

end
D1= (divisionmodes==1);
D2= (divisionmodes==2);

[P1 index] = sort(P1);
D1 =D1(index);
P1_Microsorum = zeros(floor(length(P1)/N),1);
dP1_Microsorum = zeros(size(P1_Microsorum));
r1_Microsorum = zeros(size(P1_Microsorum));
dr1_Microsorum = zeros(size(P1_Microsorum));
for i =1:floor(length(P1)/N)
    P1_Microsorum(i) = mean(P1((i-1)*N+1:i*N));
    dP1_Microsorum(i) = std(P1((i-1)*N+1:i*N));
    r1_Microsorum(i) = mean(D1((i-1)*N+1:i*N));
end


[P2 index] = sort(P2);
D2 =D2(index);
P2_Microsorum = zeros(floor(length(P2)/N),1);
dP2_Microsorum = zeros(size(P2_Microsorum));
r2_Microsorum = zeros(size(P2_Microsorum));
dr2_Microsorum = zeros(size(P2_Microsorum));
for i =1:floor(length(P2)/N)
    P2_Microsorum(i) = mean(P2((i-1)*N+1:i*N));
    dP2_Microsorum(i) = std(P2((i-1)*N+1:i*N));
    r2_Microsorum(i) = mean(D2((i-1)*N+1:i*N));
end


axes('OuterPosition',[1/2 0 1/2 1]);
hold on
x0 = .17;
a = .17;
b = .1;
phi = atan(-b/a);
theta1 = acos(x0/(sqrt(a^2+b^2))) -phi-pi;
theta2 = pi-acos(x0/(sqrt(a^2+b^2))) +phi;

theta =[theta1:pi/100:theta2 theta2];
X = x0+a*cos(theta)+b*sin(theta);
Y = x0+a*cos(theta)-b*sin(theta);
fill([X 0 X(1)],[Y 0 Y(1)],[.5 .5 .5],'EdgeColor','none','FaceAlpha',.5)

x0 = .7;
a = .3;
b = .15;
phi = atan(-b/a);
theta1 = acos((1-x0)/(sqrt(a^2+b^2))) -phi-pi;
theta2 = pi-acos((1-x0)/(sqrt(a^2+b^2))) +phi;

theta =[theta2-3*pi:pi/100:theta1+pi theta1+pi];
X = x0+a*cos(theta)+b*sin(theta);
Y = x0+a*cos(theta)-b*sin(theta);
fill([X 1 X(1)],[Y 1 Y(1)],[.5 .5 .5],'EdgeColor','none','FaceAlpha',.5)
msize =10;
errorbar2(P1_Coleochaete,dP1_Coleochaete,r1_Coleochaete,dr1_Coleochaete,...
    'sr','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize);
errorbar2(P2_Coleochaete,dP2_Coleochaete,r2_Coleochaete,dr2_Coleochaete,...
    'sr','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize);
errorbar2(P1_Dionaea,dP1_Dionaea,r1_Dionaea,dr1_Dionaea,...
    '^b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize);
errorbar2(P2_Dionaea,dP2_Dionaea,r2_Dionaea,dr2_Dionaea,...
    '^b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize);
errorbar2(P1_Zinnia,dP1_Zinnia,r1_Zinnia,dr1_Zinnia,...
    'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
errorbar2(P2_Zinnia,dP2_Zinnia,r2_Zinnia,dr2_Zinnia,...
    'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
errorbar2(P1_Microsorum,dP1_Microsorum,r1_Microsorum,dr1_Microsorum,...
    'dg','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize);
errorbar2(P2_Microsorum,dP2_Microsorum,r2_Microsorum,dr2_Microsorum,...
    'dg','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize);
plot([0 1],[0 1],'--k');
axis([0 1 0 1]);

set(gca,'FontName','Helvetica','FontSize',16)
xlabel('P_i','FontSize',20,'FontName','Helvetica');
ylabel('n_i/(n_1+...+n_N)','FontSize',20,'FontName','Helvetica');

set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
box on

% Labels
axes('Position',[0 1-labelsize labelsize labelsize]);
text(0,0,'a','FontWeight','bold','FontName','Helvetica');
axis equal off
axes('Position',[1/2 1-labelsize labelsize labelsize]);
text(0,0,'b','FontWeight','bold','FontName','Helvetica');
axis equal off
%plot2svg_beta('Figure4.svg',gcf);
