  tic;
NACC_DataFolder = pwd;
% filename0 = [NACC_DataFolder,'investigator_ftldlbd_nacc63.csv'];
% filename = [NACC_DataFolder,'/investigator_nacc64.csv'];
filename = [NACC_DataFolder,'/investigator_nacc71.csv'];
% filename = [NACC_DataFolder,'/investigator_nacc67.csv'];
% filenameMRI = [NACC_DataFolder,'investigator_scan_mriqc_nacc67.csv'];
% filenameCSF = [NACC_DataFolder,'investigator_fcsf_nacc67.csv'];
% -----------   Read in ADNI csv data files into matlab data structures -------------
T1 = readtable(filename,'ReadVariableNames',true); toc
%%
tic;
% T1MRI = readtable(filenameMRI,'ReadVariableNames',true); 
% T1CSF = readtable(filenameCSF,'ReadVariableNames',true); toc
% x2.csf = readItemLevelCSF(T1CSF); 
%% x2.im = readItemLevelImaging(T1)
tic 
x2.mmse1 = readItemLevelMMSE(T1); toc % detailed item-level MMSE data
%%
x2.moca1 = readItemLevelMOCA(T1); toc % detailed item-level MOCA data
%%
tic
x2.faq1 = readItemLevelFAQ(T1);
x2.cb1a = readItemLevelCogBat1a(T1);
x2.cb1b = readItemLevelCogBat1b(T1);
x2.cb2a = readItemLevelCogBat2a(T1);
x2.cb2b = readItemLevelCogBat2b(T1);
toc
%%
q = x2.moca1;
for k=1:length(q.idList),
    i = find(q.id==q.idList(k));
    q.ii{k} = i;
    q.iDate{k} = q.regDate(i);
    q.iDays{k} = [0; days(diff(q.iDate{k}))];
    q.iDx{k} = q.dx(i);
    q.iNtot(k) = length(q.iDx{k});
    q.iNCN(k) = sum(q.iDx{k}==0);
    q.iNMCI(k) = sum(q.iDx{k}==1);
    q.iNAD(k) = sum(q.iDx{k}==2);
    % if q.iNtot(k) > 1,
    q.iNDxChanges(k) = sum(diff(q.iDx{k})~=0);
    q.iADonset(k) = max([0,min(find(q.iDx{k}==2))]);
    q.iMCIonset(k) = max([0,min(find(q.iDx{k}==1))]);
end
    
%%
x2.im = readItemLevelImaging(T1);
dm = x2.moca1.datSum0;
dx0 = x2.moca1.dx0;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
for k=1:6,
    d = x2.im.dat0(:,k);
    i0 = (dx0 == 0 & isfinite(d));
    i1 = (dx0 == 1 & isfinite(d));
    i2 = (dx0 == 2 & isfinite(d));
    i0m = (dx0 == 0 & isfinite(d) & isfinite(dm));
    i1m = (dx0 == 1 & isfinite(d) & isfinite(dm));
    i2m = (dx0 == 2 & isfinite(d) & isfinite(dm));
    x2.im.n01(k,:) = [sum(i0), sum(i1), sum(i0+i1)]; 
    x2.im.n12(k,:) = [sum(i1), sum(i2), sum(i1+i2)]; 
    x2.im.n02(k,:) = [sum(i0), sum(i2), sum(i0+i2)]; 
    
    x2.im.auc01moca(k) = 1-roc1(dm(i0m|i1m),dx0(i0m|i1m)); 
    x2.im.auc12moca(k) = 1-roc1(dm(i1m|i2m),dx0(i1m|i2m)-1); 
    x2.im.auc02moca(k) = 1-roc1(dm(i0m|i2m),dx0(i0m|i2m)/2); 
    x2.im.auc01(k) = roc1(d(i0|i1),dx0(i0|i1)); 
    x2.im.auc12(k) = roc1(d(i1|i2),dx0(i1|i2)-1); 
    x2.im.auc02(k) = roc1(d(i0|i2),dx0(i0|i2)/2); 
end
%%
for k1=1:6, for k2=1:6,
    d = -x2.im.dat0(:,[k1,k2]);
    n2(k1,k2) = sum(isfinite(sum(d,2)));
    n2moca(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0));
    n2mocaa(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0+dx0a));
    n2moca0(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==0);
    n2moca1(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==1);
    n2moca2(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==2);
    n3mmse(k1,k2) = sum(isfinite(sum(d,2)+x2.mmse1.datSum0));
    i0 = (dx0 == 0 & isfinite(sum(d,2)));
    i1 = (dx0 == 1 & isfinite(sum(d,2)));
    i2 = (dx0 == 2 & isfinite(sum(d,2)));
    i0m = (dx0 == 0 & isfinite(sum(d,2)) & isfinite(dm));
    i1m = (dx0 == 1 & isfinite(sum(d,2)) & isfinite(dm));
    i2m = (dx0 == 2 & isfinite(sum(d,2)) & isfinite(dm));
    zz = rocMaximizeLineSearch(d(i0|i1,:),dx0(i0|i1)); x2.im.auc01b(k1,k2) = zz.auc;
    zz = rocMaximizeLineSearch(d(i1|i2,:),dx0(i1|i2)-1); x2.im.auc12b(k1,k2) = zz.auc; 
    zz = rocMaximizeLineSearch(d(i0|i2,:),dx0(i0|i2)/2); x2.im.auc02b(k1,k2) = zz.auc;
    x2.im.auc01bmoca(k1,k2) = 1-roc1(dm(i0m|i1m),dx0(i0m|i1m)); 
    x2.im.auc12bmoca(k1,k2) = 1-roc1(dm(i1m|i2m),dx0(i1m|i2m)-1); 
    x2.im.auc02bmoca(k1,k2) = 1-roc1(dm(i0m|i2m),dx0(i0m|i2m)/2); 
end; end

for k1=1:6, for k2=1:6,
    d = -x2.im.dat0(:,[k1,k2,4]);
    n3(k1,k2) = sum(isfinite(sum(d,2)));
    n3moca(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0));
    n3mocaa(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0+dx0a));
    n3moca0(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==0);
    n3moca1(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==1);
    n3moca2(k1,k2) = sum(isfinite(sum(d,2)+x2.moca1.datSum0) & dx0a==2);
    n3mmse(k1,k2) = sum(isfinite(sum(d,2)+x2.mmse1.datSum0));
    i0 = (dx0 == 0 & isfinite(sum(d,2)));
    i1 = (dx0 == 1 & isfinite(sum(d,2)));
    i2 = (dx0 == 2 & isfinite(sum(d,2)));
    i0m = (dx0 == 0 & isfinite(sum(d,2)) & isfinite(dm));
    i1m = (dx0 == 1 & isfinite(sum(d,2)) & isfinite(dm));
    i2m = (dx0 == 2 & isfinite(sum(d,2)) & isfinite(dm));
    zz = rocMaximizeLineSearch(d(i0|i1,:),dx0(i0|i1)); x2.im.auc01c(k1,k2) = zz.auc;
    zz = rocMaximizeLineSearch(d(i1|i2,:),dx0(i1|i2)-1); x2.im.auc12c(k1,k2) = zz.auc; 
    zz = rocMaximizeLineSearch(d(i0|i2,:),dx0(i0|i2)/2); x2.im.auc02c(k1,k2) = zz.auc;
    x2.im.auc01cmoca(k1,k2) = 1-roc1(dm(i0m|i1m),dx0(i0m|i1m)); 
    x2.im.auc12cmoca(k1,k2) = 1-roc1(dm(i1m|i2m),dx0(i1m|i2m)-1); 
    x2.im.auc02cmoca(k1,k2) = 1-roc1(dm(i0m|i2m),dx0(i0m|i2m)/2); 

end; end
%%
dx0 = x2.moca1.dx0;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;

i1=setdiff([1:13],2);
cbForMMSE.dat = x2.cb1b.dat0(:,i1); cbForMMSE.dat(:,7:8)=30-cbForMMSE.dat(:,7:8)/10;
cbForMMSE.dat0 = cbForMMSE.dat;
cbForMMSE.vars = x2.cb1b.vars(i1);
cbForMMSE.dx0 = dx0a; 
cbForMMSE.dx0a = dx0a; 

i1=[3,6:9,12:13];
i2=[1,3,4,7];
cbForMOCA.dat = [x2.cb1b.dat0(:,i1),x2.cb2b.dat0(:,i2)]; cbForMOCA.dat(:,4:5)=30-cbForMOCA.dat(:,4:5)/10;
cbForMOCA.dat0 = cbForMOCA.dat;
cbForMOCA.vars = {x2.cb1b.vars{i1},x2.cb2b.vars{i2}};
cbForMOCA.dx0 = dx0a; 
cbForMOCA.dx0a = dx0a; 
%%
dx0 = x2.moca1.dx0;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
x2.cb = organizeCogBatData(x2);
%%
x2.buildBase.dat0{1} = [x2.moca1.datMIS0all(:,1), 5*x2.moca1.dat0(:,[25,28,29])]; x2.buildBase.name{1}='Reweighted MoCA-8_R_O';
x2.buildBase.dat0{2} = [x2.moca1.datSum0,0*x2.moca1.datSum0]; x2.buildBase.name{2}='Standard MoCA-30';
x2.buildBase.dat0{2} = [x2.moca1.datSum0]; x2.buildBase.name{2}='Standard MoCA-30';
i2=21; 
x2.buildBase.dat0{3} = x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]).*[1,4*ones(1,29)]; x2.buildBase.name{3}='Reweighted MoCA-30';
x2.buildBase.dat0{4} = x2.moca1.datMIS0all(:,1); x2.buildBase.name{4}='MoCA-5';
x2.buildBase.dat0{5} = x2.moca1.dat0(:,20:24); x2.buildBase.name{5}='Reweighted MoCA-5';


% x2.buildBase.dat0{4} = [x2.mmse1.datSum0,0*x2.mmse1.datSum0]; x2.buildBase.name{4}='Standard MMSE-30';
% x2.buildBase.dat0{4} = [x2.mmse1.datSum0]; x2.buildBase.name{4}='Standard MMSE-30';
% x2.buildBase.dat0{4}(x2.buildBase.dat0{4}>30,:)=NaN;

%% Make paper 1 fig 4 v1 
% make Thresholded + granular single-item AUC improvement plots
tic;   
if iscell(x2.cb.MOCA1.dat0), x2.cb.MOCA1.dat0 = cell2mat(x2.cb.MOCA1.dat0); end
if iscell(x2.cb.MMSE1.dat0), x2.cb.MMSE1.dat0 = cell2mat(x2.cb.MMSE1.dat0); end
figure; set(gcf,'Position',[0 0 1160 800]);
subplot(2,5,[1 2]); plotThresholdedAUCs(x2.cb.MOCA1,[0 10 20 40 60 80 90 100], x2.moca1, 'MOCA'); toc
subplot(2,5,[6 7]); plotThresholdedAUCs(x2.cb.MMSE1,[0 10 20 40 60 80 90 100], x2.mmse1, 'MMSE'); 
subplot(2,5,[6 7]); h = title('Can any coginitive battery items supplement the MMMSE?'); h.Position = h.Position + [0, 0.001, 0]; toc
subplot(2,5,[1 2]); h = title('Can any coginitive battery items supplement the MoCA?'); h.Position = h.Position + [0, 0.001, 0];

% make granular-only multi-item-combination AUC improvement plots
GREY = 0.7*[1 1 1]; colours=colororder('glow12'); symbol1 = {'^','square','diamond'};
% make plots for 3 MOCA bases
qq = x2.cb.MOCA2s1; clear xBuild h h0; 
N=length(qq.vars); 
subplot(2,5,[4 10]); hold on; 
for k1=1:3,
    if isfield(qq,'c0'), for k=1:N, qq.dat0{k} = qq.dat0{k}.*qq.c0{k}; end; end % apply initial scaling coeffs.
    for k=0:N, q=[x2.buildBase.dat0{k1}, horzcat(qq.dat0{1:k})]; if k>1, w0=z.w; else w0=[]; end;  z=rocMaximizeLineSearchA(q,dx0a,w0); xBuild.Aroc(k+1)=z.auc; zz{k1,k+1}=z; end; disp(['MOCA combo', num2str(k1)]), z, toc
    for k=1:N, h(k) = plot(N-k+[1:2], xBuild.Aroc(N-k+[1:2]), '.-', 'Color', colours(qq.colorOrder(N-k+1),:), 'MarkerSize', 30, 'linewidth',1);  hold on; end
    h0(k1) = plot(1, xBuild.Aroc(1), '.', 'Marker', symbol1{k1},  'MarkerSize', 12, 'MarkerFaceColor', GREY, 'MarkerEdgeColor', GREY);
    % plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid; hold on;
end
% ylim([0.8 0.90]); xlim([0.5,N+1+0.5]); grid on;
% legend([h,h0],[qq.names, x2.buildBase.name(1:3)])
mocaNames = qq.names; mocaH=fliplr(h); mocaH0=fliplr(h0);   


% make plots for 1 MMSE bases (other bases not possible for NACC data that doesn't itemize MMSE well)
qq = x2.cb.MMSE2a; clear xBuild h h0; 
N=length(qq.vars); 
if isfield(qq,'c0'), for k=1:N, qq.dat0{k} = qq.dat0{k}.*qq.c0{k}; end; end
for k=0:N, q=[x2.buildBase.dat0{4}, horzcat(qq.dat0{1:k})]; if k>1, w0=z.w; else w0=[]; end;  z=rocMaximizeLineSearchA(q,dx0a,w0); xBuild.Aroc(k+1)=z.auc; end; z, toc
for k=1:N, h(k) = plot(N-k+[1:2], xBuild.Aroc(N-k+[1:2]), '.-', 'Color', colours(qq.colorOrder(N-k+1),:), 'MarkerSize', 30, 'linewidth',1);  hold on; end
h0(1) = plot(1, xBuild.Aroc(1), '.', 'Marker', 'o',  'MarkerSize', 12, 'MarkerFaceColor', GREY, 'MarkerEdgeColor', GREY);
ylim([0.75 0.90]); xlim([0.5,N+1+0.5]); grid on;
mmseNames = qq.names; mmseH=fliplr(h); mmseH0=fliplr(h0);   

legend([mocaH0, mmseH0, mocaH, mmseH(4:end)],[x2.buildBase.name(1:4), mocaNames, mmseNames(4:end)])
%legend([h,h0(1:4)],[fliplr(qq.names), x2.buildBase.name(4)])
%%
z1=rocMaximizeLineSearchA(q(:,1:3),dx0a,0.01*ones(1,2));
%%
z2=rocMaximizeLineSearchA(q(:,1:4),dx0a,0.1*ones(1,2))
z3=rocMaximizeLineSearchA(q(:,1:5),dx0a,0.1*ones(1,2))
z1a=rocMaximizeLineSearchA(q(:,1:3),dx0a,1*ones(1,2))
z2a=rocMaximizeLineSrearchA(q(:,1:4),dx0a,1*ones(1,2))
z3a=rocMaximizeLineSearchA(q(:,1:5),dx0a,1*ones(1,2))

%% Make paper 1 fig 4 v2  
% PANEL (A): make Thresholded + granular single-item AUC improvement plots
tic;   
if iscell(x2.cb.MOCA1.dat0), x2.cb.MOCA1.dat0 = cell2mat(x2.cb.MOCA1.dat0); end
if iscell(x2.cb.MMSE1.dat0), x2.cb.MMSE1.dat0 = cell2mat(x2.cb.MMSE1.dat0); end
figure; set(gcf,'Position',[0 0 1160 800]);
subplot(2,5,[1 2]); plotThresholdedAUCs(x2.cb.MOCA1,[0 10 20 40 60 80 90 100], x2.moca1, 'MOCA'); toc
subplot(2,5,[1 2]); h = title('Can any coginitive battery items help to improve on the MoCA?'); h.Position = h.Position + [0, 0.001, 0];
yt=get(gca,'yticklabel'); for k=2:2:length(yt), yt{k}=[]; yticklabels(yt); end;
%%
figure
%%

% PANEL (C): make granular-only multi-item-combination AUC improvement plots
GREY = 0.7*[1 1 1]; colours=colororder('glow12'); symbol1 = {'^','square','diamond','pentagram','.'};
% make plots for 3 MOCA bases
qq = x2.cb.MOCA2s1; clear xBuild h h0; 
N=length(qq.vars); 
subplot(2,5,[4 10]); hold on; clear zz;
for k1=1:length(x2.buildBase.dat0),   % was "3"

    if isfield(qq,'c0'), for k=1:N, qq.dat0{k} = qq.dat0{k}.*qq.c0{k}; end; end % apply initial scaling coeffs.
    
    %k1=1; qq.dat0{1}=qq.dat0{1}(:,1);
    for k=0:N, q=[x2.buildBase.dat0{k1}, horzcat(qq.dat0{1:k})]; if k>1, w0=z.w; else w0=[]; end;  z=rocMaximizeLineSearchA(q,dx0a,w0,0.2); xBuild.Aroc(k+1)=z.auc; zz{k1,k+1}=z; end; disp(['MOCA combo', num2str(k1)]), z, toc
    for k=0:N, q=[x2.buildBase.dat0{k1}, horzcat(qq.dat0{1:k})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); xBuild.Aroc(k+1)=max(z.auc,xBuild.Aroc(k+1)); zz2{k1,k+1}=z; end; disp(['w0=0 MOCA combo', num2str(k1)]), z, toc
    % for k=N-1:-1:0, q=[x2.buildBase.dat0{k1}, horzcat(qq.dat0{1:k})]; w0=zz{k1,k+1}.w; w0=w0(1:size(q,2)-1); z=rocMaximizeLineSearchA(q,dx0a,w0); xBuild.Aroc(k+1)=max(z.auc,xBuild.Aroc(k+1)); zz2{k1,k+1}=z; end; disp(['reverse MOCA combo', num2str(k1)]), z, toc
    if k1==1, LW=2.5; else LW=1; end
    for k=1:N, h(k) = plot(N-k+[1:2], xBuild.Aroc(N-k+[1:2]), '.-', 'Color', colours(qq.colorOrder(N-k+1),:), 'MarkerSize', 30, 'linewidth',LW);  hold on; end
    
    h0(k1) = plot(1, xBuild.Aroc(1), '.', 'Marker', symbol1{k1},  'MarkerSize', 12, 'MarkerFaceColor', GREY, 'MarkerEdgeColor', GREY);
    drawnow; k1, zz, zz2% pause
    % plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid; hold on;
end
mocaNames = qq.names; mocaH=h; mocaH0=fliplr(h0);   
legend([mocaH, mocaH0],[fliplr(mocaNames), fliplr(x2.buildBase.name(1:3)), ])
xlabel('Number of added items'); ylabel('Area under the ROC curve'); shg 
xlim([0.5 5.5]); xticks([1:5]); xticklabels([0:4]);

% draw inset-axis ROC curves for MoCA-8 case:
haxi = axes('Position',[.7,.2,.2,.2]); plot([0 1],[0 1],':k','LineWidth',1); axis square; hold on
clear r rq
for k=0:N, 
    q = [x2.buildBase.dat0{1}, horzcat(qq.dat0{1:k})]; 
    ii = isfinite(sum(q,2)) & (dx0a==0 | dx0a==1); 
    [fpr, tpr, r(k+1)] = roc1(-q(ii,:)*zz{1,k+1}.w1',dx0a(ii)); 
    q = q./std(q(ii,:));
    [fpr, tpr, rq(k+1)] = roc1(-q(ii,:)*[1,zz{1,k+1}.w]',dx0a(ii)); 
    if k==0, colour1 = GREY; else colour1 = colours(qq.colorOrder(k),:); end
    plot(fpr,tpr, 'Color', colour1, 'linewidth', 2); hold on;
end;
xlabel('1-Specificity'); ylabel('Sensitivity'); title('Full ROC curves')
%% end Fig 4 v2 from paper 1
tic
qqm=x2.moca1; a = testMatchingDifferentDate (qqm,qqm,[], [],'less than'); i1=a.i1; i2=a.i2; toc
%%
tic;
qqm=x2.moca1; aa = testMatchingDifferentDate (qqm,qqm,[], [],'greater than'); i1=aa.i1; i2=aa.i2; toc
ifi = find(qqm.i);
mocaMem5_ = qqm.da;
mocaMem5_20aa = NaN*qqm.i; mocaMem5_20aa(ifi(aa.i1)) = mocaMem5_(aa.i2);
%%
tic;
qqm=x2.moca1; aa = testMatchingDifferentDate (qqm,qqm,[], [],'greater than'); i1=aa.i1; i2=aa.i2; toc
mocaMem5_ = qqm.datMIS0all(qqm.i,1)/3; 
mocaMem5_20aa = NaN*qqm.i; mocaMem5_20aa(ifi(aa.i1)) = mocaMem5_(aa.i2);

%%
mocaMem5_ = qqm.datMIS0all(qqm.i,1)/3; 
mocaMem5_1 = NaN*mocaMem5_; mocaMem5_1(a.i1) = mocaMem5_(a.i1);
mocaMem5_2 = NaN*mocaMem5_; mocaMem5_2(a.i2) = mocaMem5_(a.i2);
ifi = find(qqm.i);
mocaMem5_10 = NaN*qqm.i; mocaMem5_10(ifi(a.i1)) = mocaMem5_(a.i1);
mocaMem5_10a = NaN*qqm.i; mocaMem5_10a(ifi) = mocaMem5_1;
mocaMem5_20 = NaN*qqm.i; mocaMem5_20(ifi(a.i1)) = mocaMem5_(a.i2);
%%
q1=[x2.buildBase.dat0{1}+0*mocaMem5_10+0*sum(horzcat(x2.cb.MOCA2s.dat0{:}),2)]; disp('Extra MoCA recall info!!!')
% q1=[x2.buildBase.dat0{1}]; disp('Extra MoCA recall info!!!')
q=[q1]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, -0.01*mocaMem5_10]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_10]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_20]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_20aa]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_10, mocaMem5_20]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_10, mocaMem5_20aa]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_10, mocaMem5_20, mocaMem5_20aa]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, mocaMem5_20, mocaMem5_20aa]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])

%%
z=rocMaximizeLineSearchA(mocaMem5_10,dx0a,[],0.2); disp([z.auc, z.aucPre])
z=rocMaximizeLineSearchA(mocaMem5_20,dx0a,[],0.2); disp([z.auc, z.aucPre])
z=rocMaximizeLineSearchA(mocaMem5_20aa,dx0a,[],0.2); disp([z.auc, z.aucPre])

z=rocMaximizeLineSearchA([mocaMem5_20, mocaMem5_10],dx0a,[],0.2); disp([z.auc, z.aucPre])
z=rocMaximizeLineSearchA([mocaMem5_20, mocaMem5_20aa],dx0a,[],0.2); disp([z.auc, z.aucPre])
z=rocMaximizeLineSearchA([mocaMem5_10, mocaMem5_20aa],dx0a,[],0.2); disp([z.auc, z.aucPre])
z=rocMaximizeLineSearchA([mocaMem5_10, mocaMem5_20, mocaMem5_20aa],dx0a,[],0.2); disp([z.auc, z.aucPre])

%%
q=[q1, horzcat(x2.cb.MOCA2s.dat0{:})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); disp([z.auc, z.aucPre])
q=[q1, horzcat(x2.cb.MOCA2s.dat0{:}), mocaMem5_20]; z=rocMaximizeLineSearchA(q,dx0a,[],0.02); disp([z.auc, z.aucPre])
q=[q1, horzcat(mocaMem5_20, x2.cb.MOCA2s.dat0{:})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.02); disp([z.auc, z.aucPre])
%[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))-qq.recall(a.i1(i1))/2.56+30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
%%
q1=[x2.buildBase.dat0{1}+0*mocaMem5_10];                                   z=rocMaximizeLineSearchA(q1,dx0a,[],0.02); disp([z.auc, z.aucPre])
q1=[x2.buildBase.dat0{1}]; q1(:,1)=(q(:,1)/3+mocaMem5_10)/2; z=rocMaximizeLineSearchA(q1,dx0a,[],0.02); disp([z.auc, z.aucPre])
q1=[x2.buildBase.dat0{1}]; q1(:,1)=(q(:,1)/3+mocaMem5_20)/2; z=rocMaximizeLineSearchA(q1,dx0a,[],0.02); disp([z.auc, z.aucPre])

%%





% a = x2.cb.MOCA1.dat0(dx0a==1,[5]); a=a(isfinite(a)); mean(a>=11),
% b = x2.cb.MOCA1.dat0(dx0a==1,[6]); b=b(isfinite(b)); mean(b>=11), prctile(b,100-mean(a>=11)*100)
% xh = [0:40];
% a = x2.cb.MOCA1.dat0(dx0a==1,[5]); a=a(isfinite(a)); [ya,xa]=hist

% PANEL (B): detailed specific example comparison of LF vs CF, also showing
% the effect of optimizing the threshold/sensitivity and the additional
% benefit of using the full granular data from a multi-valued test-item

qmoca1 = [x2.moca1.datSum0] - x2.moca1.dat0(:,17); % MoCA-29 data
for k=1:20, 
    z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5)>=k],dx0a); qqF(k)= z.auc;
    z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,6)>=k],dx0a); qqL(k)= z.auc;
    z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,1)>=k],dx0a); qqA(k)= z.auc;
    z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,2)>=k],dx0a); qqV(k)= z.auc;
end
[~,qqFimax] = max(qqF);
[~,qqLimax] = max(qqL);
[~,qqAimax] = max(qqA);
[~,qqVimax] = max(qqV);
clear qqq;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5)>=0],dx0a);       qqq(1)=z.auc;  % reference level for MoCA-29 AUC
z=rocMaximizeLineSearchA([qmoca1,x2.moca1.dat0(:,17)],dx0a);            qqq(2)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5)>=11],dx0a);      qqq(3)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,6)>=10],dx0a);      qqq(4)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5)>=qqFimax],dx0a); qqq(5)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,6)>=qqLimax],dx0a); qqq(6)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5)],dx0a);          qqq(7)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,6)],dx0a);          qqq(8)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,5:6)],dx0a);        qqq(9)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,1)>=qqFimax],dx0a); qqq(10)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,2)>=qqLimax],dx0a); qqq(11)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,1)],dx0a);          qqq(12)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,2)],dx0a);          qqq(13)=z.auc;
z=rocMaximizeLineSearchA([qmoca1,x2.cb.MOCA1.dat0(:,1:2)],dx0a);        qqq(14)=z.auc;

subplot(2,5,[6 7]);  
xBar=[1:3, 5:6, 8:10, 12:13 15:17]; h=bar(xBar, qqq(2:end)-qqq(1));
xticks([]);
colours=colororder('glow12'); c = colours(x2.cb.MOCA1.colorOrder([5,5,6,5,6,5,6,6,1,2,1,2,2]),:); 
c(8,:) = 0.5*mean(c(6:7,:))+0.5*[1 1 1];  % make the color for combo-bar a lightened version of the mean color
c(13,:) = 0.5*mean(c(11:12,:))+0.5*[1 1 1]; % make the color for combo-bar a lightened version of the mean color
h.FaceColor='flat'; h.CData=c; shg

ylabel ({'AUC Improvement from a 29-item version of the MoCA','(that excudes the Letter Fluency item)'})
title('Example: Comparison of Category Fluency & Letter Fluency')

% make the fairly-complicated x-axis labeling
qyspan = diff(gca().YLim); qdx=0.8; 
qdy = 0.02;
hl(1)=line([1*[1 1],3*[1 1]]   +qdx*[-1,-1,1,1],-qdy*qyspan*[0.5 1 1 0.5]); ht(1)=text(2,-qdy*qyspan,{'w/ MoCA','threshold'});
hl(2)=line([5*[1 1],6*[1 1]]   +qdx*[-1,-1,1,1],-qdy*qyspan*[0.5 1 1 0.5]); ht(2)=text(5.5,-qdy*qyspan,{'w/ optimal','threshold'});
hl(3)=line([8*[1 1],10*[1 1]]  +qdx*[-1,-1,1,1],-qdy*qyspan*[0.5 1 1 0.5]); ht(3)=text(9,-qdy*qyspan,{'for full','granular data'});
hl(4)=line([12*[1 1],13*[1 1]] +qdx*[-1,-1,1,1],-qdy*qyspan*[0.5 1 1 0.5]); ht(4)=text(12.5,-qdy*qyspan,{'w/ optimal','threshold'});
hl(5)=line([15*[1 1],17*[1 1]] +qdx*[-1,-1,1,1],-qdy*qyspan*[0.5 1 1 0.5]); ht(5)=text(16,-qdy*qyspan,{'for full','granular data'});
qdy = 0.15;
hl(6)=line([1*[1 1],10*[1 1]]  +qdx*[-1,-1,1,1],-qdy*qyspan*[0.9 1 1 0.9]); ht(6)=text(5.5,-qdy*qyspan,{'60-sec Letter Fluency','(included in MoCA)'});
hl(7)=line([12*[1 1],17*[1 1]] +qdx*[-1,-1,1,1],-qdy*qyspan*[0.9 1 1 0.9]); ht(7)=text(14.5,-qdy*qyspan,{'60-sec Category Fluency','(not in MoCA)'});

horizList = ['FFLFLFL']
for k=1:7, htBar(k)=text(xBar(k),0,horizList(k)); end
for k=1:7, set(htBar(k),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',[1 1 1]); end
vertList = {' F+L',' Animal',' Veg',' Animal',' Veg',' Animal + Veg'}
for k=8:13, htBar(k)=text(xBar(k),0,vertList(k-7),'Rotation',90); end
for k=8:13, set(htBar(k),'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',11,'FontWeight','bold','Color',[1 1 1]); end

for k=1:length(hl), set(hl(k),'Color','k','LineWidth',1.5); end
for k=1:length(ht), set(ht(k),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',11); end

ylim([0 qyspan]); set(gca,'Clipping','off');
yt=get(gca,'yticklabel'); for k=2:2:length(yt), yt{k}=[]; yticklabels(yt); end;
% exportgraphics(gcf, 'paper1_fig4b.pdf', 'ContentType', 'vector');
grid
hai = axes('Position',[.175 .195 .11 .2]);
h=bar(xBar(1:8), qqq(2:9)-qqq(1)); 
ylim([0 0.0048]); grid
xticks([]); yticklabels(num2str(get(gca,'ytick')'));
h.FaceColor='flat'; h.CData=c(1:8,:); %xlim([.15 .37]); ylim([.012 .015]); 
title('Zoomed in LF results')

%% More detailed build-up analysis: Taking a look at ordering
figure; 
%%
tic;
clear xBuild h h0; 
q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2s;  
% k=[1 1 2 4]; qq.dat0=qq.dat0(k); qq.names=qq.names(k);  qq.colorOrder=qq.colorOrder(k);
k=[1 2 3 4]; qq.dat0=qq.dat0(k); qq.names=qq.names(k);  qq.colorOrder=qq.colorOrder(k);
% qq = x2.cb.MOCA2s; 
% qq.dat0=[{mocaMem5_20}, qq.dat0, {-x2.faq1.dat0}]; 
% qq.names=[{'MOCA-5'}, qq.names, {'FAQ'}]; 
% qq.colorOrder=[6, qq.colorOrder, 10]; 
% q1 = q1+0*sum(horzcat(qq.dat0{:}),2);
% k=[1:6]; qq.dat0=qq.dat0(k); qq.names=qq.names(k);  qq.colorOrder=qq.colorOrder(k);
N=length(qq.names); 
subplot(1,1,1); hold on; 

% if isfield(qq,'c0'), for k=1:N, qq.dat0{k} = qq.dat0{k}.*qq.c0{k}; end; end % apply initial scaling coeffs.
% qq.dat0{1}=qq.dat0{1}(:,1);
for k=0:N, q=[q1, horzcat(qq.dat0{1:k})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.02); xBuild.Aroc(k+1)=z.auc; zz{k1,k+1}=z; end; disp(['MOCA combo', num2str(k1)]), z, toc
for k=0:N, q=[q1, horzcat(qq.dat0{1:k})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); xBuild.Aroc(k+1)=max(z.auc,xBuild.Aroc(k+1)); zz2{k1,k+1}=z; end; disp(['w0=0 MOCA combo', num2str(k1)]), z, toc
    % for k=N-1:-1:0, q=[x2.buildBase.dat0{k1}, horzcat(qq.dat0{1:k})]; w0=zz{k1,k+1}.w; w0=w0(1:size(q,2)-1); z=rocMaximizeLineSearchA(q,dx0a,w0); xBuild.Aroc(k+1)=max(z.auc,xBuild.Aroc(k+1)); zz2{k1,k+1}=z; end; disp(['reverse MOCA combo', num2str(k1)]), z, toc
% if k1==1, LW=2.5; else LW=1; end
for k=1:N, h(k) = plot(N-k+[1:2], xBuild.Aroc(N-k+[1:2]), '.-', 'Color', colours(qq.colorOrder(N-k+1),:), 'MarkerSize', 30, 'linewidth',1);  hold on; end
% h0(k1) = plot(1, xBuild.Aroc(1), '.', 'Marker', symbol1{k1},  'MarkerSize', 12, 'MarkerFaceColor', GREY, 'MarkerEdgeColor', GREY);
% plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid; hold on;
   
legend(h,fliplr(qq.names))
xlabel('Number of added items'); ylabel('Area under the ROC curve'); shg 
xlim([0.5 N+1.5]); xticks([0:N]+1); xticklabels([0:N]);

%% continue above analysis to systematcally find final contribution for each item!
clear xBuild
q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2x;  
N=length(qq.names); 
q=[q1, horzcat(qq.dat0{:})]; zTOT=rocMaximizeLineSearchA(q,dx0a,w0,0.2); xBuild.Aroc(1)=zTOT.auc

for k=0:N, q2=[{q1},qq.dat0]; q=[horzcat(q2{setdiff(1:N+1,k+1)})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); xBuild.Aroc(k+2)=z.auc; zz1{k+1}=z; end; disp(['w0=0.2 MOCA combo', num2str(k1)]), z, toc
for k=0:N, q2=[{q1},qq.dat0]; q=[horzcat(q2{setdiff(1:N+1,k+1)})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.02); xBuild.Aroc(k+2)=max(z.auc,xBuild.Aroc(k+2)); zz1{k+1}=z; end; disp(['w0=0.02 MOCA combo', num2str(k1)]), z, toc
figure; 
subplot(1,2,1); hold on; grid
colours1 = [GREY; colours(qq.colorOrder,:)]; 
xBuild.Aroc1 = xBuild.Aroc(1) - xBuild.Aroc;
qx = xBuild.Aroc1(2:end); [~,is] = sort(qx,'descend');
h = bar(1:N+1,diag(qx(is)),'stacked');
for k=1:length(h), h(k).FaceColor=colours1(is(k),:); end
qq.names1 = [{'MoCA-8'},qq.names]; legend(h,qq.names1(is))
xlabel('Candidate Items'); ylabel('Area under the ROC curve lost by removal'); shg 
xticks([]); yticks([0:.002:0.018]); yticklabels(num2str([0:.002:0.018]'));
xlim ([0.3 N+1+0.7]); ylim([-.001, .0195]); 
%%
figure;
%%
clear xBuild
q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2s1;
% qq.dat0=[{mocaMem5_20}, qq.dat0, {-x2.faq1.dat0}]; 
% qq.names=[{'MOCA-5 (Recall)'}, qq.names, {'FAQ (Partner)'}]; 
% qq.colorOrder=[6, qq.colorOrder, 10]; 
clear faqNames; 
for k=1:10, faqNames{k} = ['F',num2str(k)]; end 
for k=1:10, faqData{k} = -x2.faq1.dat0(:,k); end 
%qq.dat0=[qq.dat0([1,2,3]), faqData]; qq.dat0
%qq.dat0=[qq.dat0([1,2,4]), {-x2.faq1.dat0}]; qq.dat0
qq.dat0=[qq.dat0(1:3)]; 
qq.names=[qq.names(1:3), faqNames]; 
qq.names=[qq.names(1:3)]; 
%qq.names=[qq.names(1:3), {'FAQ (Partner)'}]; 
qq.colorOrder=[qq.colorOrder(1:3), 10]; 

q1visit1 = (T1.NACCVNUM==1)+0; q1visit1(q1visit1==0)=NaN;
q1 = q1 + 0*sum(horzcat(qq.dat0{:}),2) + 0*q1visit1;

N=length(qq.names); 
q=[q1, horzcat(qq.dat0{:})]; zTOT=rocMaximizeLineSearchA(q,dx0a,zTOT.w,0.2); xBuild.Aroc(1)=zTOT.auc
i = isfinite(sum(q,2)) & isfinite(dx0a) & (dx0a<2); disp([sum(i), size(w0)])
[fpr,tpr,auc] = roc1(q(i,:)*zTOT.w1',dx0a(i)); auc=1-auc
qwa=q.*zTOT.w1; disp(100*nanvar(qwa)./sum(nanvar(qwa)))
%%
for k=0:N, q2=[{q1},qq.dat0]; q=[horzcat(q2{setdiff(1:N+1,k+1)})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.2); xBuild.Aroc(k+2)=z.auc; zz1{k+1}=z; end; disp(['w0=0.2 MOCA combo', num2str(k1)]), z, toc
for k=0:N, q2=[{q1},qq.dat0]; q=[horzcat(q2{setdiff(1:N+1,k+1)})]; z=rocMaximizeLineSearchA(q,dx0a,[],0.02); xBuild.Aroc(k+2)=max(z.auc,xBuild.Aroc(k+2)); zz1{k+1}=z; end; disp(['w0=0.02 MOCA combo', num2str(k1)]), z, toc
figure;
subplot(1,2,2); hold on; grid
colours1 = [GREY; colours(qq.colorOrder,:)]; 
xBuild.Aroc1 = xBuild.Aroc(1) - xBuild.Aroc;
qx = xBuild.Aroc1(2:end); [~,is] = sort(qx,'descend');
h = bar(1:N+1,diag(qx(is)),'stacked');
for k=1:length(h), h(k).FaceColor=colours1(is(k),:); end
qq.names1 = [{'MoCA-8'},qq.names]; legend(h,qq.names1(is))
xlabel('Candidate Items'); ylabel('Area under the ROC curve lost by removal'); shg 
xticks([]); yticks([0:.002:0.018]); yticklabels(num2str([0:.002:0.018]'));
xlim ([0.3 N+1+0.7]); ylim([-.001, .0195]); 
grid

%%  Build-up analysis for  cogbat
tic;
q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2x;  
qq.dat0=[{mocaMem5_20}, qq.dat0, {-x2.faq1.dat0}, {q1}]; 
qq.names=[{'Additional MOCA-5'}, qq.names, {'FAQ (Partner)'}, {'Moca-8'}]; 
qq.colorOrder=[6, qq.colorOrder, 10, 3]; 
q1 = q1+0*sum(horzcat(qq.dat0{:}),2);
q = horzcat(qq.dat0{:}); zTOT=rocMaximizeLineSearchA(q,dx0a,w0,0.2); xBuild.Aroc(1)=zTOT.auc
N=length(qq.names); 
ii = isfinite(sum(q,2)) + 0.0; ii(ii==0) = NaN; 

remainingItems = 1:N; orderedItems =[];  qauc1 = zeros(1,N); qauc=-Inf*ones(N,N);
for k1=1:N,
    for k=1:N+1-k1,
        q11 = ii.*horzcat(qq.dat0{[orderedItems,remainingItems(k)]}); 
        if k1<=2, zzz{k}=rocMaximizeLineSearchA(q11,dx0a,[],0.2); 
        else zzz{k}=rocMaximizeLineSearchA(q11,dx0a,zp.w,0.2);
        end
        qauc(k,k1)=zzz{k}.auc;
    end
    [quac1(k1),imax] = max(qauc(:,k1)); 
    orderedItems(k1) = remainingItems(imax); %[k1,size(orderedItems)]
    remainingItems = setdiff(remainingItems, remainingItems(imax)); 
    zp = zzz{imax};
end 
figure; subplot(122); hold on; grid on;
for k1=1:N, plot(k1,quac1(k1),'.','MarkerSize',36,'Color',colours(qq.colorOrder(orderedItems(k1)),:)); end
toc
legend( qq.names(orderedItems) )
title ('Cog Battery Items +  FAQ & a 2nd MoCA recall')
qauc
q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2x;  
qq.dat0=[qq.dat0, {q1}]; 
qq.names=[qq.names,  {'Moca-8'}]; 
qq.colorOrder=[qq.colorOrder, 3]; 
q1 = q1+0*sum(horzcat(qq.dat0{:}),2);
q = horzcat(qq.dat0{:}); zTOT=rocMaximizeLineSearchA(q,dx0a,w0,0.2); xBuild.Aroc(1)=zTOT.auc
N=length(qq.names); 
ii = isfinite(sum(q,2)) + 0.0; ii(ii==0) = NaN; 

remainingItems = 1:N; orderedItems =[];  qauc1 = zeros(1,N); qauc=-Inf*ones(N,N);
for k1=1:N,
    for k=1:N+1-k1,
        q11 = ii.*horzcat(qq.dat0{[orderedItems,remainingItems(k)]}); 
        if k1<=2, zzz{k}=rocMaximizeLineSearchA(q11,dx0a,[],0.2); 
        else zzz{k}=rocMaximizeLineSearchA(q11,dx0a,zp.w,0.2);
        end
        qauc(k,k1)=zzz{k}.auc;
    end
    [quac1(k1),imax] = max(qauc(:,k1)); 
    orderedItems(k1) = remainingItems(imax); %[k1,size(orderedItems)]
    remainingItems = setdiff(remainingItems, remainingItems(imax)); 
    zp = zzz{imax};
end 
subplot(121); hold on; grid on;
for k1=1:N, plot(k1,quac1(k1),'.','MarkerSize',36,'Color',colours(qq.colorOrder(orderedItems(k1)),:)); end
toc
legend( qq.names(orderedItems) )
title ('Cog Battery Items alone')

%%

q1 = x2.buildBase.dat0{1}; 
qq = x2.cb.MOCA2s1;
% qq.dat0=[{mocaMem5_20}, qq.dat0, {-x2.faq1.dat0}]; 
% qq.names=[{'MOCA-5 (Recall)'}, qq.names, {'FAQ (Partner)'}]; 
% qq.colorOrder=[6, qq.colorOrder, 10]; 
clear faqNames; 
for k=1:10, faqNames{k} = ['F',num2str(k)]; end 
for k=1:10, faqData{k} = -x2.faq1.dat0(:,k); end 
%qq.dat0=[qq.dat0([1,2,3]), faqData]; qq.dat0
%qq.dat0=[qq.dat0([1,2,4]), {-x2.faq1.dat0}]; qq.dat0
qq.dat0=[qq.dat0(1:3)]; 
qq.names=[qq.names(1:3), faqNames]; 
qq.names=[qq.names(1:3)]; 
%qq.names=[qq.names(1:3), {'FAQ (Partner)'}]; 
qq.colorOrder=[qq.colorOrder(1:3), 10]; 

q1visit1 = (T1.NACCVNUM==1)+0; q1visit1(q1visit1==0)=NaN;
q1 = q1 + 0*sum(horzcat(qq.dat0{:}),2) + 0*q1visit1;

N=length(qq.names); 
q=[q1, horzcat(qq.dat0{:})]; zTOT=rocMaximizeLineSearchA(q,dx0a,zTOT.w,0.2); xBuild.Aroc(1)=zTOT.auc
i = isfinite(sum(q,2)) & isfinite(dx0a) & (dx0a<2); disp([sum(i), size(w0)])
[fpr,tpr,auc] = roc1(q(i,:)*zTOT.w1',dx0a(i)); auc=1-auc
qwa=q.*zTOT.w1; disp(100*nanvar(qwa)./sum(nanvar(qwa)))











% xBar=[1:3, 5:6, 8:10, 12:13 15:17]; h=bar(xBar, qqq(2:end)-qqq(1));
% xticks([]);
% colours=colororder('glow12'); c = colours(x2.cb.MOCA1.colorOrder([5,5,6,5,6,5,6,6,1,2,1,2,2]),:); 
% c(8,:) = 0.5*mean(c(6:7,:))+0.5*[1 1 1];  % make the color for combo-bar a lightened version of the mean color
% c(13,:) = 0.5*mean(c(11:12,:))+0.5*[1 1 1]; % make the color for combo-bar a lightened version of the mean color
% h.FaceColor='flat'; h.CData=c; shg
% 





% for k=1:N, h(k) = plot(N-k+[1:2], xBuild.Aroc(N-k+[1:2]), '.-', 'Color', colours(qq.colorOrder(N-k+1),:), 'MarkerSize', 30, 'linewidth',1);  hold on; end
% h0(k1) = plot(1, xBuild.Aroc(1), '.', 'Marker', symbol1{k1},  'MarkerSize', 12, 'MarkerFaceColor', GREY, 'MarkerEdgeColor', GREY);
% plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid; hold on;
   

% xlim([0.5 N+1.5]); xticks([0:N]+1); xticklabels([0:N]);

% qq = x2.cb.MOCA2s; 
% qq.dat0=[{mocaMem5_20}, qq.dat0, {-x2.faq1.dat0}]; 
% qq.names=[{'MOCA-5'}, qq.names, {'FAQ'}]; 
% qq.colorOrder=[6, qq.colorOrder, 10]; 
% q1 = q1+0*sum(horzcat(qq.dat0{:}),2);
% k=[1:6]; qq.dat0=qq.dat0(k); qq.names=qq.names(k);  qq.colorOrder=qq.colorOrder(k);



%% Analysis of the "Press sign". Can the diff between CF & LF better predict MCI than CF alone?  

c1 = [-0.1:0.01:1]; ii = isfinite(qmoca1) & (qmoca1 <=30) & (qmoca1 >0) & (dx0a < 2);
clear rr rrp rrn rrpr rrnr rrd rrdp;
for k=1:length(c1), rr(k)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1) + 0*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end
for k=1:length(c1), rrd(k)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1) - c1(k)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end
for k=1:length(c1), rrdp(k)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1) + c1(k)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end
c2 = [1:20]/100;
for k2=1:length(c2), for k=1:length(c1), rrp(k,k2)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1)+c2(k2)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end; end
for k2=1:length(c2), for k=1:length(c1), rrn(k,k2)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1)-c2(k2)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end; end
c3 = [1:9]/10;
for k2=1:length(c3), for k=1:length(c1), rrpr(k,k2)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1)+c1(k)*c3(k2)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end; end
for k2=1:length(c3), for k=1:length(c1), rrnr(k,k2)=roc1(-[qmoca1(ii)+c1(k)*x2.cb.MOCA1.dat0(ii,1)-c1(k)*c3(k2)*x2.cb.MOCA1.dat0(ii,5)],dx0a(ii)); end; end
r0 = rr(c1==0);
figure; 
subplot(131); 
hh = plot(c1,rr-r0,'k','LineWidth',2); hold on; 
hd = plot(c1,rrd-r0,'r','LineWidth',2); hold on; 
hdp = plot(c1,rrdp-r0,'b','LineWidth',2); hold on;
legend([hh,hd(1),hdp(1)],{'CF-only','CF-LF difference (''Press sign'')','CF-LF sum'})
xlim([-.1 1]); ylim([-.025, .015]); 
title('Comparison of differences and sums for combining CF and LF data')

subplot(133); 
hh = plot(c1,rr-r0,'k','LineWidth',2); hold on; 
hp = plot(c1,rrp-r0,'b'); hn = plot(c1,rrn-r0,'r'); 
legend([hh,hp(1),hn(1)],{'No LF weighting','+LF weighting','-LF weighting'})
title('Comparison of partial independent differences and sums (0.02-0.2 for LF)')

subplot(132); 
hh = plot(c1,rr-r0,'k','LineWidth',2); hold on; 
hp = plot(c1,rrpr-r0,'b'); hn = plot(c1,rrnr-r0,'r'); 
hd = plot(c1,rrd-r0,'r','LineWidth',2); hold on; 
hdp = plot(c1,rrdp-r0,'b','LineWidth',2); hold on;
legend([hh,hp(1),hn(1)],{'No LF weighting','+LF weighting','-LF weighting'})
title('Comparison of partial relative differences and sums')

for k=1:3, subplot(1,3,k)
    xlabel('CF weighting');
    if k==1, ylabel('AUC improvement (wrt MoCA-29 w/o LF)'); end
    xlim([-.1 1]); ylim([-.025, .015]); 
    if k>1, yticklabels([]); end
end

subplot(133); 
hai = axes('Position',[2.35/3 .2 .35/3 .35]);
plot(c1,rr-r0,'k','LineWidth',2); 
hold on; plot(c1,rrp-r0,'b'); hold on; plot(c1,rrn-r0,'r'); 
xlim([.15 .37]); ylim([.012 .015]); title('Zoomed in at peak')

subplot(132); 
hai = axes('Position',[1.35/3 .2 .35/3 .35]);
plot(c1,rr-r0,'k','LineWidth',2); 
hold on; plot(c1,rrpr-r0,'b'); hold on; plot(c1,rrn-r0,'r'); 
xlim([.15 .37]); ylim([.012 .015]); title('Zoomed in at peak')
% exportgraphics(gcf, 'paper1_pressSignFailure.pdf', 'ContentType', 'vector');

%%




    %%
qD=qA; qD{1} = x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]); 
xBuild.Dvars = xBuild.Avars; xBuild.Dvars{1}='MoCA-30 (reweighted)';


%%
for k=1:N, qq=horzcat(qA{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Aroc(k)=z4.auc; end
figure; plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid
xlim(0.5+[0 N]); xticks(1:N); xticklabels(xBuild.Avars(1:N));
z4.w
toc
%%
qB=qA; qB{1} = [x2.moca1.datSum0, x2.moca1.datSum0]; 
xBuild.Bvars = xBuild.Avars; xBuild.Bvars{1}='MoCA-30';
for k=1:N, qq=horzcat(qB{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Broc(k)=z4.auc; end
hold on; plot(xBuild.Broc, '.-','MarkerSize',24); 
z4.w




%%
tic;
figure; set(gcf,'Position',[0 0 1160 800])
subplot(2,5,[1 2]); plotThresholdedAUCs(cbForMOCA,[0 10 20 40 60], x2.moca1.datSum0); title('MoCA Supplementation'); toc
subplot(2,5,[6 7]); plotThresholdedAUCs(cbForMMSE,[0 10 20 40 60], x2.mmse1.datSum0); title('MMSE Supplementation'); toc

%subplot(2,5,[4 5]); plotThresholdedAUCs(cbForMOCA,[0 10 20 40 60 80 90 100], x2.moca1.datSum0); title('MoCA Supplementation'); toc
%subplot(2,5,[9 10]); plotThresholdedAUCs(cbForMMSE,[0 10 20 40 60 80 90 100], x2.mmse1.datSum0); title('MMSE Supplementation'); toc

% subplot(2,5,[1 2]); ylim1=get(gca,'ylim'); plot(0.05*[1 1],ylim1,':k'); plot(0.125*[1 1],ylim1,':k'); 
% h=text(0.05,ylim1(end)*0.95,'median sensitivity (MoCA items)','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
% h=text(0.125,ylim1(end)*0.95,'median sensitivity (MMSE items)','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
% 
% subplot(2,5,[6 7]); ylim1=get(gca,'ylim'); plot(0.05*[1 1],ylim1,':k'); plot(0.125*[1 1],ylim1,':k'); 
% h=text(0.05,ylim1(end)*0.95,'median sensitivity (MoCA items)','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
% h=text(0.125,ylim1(end)*0.95,'median sensitivity (MMSE items)','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
% shg
%%
plotItemAUC(x2.faq1); title('FAQ')
%%
plotItemAUC(x2.cb1a); title('cb1a')
plotItemAUC(x2.cb2a); title('cb2a')
%%
plotItemAUC(x2.cb1b); title('cb1b')
%%
plotItemAUC2(x2.cb1b); title('cb1b')
%%
figure; 
plotThresholdedAUCs(cbForMOCA,[10 20 40 60]); title('cb1b')
%%
plotThresholdedAUCs(x2.cb1b,[10 20 40 60]); title('cb1b')
%%
plotThresholdedAUCs(x2.cb2b,[10 20 40 60]); title('cb1b')
%%
plotItemCombinationAUC(x2.cb1b); title('cb1b')
%%
plotItemCombinationPairsAUC(x2.cb1b.dat0(:,2),x2.cb1b)
%%
qq=T1.MOCATOTS;
qq(qq<0 | qq>=88) = NaN;
plotItemCombinationPairsAUC(qq,x2.cb1b)
%%
qq=T1.MOCATOTS;
qq(qq<0 | qq>=88) = NaN;
qq = qq - 25.8*mean(x2.faq1.dat0,2) + 0.3*mean(x2.cb1b.dat0(:,6:7),2) + 0*mean(x2.cb1b.dat0(:,8:9),2);
clear qq1
qq1.dx0 = x2.cb1b.dx0; 
qq1.dat0 = [x2.cb1b.dat0, x2.cb2b.dat0, [-mean(x2.faq1.dat0,2),mean(x2.cb1b.dat0(:,6:7),2),mean(x2.cb1b.dat0(:,8:9),2)]];
qq1.vars = {x2.cb1b.vars{:}, x2.cb2b.vars{:},'FTOT','A+V','TA+TB'};

plotItemCombinationPairsAUC(qq,qq1)

%%
practiceEffectAnalysis(x2.moca1)
practiceEffectAnalysis(x2.mmse1)
%%
%qq=T1.MOCATOTS;
%qq=T1.NACCMOCA;
qq=T1.NACCMMSE; qq(qq<0 | qq>=88) = NaN;
% z=rocMaximize(-[qq,x2.cb1b.dat0(:,6:9)],x2.cb1b.dx0);
% z
dx0=x2.cb1b.dx0; 
qq1 = [qq,qq, x2.cb1b.dat0(:,6)]; %, x2.cb1b.dat0(:,3),x2.cb1b.dat0(:,1),-round(x2.cb1b.dat0(:,8:9)/10),-10*x2.faq1.dat0]; %,-10*sum(x2.faq1.dat0,2)];
% qq1 = [qq,qq, x2.cb1b.dat0(:,6:7), -round(x2.cb1b.dat0(:,8:9)/10), x2.cb1b.dat0(:,3), x2.cb2b.dat0(:,[2])]; %,-10*sum(x2.faq1.dat0,2)];
% qq1(:,4:5) = -round(qq1(:,4:5)/10);
% qq1(:,2:5) = 0*qq1(:,2:5);
% qq1(:,6:end) = -qq1(:,6:end)*10;
% for k=1:size(qq1,2),
%    z = rocPolyFit(qq1(:,k), dx0, [], 12);
% end
tic;
z = rocMaximizeLineSearch(qq1,dx0)
toc
%%
dx0=x2.cb1b.dx0; 
i = x2.moca1.i;
x2.moca1.datMIS_ = 1*x2.moca1.datMIS1+2*x2.moca1.datMIS2+3*x2.moca1.datMIS3;
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datMIS_,x2.moca1.datMIS_]; z3 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,x2.moca1.datSumNoMem]; z1 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum,x2.moca1.datSum]; z2 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem, 0*x2.moca1.datMIS_]; z1 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum, 0*x2.moca1.datMIS_]; z2 = rocMaximizeLineSearch(qq,dx0),
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem, 2*x2.moca1.datMIS_]; z1 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum, 2*x2.moca1.datMIS_]; z2 = rocMaximizeLineSearch(qq,dx0),
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem, 2*x2.moca1.datMIS_]; z1 = rocMaximizeLineSearch(qq,dx0), 
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum, 2*x2.moca1.datMIS_]; z2 = rocMaximizeLineSearch(qq,dx0),
%% Move this to readMOCA!!
i = x2.moca1.i;
qMIS = [3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1];
qq=NaN(L,3); qq(i,1:3)=qMIS; 
x2.moca1.datMIS0all = qq;
x2.moca1.datMIS0 = sum(qq,2);

%%
x2.moca1.datMIS_ = 1*x2.moca1.datMIS1+2*x2.moca1.datMIS2+3*x2.moca1.datMIS3;
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datMIS_,x2.moca1.datMIS_]; z3 = rocMaximizeLineSearch(qq,dx0), z3.maxk
%%
x2.moca1.datMIS_ = 3*x2.moca1.datMIS3+2*x2.moca1.datMIS2+1*x2.moca1.datMIS1;
%%
qq=NaN(L,4); qq(i,1:4)=[x2.moca1.datSumNoMem,3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; z4 = rocMaximizeLineSearch(qq,dx0), z4.maxk
qq=NaN(L,4); qq(i,1:4)=[x2.moca1.datSum,3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; z4 = rocMaximizeLineSearch(qq,dx0), z4.maxk
%%
i = x2.moca1.i;
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[qq,0*cbForMOCA.dat(:,9)];  z4 = rocMaximizeLineSearch(qq,dx0)
qq=NaN(L,1); qq(i,1:1)=[x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0)
qq=NaN(L,1); qq(i,1:1)=[x2.moca1.datMIS_]; qq=[cbForMOCA.dat(:,9),0*qq];   z4 = rocMaximizeLineSearch(qq,dx0)
%%
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
%% adding MIS (MIS_) or regular recall (MIS3) or CRAFT story data to MOCA or MocaNoMem
i = x2.moca1.i; 
qq=NaN(L,3); qq(i,1:3)=[x2.moca1.datSum,x2.moca1.datSum,0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[x2.moca1.datSumNoMem,3*x2.moca1.datMIS3,0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,4); qq(i,1:4)=[x2.moca1.datSumNoMem,6*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[qq,0*cbForMOCA.dat(:,9)];  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[[1 0].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[[1 1].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
%% MIS (MIS_) or regular recall (MIS3) or CRAFT story data alone
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
dx0a(T1.VISITYR>2016 | (T1.VISITYR==2016 & T1.VISITMO<3))=NaN;
i = x2.moca1.i;
%%
i = x2.moca1.i; L=length(dx0);
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qq=NaN(L,2); qq(i,1:2)=[3*x2.moca1.datSum, 3*x2.moca1.datSum];   disp('Regular MOCA (all data)'); z4 = rocMaximizeLineSearch(qq,dx0a)
z4 = rocMaximizeLineSearch([x2.moca1.datSum0,x2.moca1.datSum0],dx0a)
z4 = rocMaximizeLineSearch(x2.moca1.datSum0*[1 1],dx0a)

%%
L=length(dx0a);
qq=NaN(L,2); qq(i,1:2)=[3*x2.moca1.datMIS3, 3*x2.moca1.datMIS3];   disp('Normal recall (all data)'); z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];  disp('Normal recall (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[qq,0*cbForMOCA.dat(:,9)]; disp('Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[2*x2.moca1.datMIS_, 0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   disp('Fixed MIS (MIS available)'); z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[2*x2.moca1.datMIS_, 0*x2.moca1.datMIS_]; qq=[[0].*qq, cbForMOCA.dat(:,9)];   disp('CRAFT (MIS available)'); z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[2*x2.moca1.datMIS_, 0*x2.moca1.datMIS_]; qq=[[1].*qq, cbForMOCA.dat(:,9)];   disp('CRAFT + Fixed MIS (MIS available)'); z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,9), 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
%%
qq=NaN(L,2); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,9), 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1]), 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
%% effect of adding FAQ!!
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], -24*x2.faq1.dat0, 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], -24*x2.faq1.dat0, 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], -24*x2.faq1.dat0, 20*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], -24*x2.faq1.dat0, 20*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
%% build-up w/o FAQ
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 3*qq]; disp('CRAFT + Benton + CatFlu + Flexible MIS (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + Flexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30,17:19]), 3*qq]; disp('CRAFT + Benton + CatFlu + Flexible MIS + MocaOrientation++ (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 3*qq]; disp('CRAFT + Benton + CatFlu + Recall-only (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + Recall-only + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a)
%%
tic;
N=5; qq=horzcat(qA{1:N}); z5 = rocMaximizeLineSearch(qq,dx0a); z5, toc
N=6; qq=horzcat(qA{1:N}); z6 = rocMaximizeLineSearch(qq,dx0a); z6, toc
%%
tic
N=5; qq=horzcat(qA{1:N}); i = isfinite(sum(qq,2)) & dx0a<2; 1-roc1(qq(i,:)*[1,z5.dw1]',dx0a(i))
K=10;
z5a = rocMaximizeLineSearch([qq,K*x2.im.dat0(:,[1,3,4])],dx0a); z5a, 
z5b = rocMaximizeLineSearch([qq,K*x2.im.dat0(:,[2,4,6])],dx0a); z5b, 
z5a1 = rocMaximizeLineSearch([qq*[1,z5.dw1]',K*x2.im.dat0(:,[1,3,4])],dx0a); z5a1, 
z5b1 = rocMaximizeLineSearch([qq*[1,z5.dw1]',K*x2.im.dat0(:,[2,4,6])],dx0a); z5b1, 
%%
K=[0.1, 0.01, 0.01];
N=6; qq=horzcat(qA{1:N}); i = isfinite(sum(qq,2)) & dx0a<2; 1-roc1(qq(i,:)*[1,z6.dw1]',dx0a(i))
z6a = rocMaximizeLineSearch([qq,K.*x2.im.dat0(:,[1,3,4])],dx0a); z6a, 
z6b = rocMaximizeLineSearch([qq,K.*x2.im.dat0(:,[2,4,6])],dx0a); z6b, 
z6a1 = rocMaximizeLineSearch([qq*[1,z6.dw1]',K.*x2.im.dat0(:,[1,3,4])],dx0a); z6a1, 
z6b1 = rocMaximizeLineSearch([qq*[1,z6.dw1]',K.*x2.im.dat0(:,[2,4,6])],dx0a); z6b1, 

%%
tic;
clear qA;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qA{1} = [x2.moca1.datMIS0all(:,1), 5*x2.moca1.dat0(:,[25,28,29,30])]; xBuild.Avars{1}='MoCA-9_(_R_O_)';
qA{2} = x2.moca1.datMIS0all(:,2:3); xBuild.Avars{2}='+MoCA-MIS';
qA{3} = cbForMOCA.dat(:,8);         xBuild.Avars{3}='+Craft';
qA{4} = cbForMOCA.dat(:,2:3).*[1 2];    xBuild.Avars{4}='+CatFluency';
qA{5} = cbForMOCA.dat(:,1).*4;  xBuild.Avars{5}='+BentonFig';
qA{6} = -24*x2.faq1.dat0;  xBuild.Avars{6}='+FAQ';
N=6;
for k=1:N, qq=horzcat(qA{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Aroc(k)=z4.auc; end
figure; plot(xBuild.Aroc, '.-','MarkerSize',24, 'linewidth',2); ylim([0.8 0.92]); grid
xlim(0.5+[0 N]); xticks(1:N); xticklabels(xBuild.Avars(1:N));
z4.w
toc
%%
qB=qA; qB{1} = [x2.moca1.datSum0, x2.moca1.datSum0]; 
xBuild.Bvars = xBuild.Avars; xBuild.Bvars{1}='MoCA-30';
for k=1:N, qq=horzcat(qB{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Broc(k)=z4.auc; end
hold on; plot(xBuild.Broc, '.-','MarkerSize',24); 
z4.w

i2 = 6; z4 = rocMaximizeLineSearch(x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]),dx0a);
w(1:5)=z4.w(1:5); w(6)=1; w(7:30)=z4.w(6:29); moca_weighted = x2.moca1.dat0*w(:);
qC=qA; qC{1} = [moca_weighted, moca_weighted]; 
xBuild.Cvars = xBuild.Avars; xBuild.Cvars{1}='MoCA-30 (weighted)';
for k=1:N, qq=horzcat(qC{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Croc(k)=z4.auc; end
hold on; plot(xBuild.Croc, '.-','MarkerSize',24); 
z4.w
% 
% i2 = 6; z4 = rocMaximizeLineSearch(x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]),dx0a);
% w(1:5)=z4.dw1(1:5); w(6)=1; w(7:30)=z4.dw1(6:29); moca_weighted = x2.moca1.dat0*w(:);
% qC=qA; qC{1} = 0.25*[moca_weighted, moca_weighted]; 
% xBuild.Cvars = xBuild.Avars; xBuild.Cvars{1}='MoCA-30 (weighted)';
% for k=1:N, qq=horzcat(qC{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Croc(k)=z4.auc; end
% hold on; plot(xBuild.Croc, '.-','MarkerSize',24); 
% z4.dw1
% toc

i2 = 6; 
qD=qA; qD{1} = x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]); 
xBuild.Dvars = xBuild.Avars; xBuild.Dvars{1}='MoCA-30 (reweighted)';
for k=1:N, qq=horzcat(qD{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Droc(k)=z4.auc; end
hold on; plot(xBuild.Droc, '.--','MarkerSize',24); 
z4.w
toc
%%
tic;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qA{1} = [x2.mmse1.dat0]; xBuild.Avars{1}='MMSE-30';
qA{2} = cbForMMSE.dat(:,6:7);    xBuild.Avars{2}='+CatFluency';
qA{3} = cbForMMSE.dat(:,10).*1;  xBuild.Avars{3}='+BentonFig';
qA{4} = cbForMMSE.dat(:,1); xBuild.Avars{4}='+LM';
qA{5} = cbForMMSE.dat(:,4:5).*1;  xBuild.Avars{5}='+DigitSpan';
qA{6} = cbForMMSE.dat(:,8:9).*1;  xBuild.Avars{6}='+TrailsA+B';
%qA{7} = cbForMMSE.dat(:,10).*1;  xBuild.Avars{7}='+BNT';
qA{7} = -24*x2.faq1.dat0;  xBuild.Avars{7}='+FAQ';
N=7;
for k=1:N, qq=horzcat(qA{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Aroc(k)=z4.auc; end
%%
figure; 
hold on; plot(xBuild.Aroc, '.-','MarkerSize',24); ylim([0.75 0.92]); grid
xlim(0.5+[0 N]); xticks(1:N); xticklabels(xBuild.Avars(1:N));
z4.dw1
toc
%%
qB=qA; qB{1} = [x2.moca1.datSum0, x2.moca1.datSum0]; 
xBuild.Bvars = xBuild.Avars; xBuild.Bvars{1}='MoCA-30';
for k=1:N, qq=horzcat(qB{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Broc(k)=z4.auc; end
hold on; plot(xBuild.Broc, '.-','MarkerSize',24); 
z4.dw1

i2 = 6; z4 = rocMaximizeLineSearch(x2.moca1.dat0(:,[i2,setdiff(1:30,i2)]),dx0a);
w(1:5)=z4.dw1(1:5); w(6)=1; w(7:30)=z4.dw1(6:29); moca_weighted = x2.moca1.dat0*w(:);
qC=qA; qC{1} = [moca_weighted, moca_weighted]; 
xBuild.Cvars = xBuild.Avars; xBuild.Cvars{1}='MoCA-30 (weighted)';
for k=1:N, qq=horzcat(qC{1:k}); z4 = rocMaximizeLineSearch(qq,dx0a); xBuild.Croc(k)=z4.auc; end
hold on; plot(xBuild.Croc, '.-','MarkerSize',24); 
z4.dw1

%%
rocMaximizeLineSearch(x2.moca1.dat0(:,24:-1:20),dx0a)
%%
tic;
i2 = 6;
%z4 = rocMaximizeLineSearch(x2.moca1.dat0,dx0a), z4.dw1
%z4 = rocMaximizeLineSearch([x2.moca1.dat0(:,setdiff(1:30,i2)), 10*x2.moca1.dat0(:,i2)],dx0a), z4.dw1
z4 = rocMaximizeLineSearch([x2.moca1.dat0(:,i2), x2.moca1.dat0(:,setdiff(1:30,i2))],dx0a)
%z4 = rocMaximizeLineSearch([x2.moca1.dat0(:,i2)],dx0a); %,
%x2.moca1.dat0(:,setdiff(1:30,i2))],dx0a),
toc
%% look at the effect of removing each individual contribution 
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + Flexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(1)=z4.auc; xFinal.vars{1}='Full'; 
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[1:3]).*[2.5 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('xxCRAFT + Benton + CatFlu + Flexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(2)=z4.auc; xFinal.vars{2}='CraftIm'; 
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,2:3]).*[1 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + xxBenton + CatFlu + Flexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(3)=z4.auc; xFinal.vars{3}='Benton'; 
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1]).*[1 4], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + xxCatFlu + Flexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a);     xFinal.roc(4)=z4.auc; xFinal.vars{4}='CatFluency'; 
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + xxFlexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(5)=z4.auc; xFinal.vars{5}='MoCA-MIS'; 
qq=NaN(L,3); qq(i,1:3)=[0*x2.moca1.datMIS3,0*x2.moca1.datMIS2,0*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 10*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + xxxFlexible MIS + MocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(6)=z4.auc; xFinal.vars{6}='MoCA-Recall'; 
qq=NaN(L,3); qq(i,1:3)=[3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[cbForMOCA.dat(:,[9,1:3]).*[1 4 1 2], 0*x2.moca1.dat0(:,[25,28,29,30]), 3*qq]; disp('CRAFT + Benton + CatFlu + Flexible MIS + xxMocaOrientation (MIS available)');  z4 = rocMaximizeLineSearch(qq,dx0a); xFinal.roc(7)=z4.auc; xFinal.vars{7}='MoCA-9_R_O-Base'; 
figure;
iReorder = [6,4,5,1,2,3];
plot(xFinal.roc(1)-xFinal.roc(1+iReorder), '.-','MarkerSize',24);
xticks(1:length(xFinal.roc)-1); xticklabels(xFinal.vars(1+iReorder));
%% Better coding for above...
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qA{1} = [x2.moca1.datMIS0all(:,1), 0*x2.moca1.datMIS0all(:,1)]; xBuild.Avars{1}='MoCA-Recall_5';
qA{2} = x2.moca1.dat0(:,[25,28,29,30]).*10; xBuild.Avars{2}='MoCA-Orientation_4';
qA{3} = x2.moca1.datMIS0all(:,2:3); xBuild.Avars{3}='+MoCA-MIS';
qA{4} = cbForMOCA.dat(:,9);         xBuild.Avars{4}='+Craft';
qA{5} = cbForMOCA.dat(:,2:3).*[1 2];    xBuild.Avars{5}='+CatFluency';
qA{6} = cbForMOCA.dat(:,1).*4;  xBuild.Avars{6}='+BentonFig';
qA{7} = cbForMOCA.dat(:,4:5).*1;  xBuild.Avars{7}='+TrailsA+B';
qA{8} = cbForMOCA.dat(:,6:7).*0.1;  xBuild.Avars{8}='+VerFluency';
qA{9} = cbForMOCA.dat(:,10:11).*1;  xBuild.Avars{9}='+BentonFig';
qA{10} = cbForMOCA.dat(:,12).*0.05;  xBuild.Avars{10}='+Mint';
tic;
N=length(qA);
for k=1:N, qq=horzcat(qA{setdiff(1:N,k)}); z4 = rocMaximizeLineSearch(qq,dx0a); xUnBuild.Aroc(k)=z4.auc; disp([k,round(100*z4.dw1)]); end
k=1+N;     qq=horzcat(qA{setdiff(1:N,[])}); z4 = rocMaximizeLineSearch(qq,dx0a); xUnBuild.Aroc(k)=z4.auc; disp([k,round(100*z4.dw1)]);
xUnBuild.Aroc
toc
figure; subplot(212)
iReorder = [2,3,1,4:N];
plot(xUnBuild.Aroc(11)-xUnBuild.Aroc(iReorder), '.-','MarkerSize',24);
xticks(1:N); xticklabels(xBuild.Avars(iReorder));
%%
clear qA;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qA{1} = [x2.moca1.datMIS0all(:,1)*0.1, x2.moca1.dat0(:,[25,28,29,30])].*10; xBuild.Avars{1}='MoCA-Orientation_4';
qA{2} = x2.moca1.datMIS0all(:,2:3); xBuild.Avars{2}='+MoCA-MIS';
qA{3} = cbForMOCA.dat(:,9);         xBuild.Avars{3}='+Craft';
qA{4} = cbForMOCA.dat(:,2:3).*[1 2];    xBuild.Avars{4}='+CatFluency';
qA{5} = cbForMOCA.dat(:,1).*4;  xBuild.Avars{5}='+BentonFig';
qA{6} = cbForMOCA.dat(:,4:5).*1;  xBuild.Avars{6}='+TrailsA+B';
qA{7} = cbForMOCA.dat(:,6:7).*0.1;  xBuild.Avars{7}='+VerFluency';
qA{8} = cbForMOCA.dat(:,10:11).*1;  xBuild.Avars{8}='+cDigitSpan';
qA{9} = cbForMOCA.dat(:,12).*0.05;  xBuild.Avars{9}='+Mint';
tic;
N=length(qA);
for k=2:N, qq=horzcat(qA{[1,k]}); z4 = rocMaximizeLineSearch(qq,dx0a); x1Build.Aroc(k)=z4.auc; disp([k,round(100*z4.dw1)]); end
k=1;     qq=horzcat(qA{[1]}); z4 = rocMaximizeLineSearch(qq,dx0a); x1Build.Aroc(k)=z4.auc; disp([k,round(100*z4.dw1)]);
x1Build.Aroc
toc
%%
figure;
iReorder = [2,3,1,4:N];
subplot(211); plot(-x1Build.Aroc(1)+x1Build.Aroc(2:end), '.-','MarkerSize',24);
xticks(1:N); xticklabels(xBuild.Avars(2:end));

%%
clear qA;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
qA{1} = [x2.moca1.datMIS0all(:,1)*0.1, x2.moca1.dat0(:,[25,28,29,30])].*10; xBuild.Avars{1}='MoCA-Orientation_4';
qA{2} = x2.moca1.datMIS0all(:,2:3); xBuild.Avars{2}='+MoCA-MIS';
qA{3} = cbForMOCA.dat(:,8);         xBuild.Avars{3}='+Craft';
qA{4} = cbForMOCA.dat(:,2:3).*[1 2];    xBuild.Avars{4}='+CatFluency';
qA{5} = cbForMOCA.dat(:,1).*4;  xBuild.Avars{5}='+BentonFig';
qA{6} = cbForMOCA.dat(:,4:5).*1;  xBuild.Avars{6}='+TrailsA+B';
qA{7} = cbForMOCA.dat(:,6:7).*0.1;  xBuild.Avars{7}='+VerFluency';
qA{8} = cbForMOCA.dat(:,9:10).*1;  xBuild.Avars{8}='+cDigitSpan';
qA{9} = cbForMOCA.dat(:,11).*0.05;  xBuild.Avars{9}='+Mint';
tic;
N=length(qA);
for k=2:N, qq=horzcat(qA{[1,k]}); z4 = rocMaximizeLineSearch(qq,dx0a); x1Build.Aroc(k)=z4.auc; disp([k,round(100*z4.w)]); end
k=1;     qq=horzcat(qA{[1]}); z4 = rocMaximizeLineSearch(qq,dx0a); x1Build.Aroc(k)=z4.auc; disp([k,round(100*z4.w)]);
x1Build.Aroc
toc
%%
figure;
iReorder = [2,3,1,4:N];
subplot(211); plot(-x1Build.Aroc(1)+x1Build.Aroc(2:end), '.-','MarkerSize',24);
xticks(1:N); xticklabels(xBuild.Avars(2:end));



%%
i = x2.moca1.i;
qq=NaN(L,3); qq(i,1:3)=[x2.moca1.datSum,x2.moca1.datSum,0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,3); qq(i,1:3)=[x2.moca1.datSumNoMem,3*x2.moca1.datMIS3,0*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,4); qq(i,1:4)=[x2.moca1.datSumNoMem,3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[qq,0*cbForMOCA.dat(:,9)];  z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[[1 0].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSumNoMem,2*x2.moca1.datMIS_]; qq=[[1 1].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0a)
%%
i = x2.moca1.i;
qq=NaN(L,4); qq(i,1:4)=[x2.moca1.datSum,3*x2.moca1.datMIS3,2*x2.moca1.datMIS2,1*x2.moca1.datMIS1]; qq=[qq,0*cbForMOCA.dat(:,9)];  z4 = rocMaximizeLineSearch(qq,dx0)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum,2*x2.moca1.datMIS_]; qq=[qq,0*cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum,2*x2.moca1.datMIS_]; qq=[[1 0].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0)
qq=NaN(L,2); qq(i,1:2)=[x2.moca1.datSum,2*x2.moca1.datMIS_]; qq=[[1 1].*qq, cbForMOCA.dat(:,9)];   z4 = rocMaximizeLineSearch(qq,dx0)


%% Paper 1 Fig 4c: 2/29/2024
% Make scatter plot of AUC improvement over MMSE/MOCA alone vs item AUC for adding each CogBat test!  
% Also show correlation coeff value a scatter point size, use color for each item name (colororder)

% i = dx0==1; %i01a = (dx0==0 | dx0==1);
x2.mmse1.totRaw0 = NaN*dx0; 
x2.mmse1.totRaw0(x2.mmse1.i)=x2.mmse1.totRaw;

for k1=1:2,
    if k1==1; q = T1.NACCMMSE; q(q<0 | q>=88) = NaN; q1=cbForMMSE; end
    if k1==2; q = T1.NACCMOCA; q(q<0 | q>=88) = NaN; q1=cbForMOCA; end

    for k=1:size(q1.dat,2),
        qq=[q,q1.dat(:,k)];
        i = ( (dx0==0 | dx0==1) & isfinite(sum(qq,2)) );
        q1.corr(k) = corr(qq(i,1),qq(i,2),'rows','pairwise');
        z0 = roc1(qq(i,1),dx0(i));
        z1 = roc1(qq(i,2),dx0(i));
        z2 = rocMaximizeLineSearch(qq,dx0);
        q1.rocSingle(k) = z1;
        q1.rocDouble(k) = z2.auc;
        q1.rocPre(k) = z2.aucPre;
        q1.rocPre0(k) = z0;
        q1.rocDelta(k) = z2.auc - z2.aucPre;
    end
    if k1==1; cbForMMSE=q1; figure; end
    if k1==2; cbForMOCA=q1; end
    scatter(0.5+abs(q1.rocSingle-0.5),q1.rocDelta,200*abs(q1.corr.^2),'filled'); hold on;
end

%for k=1:size(cbForMOCA.dat,2), cbForMOCA.corr(k) = corr(cbForMOCA.dat(i,k),x2.moca1.totRaw(i),'rows','pairwise'); end


%Make scatter here
%%
z2 = rocMaximizeLineSearch(qq,dx0);z2
    
%% MoCA-MIS analysis: compare AUC's for: MOCA, mocaNoMem, MIS,  moca+mis, mocanomem+mis, also dissect mis into 3 vars


%% Make a plot showing # of MCI participants w/ data for each CogBat test AND MoCA (or MMSE)
figure; clear L;
for k1=1:3,
    if k1==1, qq=T1.NACCMMSE; end
    if k1==2, qq=T1.NACCMOCA; end
    if k1==3, qq=T1.MOCATOTS; end
    qq(qq<0 | qq>=88) = NaN;
    L(1,k1) = sum(isfinite(qq) & dx0==1); 
    for k=1:size(x2.cb1b.dat0,2), L(k+1,k1) = sum(isfinite(qq+x2.cb1b.dat0(:,k)) & dx0==1); end
end 
subplot(121); plot(L); xticks(1:k+1); xticklabels({'Alone',x2.cb1b.vars{:}}); title('x2.cb1b tests'); grid
clear L;
for k1=1:3,
    if k1==1, qq=T1.NACCMMSE; end
    if k1==2, qq=T1.NACCMOCA; end
    if k1==3, qq=T1.MOCATOTS; end
    qq(qq<0 | qq>=88) = NaN;
    L(1,k1) = sum(isfinite(qq) & dx0==1); 
    for k=1:size(x2.cb2b.dat0,2), L(k+1,k1) = sum(isfinite(qq+x2.cb2b.dat0(:,k)) & dx0==1); end
end
subplot(122); plot(L); xticks(1:k+1); xticklabels({'Alone',x2.cb2b.vars{:}}); title('x2.cb2b tests'); grid; legend({'MMSE','MoCA'})


% z = rocPolyFit(qq, dx0, 0:30, 12);


%%
qq=T1.MOCATOTS;
qq(qq<0 | qq>=88) = NaN;
% z=rocMaximize(-[qq,x2.cb1b.dat0(:,6:9),x2.faq1.dat0],x2.cb1b.dx0);
% z=rocMaximize(-[qq,1000*mean(x2.faq1.dat0,2)],x2.cb1b.dx0);
% z=rocMaximize(-[qq,1000 * x2.faq1.dat0],x2.cb1b.dx0,-0.024148*0.1*ones(1,10));
z=rocMaximize(-[qq,x2.cb1b.dat0(:,6:9),1000*mean(x2.faq1.dat0,2)],x2.cb1b.dx0,[0.12,0.11, -0.015, -0.015, -0.024148]);

z
%%
qq=T1.NACCMMSE;
qq(qq<0 | qq>=88) = NaN;
plotItemCombinationPairsAUC(qq,x2.cb1b)
%%

plotItemAUC(x2.cb2a); title('cb2a')
plotItemAUC(x2.cb2b); title('cb2b')

%%
plotItemCorr(x2.cb1a); title('cb1a');
%%
plotItemCorr(x2.cb1b); title('cb1b');

%% education Split Analaysis for MMSE & MoCA
mmseDx = x2.mmse1.dx; 
mocaDx = x2.moca1.dx; 

edSplitList = 13:20; clear r1 r2
for k=1:length(edSplitList),
    edSplit = edSplitList(k);
    roc = educationSplitAnalysis(x2.moca1.totRaw,x2.moca1.education,mocaDx,[0 1],edSplit); r1(k,:)=[roc.lowEd,roc.highEd,roc.tot];
    roc = educationSplitAnalysis(x2.mmse1.totRaw,x2.mmse1.education,mmseDx,[0 1],edSplit); r2(k,:)=[roc.lowEd,roc.highEd,roc.tot];
end
% figure; 
% subplot(121); plot(edSplitList,r1); title('MoCA - Education Split'); legend({'Low','High','Overall'})
% subplot(122); plot(edSplitList,r2); title('MMSE - Education Split'); legend({'Low','High','Overall'})
figure; subplot(111); hold on;
plot(edSplitList,r1,'-','linewidth',4); title('MoCA & MMSE - Education Split'); legend({'Low','High','Overall'})
plot(edSplitList,r2,'--','linewidth',2); title('MMSE - Education Split'); 
legend({'MoCA Low','MoCA High','MoCA Overall','MMSE Low','MMSE High','MMSE Overall'})
c = loadColors(); colororder([c.BLUE2; c.RED; c.LTPURPLE]) 


%% education Bin Analaysis for MMSE & MoCA
mmseDx = x2.mmse1.dx; 
mocaDx = x2.moca1.dx; 

edBins = [0, 11.5,12.5, 15.5, 16.5, 19.5, 20.5,50]; clear r1 r2
edBinLabels = {'<12','12','13-15','16','17-19','20','>20'};
r1 = educationBinAnalysis(x2.moca1.totRaw,x2.moca1.education,mocaDx,[0 1],edBins); 
r2 = educationBinAnalysis(x2.mmse1.totRaw,x2.mmse1.education,mmseDx,[0 1],edBins); 

figure; subplot(111); hold on;
n=length(edBins)-1;
plot(1:n,r1.bins,'.-','linewidth',4); plot(1:n,r1.tot*ones(1,n),'-','linewidth',4); 
plot(1:n,r2.bins,'.--','linewidth',2); plot(1:n,r2.tot*ones(1,n),'--','linewidth',2); 
title('MoCA & MMSE - Education Bins'); 
legend({'MoCA','MoCA Overall','MMSE','MMSE Overall'})
xticks(1:n); xticklabels(edBinLabels); xlabel('Education level (years)')
%%
c = loadColors(); colororder([c.BLUE2; c.LTPURPLE]) 
%%

qq=x2.moca1; i=ismember(qq.dx,[0 1]); [~,~,rocAUC] = roc1(30-qq.datSum(i),qq.dx(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.mmse1; i=ismember(qq.dx,[0 1]); [~,~,rocAUC] = roc1(30-qq.datSum(i),qq.dx(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.faq1;  i=ismember(qq.dx,[0 1]); [~,~,rocAUC] = roc1(30 - (30-qq.datSum(i)),qq.dx(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
toc
%%
qq=x2.moca1; i=ismember(qq.dxa,[0 1]); [~,~,rocAUC] = roc1(30-qq.totRaw(i),qq.dxa(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.mmse1; i=ismember(qq.dxa,[0 1]); [~,~,rocAUC] = roc1(30-qq.totRaw(i),qq.dxa(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.faq1;  i=ismember(qq.dxa,[0 1]); [~,~,rocAUC] = roc1(30 - (30-qq.totRaw(i)),qq.dxa(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
toc
%%
qq=x2.moca1; i=ismember(qq.dx0a,[0 1])&isfinite(qq.datSum0); [~,~,rocAUC] = roc1(30-qq.datSum0(i),qq.dx0a(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.mmse1; i=ismember(qq.dx0a,[0 1])&isfinite(qq.datSum0); [~,~,rocAUC] = roc1(30-qq.datSum0(i),qq.dx0a(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
%%
qq=x2.moca1; i=ismember(qq.dx0a,[0 1])&isfinite(qq.datSum0)&(T1.NACCVNUM==1); [~,~,rocAUC] = roc1(30-qq.datSum0(i),qq.dx0a(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
qq=x2.mmse1; i=ismember(qq.dx0a,[0 1])&isfinite(qq.datSum0)&(T1.NACCVNUM==1); [~,~,rocAUC] = roc1(30-qq.datSum0(i),qq.dx0a(i),[0:1:30],0); disp([rocAUC, sum(i)/1000])
toc
%%
tic
qq=x2.moca1; a = testMatchingDifferentDate (qq,qq,[], [],'less than'); i1=a.i1; i2=a.i2; 
toc
% [~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))-qq.recall(a.i1(i1))/2.56+30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])

i1=ismember(qq.dx(a.i1),[0 1]); 
i2=ismember(qq.dx(a.i2),[0 1]); 
%%  paper1 - "extra" moca-recall analyis shows improvement from 0.82 to 0.86 from adding additional (weighted) MOCA-recall info!!!
disp('MoCA + MoCA test combination (single/combined/single1/single2)')
rocAUC = roc1(30-qq.totRaw,qq.dx,[0:1:30],0); disp(rocAUC)
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.totRaw(a.i2(i2)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i2(i2)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])

qq.dx0a_ = dx0a(qq.i);
disp('MoCA + MoCA test combination (for aMCI) (single/combined/single1/single2)')
rocAUC = roc1(30-qq.totRaw,qq.dx0a_,[0:1:30],0); disp(rocAUC)
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.totRaw(a.i2(i2)),qq.dx0a_(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1)),qq.dx0a_(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i2(i2)),qq.dx0a_(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])

qq.recall = qq.datMIS0all(qq.i,1)/3;
disp('MoCA + (MoCA-Recall) test combination (for aMCI) (single/combined/single1/single2)')
rocAUC = roc1(30-qq.totRaw,qq.dx0a_); disp(rocAUC)
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
%%
qq.recall = qq.datMIS0all(qq.i,1)/1.83;
disp('MoCA + (MoCA-Recall [weighted]) test combination (for aMCI) (single/combined/single1/single2)')
rocAUC = roc1(30-qq.totRaw,qq.dx0a_); disp(rocAUC)
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])

qq.recall = qq.datMIS0all(qq.i,1)/1.83;
disp('MoCA + (MoCA-Recall [weighted+]) test combination (for aMCI) (single/combined/single1/single2)')
rocAUC = roc1(30-qq.totRaw,qq.dx0a_); disp(rocAUC)
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))-qq.recall(a.i1(i1))/2.56+30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))-qq.recall(a.i1(i1))/2.56,qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.recall(a.i2(i2)),qq.dx0a_(a.i1(i1))); disp([rocAUC, length(a.i1(i1))/1000])

%%
tic
qq=x2.mmse1; a = testMatchingDifferentDate (qq,qq,[], [],'less than'); i1=a.i1; i2=a.i2; 
toc
i1=ismember(qq.dx(a.i1),[0 1]); 
i2=ismember(qq.dx(a.i2),[0 1]); 
disp('MMMSE + MMMSE test combination')
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.totRaw(a.i2(i2)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1))+30-qq.totRaw(a.i2(i2)),qq.dx(a.i2(i2)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i1(i1)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
[~,~,rocAUC] = roc1(30-qq.totRaw(a.i2(i2)),qq.dx(a.i1(i1)),[0:1:30],0); disp([rocAUC, length(a.i1(i1))/1000])
toc
%%

subplot(221); ax=buildUpPlotOpt([qq.dat(a.i1,:),qq.dat(a.i2,:)],qq.dx(a.i1),[0 1],{'MOCA','MOCA'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+MOCA')
buildUpPlotOptSimple(qq.dat(a.i1,:),qq.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq.dat(a.i2,:),qq.dx(a.i2),[0 1],[0 0 0.5]);  



d.dx0 = T1.NACCUDSD; d.dx=d.dx0; d.dx(d.dx0==2)=1;

id = str2double(cellstr(erase(T1.NACCID,'NACC')));
uid = unique(id);


% how common is each Dx, and how common is progression to a new Dx (specific progressions like CN->MCI analyzed)
d.visitDelta_ = [];
for k=1:length(uid),
    i = find(id==uid(k));
    d.visitSpan(k) = max(T1.NACCDAYS(i));
    visitDelta = diff(T1.VISITYR(i)*12+T1.VISITMO(i));
    d.visitDelta_ = [d.visitDelta_;visitDelta];
    % if length(i) > 8, keyboard; end
    % d.visitYear(k) = T1.VISITYR(i);
    d.nVisits(k) = length(i);
    d.nDx(k) = length(unique(d.dx(i)));
    d.n1(k) = (sum(d.dx(i)==1) > 0);
    d.n1a(k) = (sum(d.dx(i)==1) > 1);
    
    d.n13(k) = (sum(d.dx(i)==1) > 0) & (sum(d.dx(i)==3) > 0);
    d.n14(k) = (sum(d.dx(i)==1) > 0) & (sum(d.dx(i)==4) > 0);
    d.n12(k) = (sum(d.dx(i)==1) > 0) & (sum(d.dx(i)==2) > 0);
    d.n34(k) = (sum(d.dx(i)==3) > 0) & (sum(d.dx(i)==4) > 0);
    d.n134(k) = (sum(d.dx(i)==1) > 0) & (sum(d.dx(i)==3) > 0) & (sum(d.dx(i)==4) > 0);
    
    d.n13a(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==3) > 1);
    d.n14a(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==4) > 1);
    d.n12a(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==2) > 1);
    d.n34a(k) = (sum(d.dx(i)==3) > 1) & (sum(d.dx(i)==4) > 1);
    d.n134a(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==3) > 1) & (sum(d.dx(i)==4) > 1);
    
    d.n13b(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==3) > 0);
    d.n14b(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==4) > 0);
    d.n12b(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==2) > 0);
    d.n34b(k) = (sum(d.dx(i)==3) > 1) & (sum(d.dx(i)==4) > 0);
    d.n134b(k) = (sum(d.dx(i)==1) > 1) & (sum(d.dx(i)==3) > 0) & (sum(d.dx(i)==4) > 0);
end
[sum(d.nDx==1),sum(d.nDx==2),sum(d.nDx==3),sum(d.nDx==4)]
[sum(d.n12),sum(d.n13),sum(d.n14),sum(d.n34),sum(d.n134)]  
[sum(d.n12a),sum(d.n13a),sum(d.n14a),sum(d.n34a),sum(d.n134a)]  
[sum(d.n12b),sum(d.n13b),sum(d.n14b),sum(d.n34b),sum(d.n134b)]  

figure; 
d.iMMSE1=(T1.MMSECOMP==1);
d.iMOCA1=(T1.MOCACOMP==1);
d.iMMSE2=(T1.NACCMMSE<=30 & T1.NACCMMSE>=0);
d.iMOCA2=(T1.NACCMOCA<=30 & T1.NACCMOCA>=0);

subplot(511); hist(T1.VISITYR,2001:2024); 
subplot(512); hist(T1.VISITYR(d.iMMSE1),2001:2024); 
subplot(513); hist(T1.VISITYR(d.iMOCA1),2001:2024)
subplot(514); hist(T1.VISITYR(d.iMMSE2),2001:2024); 
subplot(515); hist(T1.VISITYR(d.iMOCA2),2001:2024)
for k=1:5, subplot(5,1,k); ylim([0 12000]); grid; end; shg


dx1=(d.dx-1)/2; i=ismember(dx1,[0 1])&d.iMMSE2; [~,~,rocAUC] = roc1(30-T1.NACCMMSE(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]
dx1=(d.dx-1)/3; i=ismember(dx1,[0 1])&d.iMMSE2; [~,~,rocAUC] = roc1(30-T1.NACCMMSE(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]
dx1=(d.dx-1)/2; i=ismember(dx1,[0 1])&d.iMOCA2; [~,~,rocAUC] = roc1(30-T1.NACCMOCA(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]
dx1=(d.dx-1)/3; i=ismember(dx1,[0 1])&d.iMOCA2; [~,~,rocAUC] = roc1(30-T1.NACCMOCA(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]
dx1=(d.dx-1)/2; i=ismember(dx1,[0 1])&d.iMOCA2; [~,~,rocAUC] = roc1(30-T1.MOCATOTS(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]
dx1=(d.dx-1)/3; i=ismember(dx1,[0 1])&d.iMOCA2; [~,~,rocAUC] = roc1(30-T1.MOCATOTS(i),dx1(i),[0:1:30],0); [rocAUC, sum(i)/1000]

%% single-test build-up plot for MOCA !!
%mmseDx = x.mmseN.dx; 
mocaDx = x.mocaN.dx; 

figure; set(gcf,'position',[0 0 1200 800])
subplot(221);  buildUpPlotOpt(x.mocaN.dat,mocaDx,[0 1],'MOCA','ROC',[20 20],[30 6],[0.5 0.85])
% subplot(222);  buildUpPlotOpt(x.mmse1.dat,mmseDx,[0 1],'MMSE','ROC',[20 20],[30 6],[0.5 0.78]);
subplot(223);  buildUpPlotOpt(x.mocaN.dat,mocaDx,[0 2],'MOCA','ROC',[20 20],[30 6],[0.5 1.02])
% subplot(224);  buildUpPlotOpt(x.mmse1.dat,mmseDx,[0 2],'MMSE','ROC',[20 20],[30 6],[0.5 1.02])

%% double-test build-up plot for MOCA !!

figure; set(gcf,'position',[0 0 1200 800])
% 
% qq=x.mmse1; a = testMatchingDifferentDate (qq,qq,qq.nSameDayMoca==0, qq.nSameDayMoca==0,'less than'); i1=a.i1; i2=a.i2; 
% subplot(222); ax=buildUpPlotOpt([qq.dat(a.i1,:),qq.dat(a.i2,:)],qq.dx(a.i1),[0 1],{'MMSE','MMSE'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
% ylabel(ax(1),'MMSE+MMSE')
% buildUpPlotOptSimple(qq.dat(a.i1,:),qq.dx(a.i1),[0 1],[0 0 1]);  
% buildUpPlotOptSimple(qq.dat(a.i2,:),qq.dx(a.i2),[0 1],[0 0 0.5]);  

tic;
qq=x2.moca1; a = testMatchingDifferentDate (qq,qq,[], [],'less than'); i1=a.i1; i2=a.i2; 
disp([num2str(toc),' <-- MOCA+MOCA'])
subplot(221); ax=buildUpPlotOpt([qq.dat(a.i1,:),qq.dat(a.i2,:)],qq.dx(a.i1),[0 1],{'MOCA','MOCA'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+MOCA')
buildUpPlotOptSimple(qq.dat(a.i1,:),qq.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq.dat(a.i2,:),qq.dx(a.i2),[0 1],[0 0 0.5]);  
tic; 
qq=x2.faq1; qq.dat = (qq.dat==0); a = testMatchingDifferentDate (qq,qq,[], [],'less than'); i1=a.i1; i2=a.i2; 
disp([num2str(toc),' <-- FAQ+FAQ'])

subplot(222); ax=buildUpPlotOpt([qq.dat(a.i1,:),qq.dat(a.i2,:)],qq.dx(a.i1),[0 1],{'FAQ','FAQ'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'FAQ+FAQ')
buildUpPlotOptSimple(qq.dat(a.i1,:),qq.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq.dat(a.i2,:),qq.dx(a.i2),[0 1],[0 0 0.5]);  
toc
tic
qq1=x2.faq1; qq2=x2.moca1; qq1.dat = (qq1.dat==0); a = testMatchingDifferentDate (qq1,qq2,[], [],'less than'); i1=a.i1; i2=a.i2; 
disp([num2str(toc),' <-- FAQ+MOCA'])
subplot(223); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'FAQ','MOCA'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'FAQ+MOCA')
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]);  
toc
tic
qq1=x2.moca1; qq2=x2.faq1; qq2.dat = (qq2.dat==0); a = testMatchingDifferentDate (qq1,qq2,[], [],'less than'); i1=a.i1; i2=a.i2; 
disp([num2str(toc),' <-- MOCA+FAQ'])
subplot(224); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','FAQ'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+FAQ')
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]);  
toc
%%

tic
figure; 
qq1=x2.moca1; qq2=x2.faq1; qq2.dat = (qq2.dat==0); a = testMatchingDifferentDate (qq1,qq2,[], [],'equal'); i1=a.i1; i2=a.i2; 
disp([num2str(toc),' <-- MOCA+FAQ'])
subplot(224); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','FAQ'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+FAQ');
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]); 
toc
%%
figure;
tic
qq1=x2.moca1; qq2=x2.bc1a; qq2.dat = (qq2.dat==0); [a.i,a.i1,a.i2] = intersect(find(qq1.i),find(qq2.i));
disp([num2str(toc),' <-- MOCA+CB1a'])
subplot(223); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','CB1a'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+CB1a');
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]); 
toc
%%
tic
qq1=x2.moca1; qq2=x2.BC1b; qq2.dat = (qq2.dat==0); [a.i,a.i1,a.i2] = intersect(find(qq1.i),find(qq2.i));
disp([num2str(toc),' <-- MOCA+FAQ'])
subplot(223); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','FAQ'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+FAQ');
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]); 
toc

tic
qq1=x2.moca1; qq2=x2.BC2a; qq2.dat = (qq2.dat==0); [a.i,a.i1,a.i2] = intersect(find(qq1.i),find(qq2.i));
disp([num2str(toc),' <-- MOCA+FAQ'])
subplot(223); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','FAQ'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+FAQ');
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]); 
toc

tic
qq1=x2.moca1; qq2=x2.BC2b; qq2.dat = (qq2.dat==0); [a.i,a.i1,a.i2] = intersect(find(qq1.i),find(qq2.i));
disp([num2str(toc),' <-- MOCA+FAQ'])
subplot(223); ax=buildUpPlotOpt([qq1.dat(a.i1,:),qq2.dat(a.i2,:)],qq1.dx(a.i1),[0 1],{'MOCA','BC2B'},'ROC',[20 20],[],[0.5 0.9]);  ylabel(ax(1),'MoCA & eCog Combined')
ylabel(ax(1),'MOCA+FAQ');
buildUpPlotOptSimple(qq1.dat(a.i1,:),qq1.dx(a.i1),[0 1],[0 0 1]);  
buildUpPlotOptSimple(qq2.dat(a.i2,:),qq1.dx(a.i1),[0 1],[0 0 0.5]); 
toc


%%

k1=find(d.n134);
for k2=1:20,
disp(d.dx(find(id==uid(k1(k2))))');
end
%%
% 
% x.mmse1=readItemLevelMMSE(filename); % detailed item-level MMSE data
% disp('mmse'); toc
% disp('moca'); toc

%% Score Change analaysis: Examine the trajectories of test score changes over time for MoCA and CogBat items
%Setup the data for this analysis
dx0 = x2.moca1.dx0;
dx0a = dx0; dx0a(T1.NACCTMCI>2 & dx0==1)=NaN;
x2.cb = organizeCogBatData(x2);
clear q;
i = isfinite(sum(horzcat(x2.cb.MOCA1.dat0{:})+x2.moca1.datSum0+dx0a,2)); length(i), sum(i)
q = loadBasicInfoNACC(T1,i);
CDRSOB = T1.CDRSUM; CDRSOB(CDRSOB==99)=NaN;
CDR = T1.CDRGLOB; CDR(CDR==99)=NaN;
% q.dat0 = [x2.moca1.datSum0, sum(x2.moca1.dat0(:,20:24),2), horzcat(x2.cb.MOCA1.dat0{:})];  % Add Moca-30 & Moca-5(mem)
q.dat0 = [-CDRSOB, sum(x2.moca1.dat0(:,20:24),2), horzcat(x2.cb.MOCA1.dat0{:})];             % Add CDR-SB  & Moca-5(mem)
q.dat = q.dat0(i,:);
q.names = [{'CDR-SOB'},{'MoCA-5'},x2.cb.MOCA1.names];
q.colours = [0 0 0; 1 0 0; colororder('glow12')]; 
q.colorOrder = [1, 2, 2+x2.cb.MOCA1.colorOrder];
%%

iCB2 = [1:4,7:12];
i = isfinite(sum(horzcat(x2.cb.MMSE1.dat0{iCB2})+dx0a,2)); length(i), sum(i)  % remove last two vars (BVRT) because of limited n for that data
q2 = loadBasicInfoNACC(T1,i);
q2.dat0 = [-CDRSOB, sum(x2.mmse1.datSum0,2), horzcat(x2.cb.MMSE1.dat0{iCB2})];             % Add CDR-SB  & Moca-5(mem)
q2.dat0 = [-CDRSOB, -CDR, horzcat(x2.cb.MMSE1.dat0{iCB2})];             % Add CDR-SB  & Moca-5(mem)
q2.dat = q2.dat0(i,:);
q2.names = [{'CDR-SOB'},{'MMSE-30'},x2.cb.MMSE1.names{iCB2}];
q2.colours = [0 0 0; 1 0 0; colororder('glow12')]; 
q2.colorOrder = [1, 2, 2+x2.cb.MOCA1.colorOrder(iCB2)];
%% Run the Score Change analysis for UDSv3 & UDSv2 items (MoCA-associated data) 
tic; clear z
for k=1:size(q.dat,2), z{k} = scoreChangeAnalysis(q,k); disp([k, toc]), end; disp('UDSv3 runtime:'); toc
%%
tic; clear z2
for k=1:size(q2.dat,2), z2{k} = scoreChangeAnalysis(q2,k); disp([k, toc]), end; disp('UDSv2 runtime:'); toc
%% Construct matrix m (z{k1}.m{k,kdx}) to summarize serial visit daat to facilitate further analysis
for k1=1:15, z{k1}.m=[]; for k=2:7, for kdx=1:3, 
            iii=find(z{k1}.n>=k & z{k1}.ndxi==(kdx-1)); 
            for kk=1:length(iii), z{k1}.m{k,kdx}(kk,:)=z{k1}.ny{iii(kk)}(1:k)'; end; 
end; end; end
%% Plot the results of the Score Change analysis for UDSv3 items
scoreChangeAnalysisPlotAUC(q,z)
%%
scoreChangeAnalysisPlotSNRvsVisits_CN_MCI_DEM(q,z)
%%
yy = scoreChangeAnalysisPlotSNRvsVisits_CN_MCI_DEM2(q,z);
%%
scoreChangeAnalysisPlotSNRvsVisits_MCI_DEM(q,z,1)

%% Run the Score Change analysis for UDSv2 items (MMSE-associated data) 
tic; clear z2
for k=1:size(q2.dat,2), z2{k} = scoreChangeAnalysis(q2,k); disp([k, toc]), end; toc
%% Plot the results of the Score Change analysis for UDSv2 items
scoreChangeAnalysisPlotAUC(q2,z2)
scoreChangeAnalysisPlotSNRvsVisits_CN_MCI_DEM(q2,z2)
scoreChangeAnalysisPlotSNRvsVisits_MCI_DEM(q2,z2)

%%  variability analysis for U01 proposal Fig 4
m33 = z{1}.m{3,3};
i7 = find(m33(:,1)>-7);
m33a = m33(i7,:);         % remove the most-severe 10% of AD (n: 453 --> 405)
m330 = m33a - mean(m33a); % difference from population average 
m3301 = m330 - m330(:,1); % change from self-baseline - not currently used
% figure; for k=1:15, plot(m3301((k-1)*10+[1:30],:)','-o'); title(k); pause; end

% plot score change in 1st year vs 2nd
figure; set(gcf,'Position',[655 648 1066 267]);
k=9; 
subplot(161); plot([0:2],m33a((k-1)*10+[1:10],:)','-o'); xlim([-.5 2.5]); xlabel('Visit year'); ylabel('CDR-SB score'); title('Raw scores')
subplot(162); plot([0:2],m330((k-1)*10+[1:10],:)','-o'); xlim([-.5 2.5]); xlabel('Visit year'); ylabel('CDR-SB score (relative to group mean)'); title('Relative to mean decline')

subplot(1,6,[3,4]); qd=(diff(m33a')'); plot(qd(:,1)+(rand(size(qd(:,1)))-0.5)/4,qd(:,2)+(rand(size(qd(:,1)))-0.5)/4,'.')
% hold on; plot([-6,4],[-6,4]); shg
qr = corr(qd(:,1),qd(:,2)); 
title({'Lack of relationship between score progression','during 1st year and 2nd years'})
xlabel('CDR-SB progression from year 0 to year 1')
ylabel('CDR-SB progression from year 1 to year 2')
text(-5,4,['r = ',num2str(qr,2)],'Fontname','Arial','Fontsize',12)
text(1,4,['n = ',num2str(length(qd),3)],'Fontname','Arial','Fontsize',12)
xlim([-7,4.5])
ylim([-9,5.2])
% noise estimation
b = regress(mean(m33a)',[[0:2]',[1 1 1]']); 
decayPattern1 = b(1)*[0:2]';
decayPattern2 = [mean(m33a)-mean(m33a(:,1))]';
% for k=1:length(m330), [~,~,resid1(:,k)]=regress(m33a(k,:)'-decayPattern1,[[1 1 1]']);  end
% for k=1:length(m330), [~,~,resid2(:,k)]=regress(m33a(k,:)',[[0:2]',[1 1 1]']);  end
for k=1:length(m330), [~,~,resid1(:,k)]=regress(m33a(k,:)'-decayPattern1,[[1 1 1]']);  end
for k=1:length(m330), [~,~,resid2(:,k)]=regress(m33a(k,:)',[decayPattern1,[1 1 1]']);  end
% qvar.TotRaw = mean(var(m330'));
qvar.TotRaw1 = (sum(resid1.^2)/2);
qvar.NoiseRaw1 = (sum(resid2.^2));
qvar.TotRaw = mean(sum(resid1.^2)/2);
qvar.NoiseRaw = mean(sum(resid2.^2));
qvar.Noise = qvar.NoiseRaw/qvar.TotRaw;
qvar.SlopeVar = 1-qvar.NoiseRaw/qvar.TotRaw;
qvar
% Bootstrap the noise / slopevar analysis to get CIs
N=length(m33a);
for k=1:1e4, i=randi(N,1,N); qvar.Noise_(k) = mean(qvar.NoiseRaw1(i))/mean(qvar.TotRaw1(i)); end
disp([qvar.Noise, mean(qvar.Noise_), std(qvar.Noise_)])
m = 100*mean([qvar.Noise_',1-qvar.Noise_'])';
se = 100*std([qvar.Noise_',1-qvar.Noise_'])';
% subplot(165); bar(100*[qvar.Noise, qvar.SlopeVar]);  
subplot(165); hold on
h = bar([1 2]', m'); 
h.FaceColor='flat'; h.CData=[0 0 1; 1 0 0]; 
plot([1 2; 1 2],m'+se'.*[1;-1],'-k','Linewidth',2); ylim([0 67])
ylabel('Percentage of overall variability')
title('Variability source analyis')
% xticklabels({'Due to noise','Due to slope variation'})
text(1,2,'Random noise','Color','w','FontName','Arial','FontSize',12,'FontWeight','bold','Rotation',90);
text(2,2,'Slope variation','Color','w','FontName','Arial','FontSize',12,'FontWeight','bold','Rotation',90);

ha = findobj(gcf,'type','axes'); 
ha = flipud(ha);
xShift = 0.3*[0 .05 .15 .25];
for k=1:length(ha),
    ha(k).XRuler.TickLabelGapOffset = -1;
    ha(k).Position = ha(k).Position + 0.05*[0 1 0 -2];
    ha(k).Position = ha(k).Position + [xShift(k) 0 0 0];
    if k==4, xticks([]); box on; end
    if k==3, xticks([]); box on; end
end



%%
corr(diff(m330')')
corr(diff(z{1}.m{3,3}')')
qd=(diff(z{1}.m{3,3}')'); figure; plot(qd(:,1),qd(:,2),'.')
qd=(diff(z{1}.m{3,3}')'); figure; plot(qd(:,1)+(rand(size(qd(:,1)))-0.5)/4,qd(:,2)+randn((size(qd(:,1)))-0.5)/4,'.')
corr(diff(z{1}.m{3,3}')')
corr(diff(z{1}.m{2,3}')')
corr(diff(z{1}.m{3,2}')')
qr=randn(1000,3); qr=diff(qr')'; corr(qr(:,1),qr(:,2))
qr=randn(100000,3); qr=diff(qr')'; corr(qr(:,1),qr(:,2))
mean(std(m330'))
sqrt(mean(var(m330')))
for k=1:length(m330), [~,resid(:,k)]=regress(m330(k,:)',[[0:2]',[1 1 1]']);  end
k
help regress
sqrt(mean(sum(resid.^2)))
sqrt(mean(var(m330')))
(mean(var(m330')))
(mean(sum(resid'.^2)))
hist(z{1}.m{3,3}(1,:),30)
shg
hist(z{1}.m{3,3}(:,1),30); shg
mean(z{1}.m{3,3})
mean(z{1}.m{3,3}(z{1}.m{3,3}(:,1)>-7,:))
mean(z{1}.m{3,3}(z{1}.m{3,3}(:,1)>-4,:))
mean(z{1}.m{3,3}(z{1}.m{3,3}(:,1)>-5,:))
(mean(var(m330')))
(mean(sum(resid'^2)))
(mean(sum(resid.^2)))
(mean(sum(resid.^2)))/(mean(var(m330')))
1-(mean(sum(resid.^2)))/(mean(var(m330')))
%%

%% Time-Efficiency plots: SNR vs administration time (linear & squared)
colour1 = q.colours;
% create SNRdata for SNR plots
for k=1:length(z)
    m33 = z{k}.m{3,3};
    if k==1, i7 = find(m33(:,1)>-7); end
    m33a = m33(i7,:);         % remove the most-severe 10% of AD (n: 453 --> 405)
    b = regress(mean(m33a)',[[0:2]',[1 1 1]']);
    decayPattern1 = b(1)*[0:2]';
    % decayPattern2 = [mean(m33a)-mean(m33a(:,1))]';
    for k1=1:length(m33a), [~,~,resid1(:,k1)]=regress(m33a(k1,:)'-decayPattern1,[[1 1 1]']);  end
    for k1=1:length(m33a), [~,~,resid2(:,k1)]=regress(m33a(k1,:)',[decayPattern1,[1 1 1]']);  end
    qvar.TotRaw = mean(sum(resid1.^2)/2);
    qvar.NoiseRaw = mean(sum(resid2.^2));
    SIGNALdata(k,:) = -[mean(m33a(:,2:3))-mean(m33a(:,1))];
    SIGNALdata2 = -[mean(m33a(:,2:3))-mean(m33a(:,1))];
    SNRdata_(k,:) = [SIGNALdata(k,:)./sqrt(qvar.NoiseRaw), SIGNALdata(k,:)./sqrt(qvar.TotRaw)];
    SNRdata2_(k,:) = [SIGNALdata2./sqrt(qvar.NoiseRaw), SIGNALdata2./sqrt(qvar.TotRaw)];
    SNRdata(k) = SNRdata_(k,2);
end
%%

timeList = [30, 1.5, 1.5, 1.5, 1, 3, 1.5, 1.5, 2.5, 3, 1, 1, 7, 2, 2.5]; % admin durations for UDSv3
figure;  
colororder(colour1);
clear h;
% two  versions here 
% 1st: linear SNR w/ sqrt isoclines
subplot (1,6,[1 2]); hold on;
N=15; g=[0:(N-1)]/(N-1); g = round(2*g)/2; colourIso = [g',0*g',1-g'];
for k=1:N, plot([0:0.01:1]*30*cosd((k/(N+1))^2*90)*2,sqrt([0:0.01:1]*0.7*sind((k/(N+1))^2*90)*2),':','color',colourIso(k,:),'linewidth',1.5); end
for k=1:length(timeList), h(k) = plot(timeList(k),SNRdata(k),'.','markersize',30); end 
xlabel('Administration time (min)'); ylabel('Linear SNR')
ylim([0 0.7]); xlim([0 33])
[~,is] = sort(SNRdata,'descend'); 
legend(h(is),q.names(is));  
title('Several Cognitive Battery items display greater efficiency than the CDR-SB')

% 2nd: squared SNR w/ linear isoclines
subplot (1,6,[4 5]); hold on;
for k=1:N, plot([0 30*cosd(k*90/(N+1))]*2,[0 0.5*sind(k*90/(N+1))]*2,':','color',colourIso(k,:),'linewidth',1.5); end
for k=1:length(timeList), h(k) = plot(timeList(k),SNRdata(k)^2,'.','markersize',30); end 
xlabel('Administration time (min)'); ylabel('Squared SNR')
ylim([0 0.5]); xlim([0 33])
legend(h(is),q.names(is));  
title('Several Cognitive Battery items display greater efficiency than the CDR-SB')

%% Squared-only Time-Efficiency plot w/ 4-way-combo: SNR vs administration time 

timeList = [30, 1.5, 1.5, 1.5, 1, 3, 1.5, 1.5, 2.5, 3, 1, 1, 7, 2, 2.5];
figure;  
colororder(colour1);
clear h;

N=15; g=[0:(N-1)]/(N-1); g = round(2*g)/2; colourIso = [g',0*g',1-g'];
 
% squared SNR w/ linear isoclines
subplot (1,4,[1 3]); hold on; clear h
for k=1:N, plot([0 30*cosd(k*90/(N+1))]*2,[0 0.5*sind(k*90/(N+1))]*2,':','color',colourIso(k,:),'linewidth',1.5); end
for k=1:length(timeList), h(k) = plot(timeList(k),SNRdata(k)^2,'.','markersize',30); end 
xlabel('Administration time (min)'); ylabel('Squared SNR')
ylim([0 0.799]); xlim([0 33])
title({'Combining CogBat items yields','both greater efficiency & greater SNR than the CDR-SB'})

% add combos
kk=5;
i = (z{k}.ndxf==1) & (z{k}.n >= kk); ii=find(i);

itemList = [6,4,15,2]; 
g2=[0.6, 0.4, 0.15, 0.15]'; clear g;
for k=1:length(itemList), g(k,1) = 1/nanmean(z{itemList(k)}.nreg(ii,1,5)); end
nx = (1:kk)'-1;
dat = [];
for k1=1:length(ii),
    for ki=1:length(itemList),
        if ki==1, dat = z{itemList(ki)}.ny{ii(k1)}(1:kk);
        else dat = [dat,z{itemList(ki)}.ny{ii(k1)}(1:kk)];
        end
        [b,bint,res,~,stats] = regress(dat*(g(1:ki).*g2(1:ki)), [nx,nx*0+1]);
        slope(k1,ki) = b(1);
        slopeSD(k1,ki) = diff(bint(1,:))/2/1.96;
    end
end
SNR1 = mean(slope)./mean(slopeSD);
SNR2 = SNR1.^2;
for ki=1:length(itemList),
    totTime(ki) = sum(timeList(itemList(1:ki)));
    if ki > 1,
        plot([totTime(ki-1),totTime(ki)],[SNR2(ki-1),SNR2(ki)],'-','color',colour1(itemList(ki),:),'linewidth',1);
        h(ki+length(timeList)) = plot(totTime(ki),SNR2(ki),'s','color',colour1(itemList(ki),:),'markersize',15,'linewidth',2);
    end
end 
[~,is] = sort(SNRdata,'descend'); legend(h([(length(timeList)+[ki:-1:2]),is]),q.names([itemList(end:-1:2),is])); 


%%
ii = [6 4 2]; clear qz;
for k=1:length(ii), qz.dat(:,k) = sum(q.dat(:,1:k)./nanstd(q.dat(:,1:k)),2); end
qz.idList = q.idList;
qz.id = q.id;
qz.regDate = q.regDate;
qz.dx = q.dx;
qz.names = q.names(ii); for k=1:length(ii), qz.names{k}=['+',qz.names{k}]; end

tic;
for k=1:size(qz.dat,2), zc{k} = scoreChangeAnalysis(qz,k); end; toc
%%
tic;
for k=1:1, zc{k} = scoreChangeAnalysis(qz,k); end; toc
%%

scoreChangeAnalysisPlotIndividualExamples_MCI_DEM(z,1)
%%
scoreChangeAnalysisPlotIndividualExamples_MCI_DEM(z2,1,[13,15,17,21])
%%
scoreChangeAnalysisPlotIndividualExamples_MCI_DEM(z,1,[5:8, 13,15,17,21,30:33])

%% prediction plot!!
clear cc
y0 =[0.9 0.78 0.76];
y0 = sig1(4,is([1 2 8 10]))./noi1(4,is([1 2 8 10]))
y0 = sig1(4,is([1 2 8]))./noi1(4,is([1 2 8]))
y1a = y0*sqrt(18/48)*sqrt(1/4);  % = sqrt(36/144)=6/12=1/2; % 48-->18mo  & 5-->3 measurements (6mo interval)
y1b = y0*sqrt(18/48)*sqrt(18/4);  % = 18/12=3/2;            % 48-->18mo  & 5-->19 measurements (1mo interval)
y1c = y0*sqrt(18/48)*sqrt(384/4);  % = sqrt(36/144=6/12=1/2; & 48-->18mo  & 5-->3 measurements (5d/wk interval)
figure; 
LW = 1.5; MS = 20
cc([1:4],1:3) = 0;
cc(2:3,:) = colours(x2.cb.MOCA1.colorOrder(is([2 8])-2),:);
M=0.1;
h1 = plot([M 1]',[y0; y1a],'.--','markersize',MS); for k=1:3, h1(k).Color=cc(k,:); h1(k).LineWidth=LW; end; hold on; 
h2 = plot([M 1]',[y0; y1b],'.-.','markersize',MS);  for k=1:3, h2(k).Color=cc(k,:); h2(k).LineWidth=LW; end
h3 = plot([M 1]',[y0; y1c],'.-','markersize',MS);  for k=1:3, h3(k).Color=cc(k,:); h3(k).LineWidth=LW; end
h4 = plot([1:800]', y1a.*sqrt([1:800]'),'--');  for k=1:3, h4(k).Color=cc(k,:); h4(k).LineWidth=LW; end
h5 = plot([1:800]', y1b.*sqrt([1:800]'),'-.');  for k=1:3, h5(k).Color=cc(k,:); h5(k).LineWidth=LW; end
h6 = plot([1:800]', y1c.*sqrt([1:800]'),'-');  for k=1:3, h6(k).Color=cc(k,:); h6(k).LineWidth=LW; end
set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([0.05 30000]); 
xticks([0.1, 1, 10, 100, 1000]); xtl = {'0','1','10','100','1000'}; xtl{1}=sprintf('Estimate\\newlinebased on\\newlineannual testing\\newlinein 4-year data'); set(gca,'xticklabel',xtl)
yticks([0.1, 1, 10, 100])
plot(795,14.52,'sk','markersize',15); %'facecolor','k')
text(795,14.52,{'    (n=795, SNR=14.5)','     observed for CDR-SOB','    in the Lecanemab trial'},'VerticalAlignment','middle','HorizontalAlignment','left')
ylim([0.1 400])
legend(h6, {'CDR-SOB','Trails B','Benson delayed reproduction'})
% break_axis('axis', 'x','position', 1,'length', 0.2);

%% Let's try combining items!!!
kk=5;
k=6; % Trails b
i = (z{k}.ndxf==1) & (z{k}.n >= kk); ii=find(i);
sig1 = -squeeze(nanmean(z{k}.nreg(i,1,5)));
noi1 = squeeze(nanmean(z{k}.nregSD(i,1,5)));
disp([k, sig1, noi1, sig1/noi1])
k=4; % CatFlu (veg)
sig1 = -squeeze(nanmean(z{k}.nreg(i,1,5)));
noi1 = squeeze(nanmean(z{k}.nregSD(i,1,5)));
disp([k, sig1, noi1, sig1/noi1])

k1_=1; k2_=6;
g1 = 1/nanmean(z{k1_}.nreg(ii,1,5));
g2 = 1/nanmean(z{k2_}.nreg(ii,1,5));
nx = (1:kk)'-1;
alphaList = [0:0.1:1];
for ka=1:length(alphaList), 
    alpha=alphaList(ka);
    for k1=1:length(ii),
        [b,bint,res,~,stats] = regress(alpha*g1*z{k1_}.ny{ii(k1)}(1:kk) + (1-alpha)*g2*z{k2_}.ny{ii(k1)}(1:kk), [nx,nx*0+1]);
        slope(k1,ka) = b(1);
        slopeSD(k1,ka) = diff(bint(1,:))/2/1.96;
    end
end

figure; plot(alphaList,mean(slope)./mean(slopeSD)); 
xlabel('alpha'); ylabel('SNR^1'); title(['2 variable combination: ',q.names{k2_},' + ',q.names{k1_}]);
%% combine 4 items!

k1_=4; k2_=6; k3_=2;k4_=15;
g1 = 1/nanmean(z{k1_}.nreg(ii,1,5));
g2 = 1/nanmean(z{k2_}.nreg(ii,1,5));
g3 = 1/nanmean(z{k3_}.nreg(ii,1,5));
g4 = 1/nanmean(z{k4_}.nreg(ii,1,5));
nx = (1:kk)'-1;
alphaList = [0:0.1:1];
for ka=1:length(alphaList), 
    alpha=alphaList(ka);
    for k1=1:length(ii),
        [b,bint,res,~,stats] = regress(alpha*g1*z{k1_}.ny{ii(k1)}(1:kk) + (1-alpha)*g2*z{k2_}.ny{ii(k1)}(1:kk)+ 0.15*g3*z{k3_}.ny{ii(k1)}(1:kk)+ 0.15*g4*z{k4_}.ny{ii(k1)}(1:kk), [nx,nx*0+1]);
        slope(k1,ka) = b(1);
        slopeSD(k1,ka) = diff(bint(1,:))/2/1.96;
    end
end

figure; plot(alphaList,mean(slope)./mean(slopeSD)); 
xlabel('alpha'); ylabel('SNR^1'); title(['2 variable combination: ',q.names{k2_},' + ',q.names{k1_},' + ',q.names{k3_},' + ',q.names{k4_}]);



%% ---------------------- function list!!!! -----------------------
function z=scoreChangeAnalysis(y,kk)
% needed fields: id, idList, regDate, dx, dat
if nargin<2, kk=1; end
for k=1:length(y.idList); 
    i = y.idList(k)==y.id;
    ii = find(i);
    z.n(k)=sum(i);
    [datemin,imin] = min(y.regDate(i));
    [datemax,imax] = max(y.regDate(i));
    z.nd{k}=years(y.regDate(i)-datemin)';
    z.ndx{k} = y.dx(i);
    z.ndxi(k) = y.dx(ii(imin));
    z.ndxf(k) = y.dx(ii(imax));
    z.ny{k} = y.dat(ii,kk);
    for k1=2:6,
        if (z.n(k) >= k1), 
            ny = y.dat(ii(1:k1),kk);
            nx = (1:k1)'-1;
            [b,bint,res,~,stats] = regress(ny,[nx,ones(k1,1)]); % if prod(b)==0, b, keyboard; end  %check for error
            %keyboard
            z.nreg(k,:,k1-1) = b;                          % coeffs for linear fit (slope & offset)
            z.nregSD(k,:,k1-1) = diff(bint')'/1.96/2;                           % coeffs for linear fit (slope & offset)
            z.nregres(k,k1-1) = std(res);                   % sigma for the residuals of the linear fit
            z.nregres2(k,k1-1) = std(res)/sqrt(1-stats(1)); % estimate of sigma for the data
            z.nregcorr(k,k1-1) = corr(nx,ny);
            z.nregdxf(k,k1-1) = y.dx(ii(k1)); 
            z.nregdxi(k,k1-1) = y.dx(ii(1)); 
            z.nregmean(k,k1-1) = mean(ny); 
            z.nregf(k,k1-1) = ny(k1); 
            z.nregfest(k,k1-1) = mean(ny)+((k1-1)/2)*b(1); 
        else 
            z.nreg(k,:,k1-1) = NaN(1,2);
            z.nregres(k,k1-1) = NaN;
            z.nregres2(k,k1-1) = NaN; 
            z.nregcorr(k,k1-1) = NaN;  
            z.nregdxf(k,k1-1) = NaN; 
            z.nregdxi(k,k1-1) = NaN; 
            z.nregmean(k,k1-1) = NaN; 
            z.nregf(k,k1-1) = NaN; 
            z.nregfest(k,k1-1) = NaN; 
        end
    end
end
%%
% [X,Y,T,AUC] = perfcurve(labels,scores,posclass) % matlab's native roc fn

for k1=2:6, 
i=(z.n'>=k1) & (z.nregdxf(:,k1-1)<2); dx=z.nregdxf(i,k1-1);
dat = z.nregf(i,k1-1);    z.roc11(k1-1,1) = 1-roc1(dat, dx); [A(k1-1,2),ci] = auc([dx,dat],0.05,'hanley'); Aci(k1-1,1)=diff(ci)/1.96/2; %disp([k1,sum(i),z.roc11(k1-1,1)]), % single (final) score only
dat = z.nregfest(i,k1-1); z.roc11(k1-1,2) = 1-roc1(dat, dx); [A(k1-1,2),ci] = auc([dx,dat],0.05,'hanley'); Aci(k1-1,2)=diff(ci)/1.96/2; %disp([k1,sum(i),z.roc11(k1-1,2)]), % regression-estimate of final
dat = z.nregmean(i,k1-1); z.roc11(k1-1,3) = 1-roc1(dat, dx); [A(k1-1,3),ci] = auc([dx,dat],0.05,'hanley'); Aci(k1-1,3)=diff(ci)/1.96/2; %disp([k1,sum(i),z.roc11(k1-1,3)])  % mean of scores
dat = z.nreg(i,1,k1-1);   z.roc11(k1-1,4) = 1-roc1(dat, dx); [A(k1-1,4),ci] = auc([dx,dat],0.05,'hanley'); Aci(k1-1,4)=diff(ci)/1.96/2; %disp([k1,sum(i),z.roc11(k1-1,4)])   % slope of scores
%[A,Aci] = auc([t,y],0.05,'hanley')
%Acid = 
end
for k2=1:3,
    for k1=2:6,
        i = (z.ndxf==k2-1) & (z.n >= k1); 
        d=z.nreg(i,:,k1-1);
        dn=sum(isfinite(d));
        SDr = z.nregres(z.ndxf==k2-1,k1-1);
        SDs = z.nregres2(z.ndxf==k2-1,k1-1);
        cc = z.nregcorr(z.ndxf==k2-1,k1-1);
        z.regsummary(k1-1,:,k2) = [dn, nanmean(d), nanstd(d), nanstd(d)./sqrt(dn), nanmean(SDr), nanmedian(SDr)/nanmedian(1-cc.^2), nanmedian(SDs), nanmedian(cc)];
    end
end
%num2str(z.regsummary,2)
end


function scoreChangeAnalysisPlotAUC(q,z)
figure; 
for k=1:size(q.dat,2), subplot(3,6,k); plot(2:6,z{k}.roc11,'.-','markersize',12); title(q.names{k}); end
for k=1:size(q.dat,2), subplot(3,6,k); ylim([0.5 0.94]); xlim([1.5, 6.5]); grid on; end;
subplot(3,6,7); ylabel('ROC AUC')
subplot(3,6,15); xlabel('Number of visits analyzed')
legend({'Current visit only', 'Current estimated by regression','Mean across visits','Decline across visists (slope)'})
subplot(3,6,18); plot(2:6,z{k}.regsummary(:,1,2),'k'); xlim([1.5, 6.5]); ylabel('N_M_C_I'); xlabel('Number of visits analyzed'); title ('Number of Participants')
end


function scoreChangeAnalysisPlotSNRvsVisits_CN_MCI_DEM(q,z)

%% Different plot of the results of the Score Change analysis
figure; set(gcf,'position', [0 0 1400 500])
titleList = {'CN','MCI','DEM'};
subplot(1,6,1); hold on;
for k=1:size(q.dat,2), 
    y=z{k}.roc11(:,4); yf(k)=y(end); 
    h(k)=plot(2:6,y,'.-','markersize',24,'linewidth',2.5,'Color',q.colours(q.colorOrder(k),:)); 
end
ylim([0.45 0.799]); xlim([1.5, 6.5]); grid on; 
ylabel('ROC AUC')
xlabel('Number of visits analyzed'); title('CN vs MCI')
[~,is] = sort(yf,'descend'); legend(h(is),q.names(is))

for k2=1:3,
    subplot(1,6,2+k2); hold on;
    for k=1:size(q.dat,2),
        for k1=2:6,
            i = (z{k}.ndxf==k2-1) & (z{k}.n >= k1); 
            sig(k1-1) = -squeeze(nanmean(z{k}.nreg(i,1,k1-1)));  sig1(k1-1,k)=sig(k1-1); if ~isfinite(sig), [k,k1, sig], keyboard; end
            noi(k1-1) = squeeze(nanmean(z{k}.nregSD(i,1,k1-1))); noi1(k1-1,k)=noi(k1-1); % if any(noi < 1e-6), k, keyboard; end
        end
        y(:,1) = sig./noi; yf(k)=y(end);
        h(k)=plot(2:6,y,'.-','markersize',24,'linewidth',2.5,'Color',q.colours(q.colorOrder(k),:)); 
    end
    ylim([-0.5 1.9]);
    xlim([2.5, 6.5]); grid on; title(titleList{k2})
    if k2==1, ylabel('Mean Signal-to-Noise Ratio (SNR) for individual participant slope estimates'); end
    xlabel('Number of visits analyzed')
    [~,is] = sort(yf,'descend');
    if k2==2, legend(h(is),q.names(is)); end
end
end

function yy = scoreChangeAnalysisPlotSNRvsVisits_CN_MCI_DEM2(q,z)

%% Different plot of the results of the Score Change analysis
figure; set(gcf,'position', [0 0 600 500])
titleList = {'CN','MCI','DEM'};
% subplot(1,6,1); hold on;
% for k=1:size(q.dat,2), 
%     y=z{k}.roc11(:,4); yf(k)=y(end); 
%     h(k)=plot(2:6,y,'.-','markersize',24,'linewidth',2.5,'Color',q.colours(q.colorOrder(k),:)); 
% end
% ylim([0.45 0.799]); xlim([1.5, 6.5]); grid on; 
% ylabel('ROC AUC')
% xlabel('Number of visits analyzed'); title('CN vs MCI')
% [~,is] = sort(yf,'descend'); legend(h(is),q.names(is))

for k2=1:3,
    subplot(1,1,1); hold on;
    for k=1:size(q.dat,2),
        for k1=2:6,
            i = (z{k}.ndxf==k2-1) & (z{k}.n >= k1); 
            sig(k1-1) = -squeeze(nanmean(z{k}.nreg(i,1,k1-1)));  sig1(k1-1,k)=sig(k1-1); if ~isfinite(sig), [k,k1, sig]; end
            noi(k1-1) = squeeze(nanmean(z{k}.nregSD(i,1,k1-1))); noi1(k1-1,k)=noi(k1-1); % if any(noi < 1e-6), k, keyboard; end
        end
        y(:,1) = sig./noi; yf(k)=y(end);
        yy(:,k,k2) = y(:,1);
        if k2==3,
            yy1(:,k) = yy(:,k,2)*0.8 + yy(:,k,3)*0.2;
            h(k)=plot(1:2,yy1(3:4,k),'.-','markersize',24,'linewidth',2.5,'Color',q.colours(q.colorOrder(k),:));
        end
    end
    if k2==3,
        % ylim([-0.5 1.9]);
        % xlim([2.5, 6.5]); 
        grid on; title(titleList{k2})
        ylabel('Mean Signal-to-Noise Ratio (SNR)');
        xlabel('Number of visits analyzed')
        [~,is] = sort(yy1(end,:),'descend');
        % keyboard
        legend(h(is),q.names(is));
        xlim([.5 2.5]);
        xtick ([1 2])
    end
end
end


function scoreChangeAnalysisPlotSNRvsVisits_MCI_DEM(q,z,flagForExtraPlot)
%% version of the above w/o CN (because we care about progression for those w/ disease!)
if ((nargin > 2) & (flagForExtraPlot==1)), nPlotRows=4; else nPlotRows = 0; end
figure; 
titleList = {'MCI','Mild-AD'}; % clear sig sig1 noi noi1 y
for k2=2:3,
    subplot(nPlotRows,5,2*(k2-2)+[1 2]); hold on;
    for k=1:size(q.dat,2),
        for k1=2:5,
            i = (z{k}.ndxf==k2-1) & (z{k}.n >= k1); 
            sig(k1-1) = -squeeze(nanmean(z{k}.nreg(i,1,k1-1)));  sig1(k1-1,k)=sig(k1-1); if ~isfinite(sig), [k,k1, sig], keyboard; end
            CI = 1.96 * z{k}.nregSD(i,1,k1-1);
            noi(k1-1) = sqrt( mean( (CI/tinv(0.975,k1-2)).^2 ));   noi1(k1-1,k) = noi(k1-1); 
            sigmaY(k1-1) = noi(k1-1) * sqrt(k1-1)*std(1:k1);       sigmaY1(k1-1,k) = sigmaY(k1-1); 
            sigmaYY(k1-1) = sqrt(mean( z{k}.nregres(i,k1-1).^2 )); sigmaYY1(k1-1,k) = sigmaYY(k1-1); 

            % sqrt(mean(bi.^2))*sqrt(xxn-1)*std(xx)/tinv(0.975,xxn-2)
            % noi(k1-1) = squeeze(nanmean(z{k}.nregSD(i,1,k1-1))); noi1(k1-1,k)=noi(k1-1); % if any(noi < 1e-6), k, keyboard; end
        end
        y(:,1) = sig./noi; yf(k)=y(end);
        sig1(:,k) = sig; % normalized signal
        noi1(:,k) = noi;            % noise
        h(k)=plot(2:5,y,'.-','markersize',24,'linewidth',2.5,'Color',q.colours(q.colorOrder(k),:));
        yy(:,k,k2-1) = y(:,1); 
    end
    ylim([-0.04 1.49]);
    xlim([2.8, 5.2]); grid on; title(titleList{k2-1})
    if k2==2, ylabel('Mean Signal-to-Noise Ratio (SNR) for individual participant slope estimates'); end
    xlabel('Number of visits analyzed'); xticks(3:5)
    [~,is] = sort(yf,'descend');
    if k2==3, legend(h(is),q.names(is)); end
    if nPlotRows==4,
        subplot(nPlotRows,5,5+2*(k2-2)+[1 2]); hold on;
        colororder(q.colours(q.colorOrder,:));
        plot([2:5]',sig1./sigmaY1,'.-','markersize',24,'linewidth',2.5); ylabel('SNR/year')
        subplot(nPlotRows,5,10+2*(k2-2)+[1 2]); hold on;
        colororder(q.colours(q.colorOrder,:));
        plot([3:5]',sig1(2:end,:)./sigmaYY1(2:end,:),'.-','markersize',24,'linewidth',2.5); ylabel('SNR/year')
        subplot(nPlotRows,5,15+2*(k2-2)+[1 2]); hold on;
        colororder(q.colours(q.colorOrder,:));
        plot([2:5]',noi1./nanmean(noi1),'.-','markersize',24,'linewidth',2.5); ylabel('Noise in slop')
        
        
    end
end
end

function scoreChangeAnalysisPlotIndividualExamples_MCI_DEM(z,k1,iIndiv)
if nargin < 3, iIndiv = [5:8]; end
figure;
for dx=1:2,
    z{k1}.m7 =[];
    z{k1}.m6 =[];
    i = find(z{k1}.n>=7 & z{k1}.ndxi==dx); for k=1:length(i), z{k1}.m7(k,:)=z{k1}.nregf(i(k),:); end
    i = find(z{k1}.n>=5 & z{k1}.ndxi==dx); for k=1:length(i), z{k1}.m6(k,:)=z{k1}.ny{i(k)}(1:5); end

    subplot(1,2,dx); hold on;
    nIndiv = length(iIndiv); offset = ((1:nIndiv)-mean(1:nIndiv))*0.05*0.0;
    colour1 = [0 0 0]; colourMean = [1 1 1]*0.9;
    colours = colororder('glow12');
    lineStyleList = {'-',':','--','-.'}; lineStyleList = lineStyleList(repmat(1:4,1,5));
    lineWidthList = [1,2,1,1];           lineWidthList = lineWidthList(repmat(1:4,1,5));  MS=24;

    dat7 = -z{k1}.m7'; % note dat7 is actually a subset of dat6 - so use dat6 for the mean
    dat6 = -z{k1}.m6';
    h(1) = plot(mean(dat6'),'color',colourMean,'linewidth',4);
    % for k=1:nIndiv, h(k+1) = plot(dat7(:,iIndiv(k))+offset(k),lineStyleList{k},'color',colour1,'linewidth',lineWidthList(k)); end
    % for k=1:nIndiv, h(k+1+nIndiv) = plot(dat7(:,iIndiv(k))+offset(k),'.','color',colour1,'markersize',MS); end
    for k=1:nIndiv, h(k+1) = plot(dat6(:,iIndiv(k))+offset(k),lineStyleList{k},'color',colours(k,:),'linewidth',lineWidthList(k)); end
    for k=1:nIndiv, h(k+1+nIndiv) = plot(dat6(:,iIndiv(k))+offset(k),'.','color',colour1,'markersize',MS); end

    %leg1 = {'Population mean'}; 
    for k=1:nIndiv, leg1{k} = ['Individual ', num2str(k)]; end
    legend(h(1:1+nIndiv),leg1);
    xlim([0.5 5.5]); xticks(1:5); xticklabels(0:4); xlabel('Visit number (year)'); ylabel('CDR-SB Score'); shg
    %ylim(0.3*[-1,1]+[0, max(max(dat7(:,iIndiv)))]); yticks(0:.5:3.5); 
    ylim(0.3*[-1,1]+[0, max(max(dat6(:,iIndiv)))]); yticks(0:.5:3.5); 
    ylim([-0.3 7.8]); yticks(0:7); 
    if dx==1, title({'Gold-standard measurement of disease state trajectory:','Mean trend and randomly-selected example individuals (MCI)'}); end
    if dx==2, title({'Gold-standard measurement of disease state trajectory:','Mean trend and randomly-selected example individuals (AD)'}); end
    text(h(1).XData(3),h(1).YData(3)+0.05,{['Mean Trend (n=', num2str(size(dat6,2)),')']},'VerticalAlignment','middle','HorizontalAlignment','left','Rotation',20);
    text(0.7,6,{['Avg SD in these 4 examples: ', num2str(mean(std(dat7(:,iIndiv))),2)],...
                ['Avg SD in dataset (n=', num2str(size(dat6,2)),'): ', num2str(mean(std(dat6)),2)],...
                ['SD for mean trend in dataset: ', num2str(std(mean(dat6')),2)]...
                });
    % if dx==2, title('Mean trend and example individuals with AD'); end
end
end


function y = loadBasicInfoNACC(T,i)

dx0 = T.NACCUDSD; % read NACC dx info
dx=NaN*dx0; dx(dx0==1)=0; dx(dx0==2)=0.5; dx(dx0==3)=1; dx(dx0==4)=2;  % convert NACC dx coding to ADNI coding (but with dx=0.5 for dx0=2)
y.dx0_ = dx;     % original full NACC dx coding including dx=0.5 for not normal but not (yet) MCI
dx=NaN*dx0; dx(dx0==1)=0; dx(dx0==2)=0; dx(dx0==3)=1; dx(dx0==4)=2;  % convert NACC dx coding to ADNI coding
y.dx0 = dx;
y.dx0a = y.dx0; y.dx0a(T.NACCTMCI>2 & y.dx0==1)=NaN;
y.dx = y.dx0(i);
y.dx_ = y.dx0_(i);
y.dxa = y.dx0a(i);
y.id = str2double(cellstr(erase(T.NACCID(i),'NACC')));
y.idList = unique(y.id);
y.i = i;
y.regDate = datetime(T.VISITYR(i),T.VISITMO(i),T.VISITDAY(i));
end

function dispTestInfo1(y)
fprintf('Dataset size (Raw / Retained / ID / Unique ID''s):  ');
disp([length(y.i), size(y.dat), length(y.id), length(y.idList)]);
fprintf('Dataset size by dx (CN / MCI / DEM / non-MCI):  ');
disp([sum(y.dx==0),sum(y.dx==1),sum(y.dx==2),sum(y.dx_==0.5)]);
end


function y = loadVarData (y, cb1, cb2)
y.dat0=[];
cb1.dat0(:,strcmp(cb1.vars,'TRAILA')) = 150 - cb1.dat0(:,strcmp(cb1.vars,'TRAILA')); % flip the Trails A & B data
cb1.dat0(:,strcmp(cb1.vars,'TRAILB')) = 300 - cb1.dat0(:,strcmp(cb1.vars,'TRAILB'));
for k=1:length(y.vars),
    vk = y.vars{k};
    if iscell(vk), for k1=1:length(vk), y.dat0{k}(:,k1)= [cb1.dat0(:,strcmp(cb1.vars,vk{k1})), cb2.dat0(:,strcmp(cb2.vars,vk{k1}))]; end
    else y.dat0{k} = [cb1.dat0(:,strcmp(cb1.vars,vk)), cb2.dat0(:,strcmp(cb2.vars,vk))];
    end 
end
y.dx0 = cb1.dx0;
y.dx0a = cb1.dx0a;
y.dx0_ = cb1.dx0_;
end


function q = addColorOrderInfo (q, q1)
for k=1:length(q.vars), 
    q2=q.vars{k}; 
    if iscell(q2), q2=q2{2}; end; 
    size(ii)
    ii=strcmp(q1.vars,q2); disp([k,size(ii),size(q1.colorOrder(ii))])
    q.colorOrder(k)=q1.colorOrder(ii); 
end
end


function z = rocMaximize(x,dx,w0)
dx(~isfinite(sum(x,2)) | dx>1) = NaN;
i = isfinite(dx);
dx = dx(i); x = x(i,:);
N2 = size(x,2); 
rocWeighted = @(w) (1-roc1(x*[1,w]',dx));

if nargin < 3, w0 = 0.01*randn(1,N2-1); end
disp(1 - rocWeighted(w0))
% keyboard
% options = optimoptions(@lsqnonlin,'OptimalityTolerance',1e-6,'Display','iter','FiniteDifferenceStepSize',1e-3);
options = optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',3e-3);
w = lsqnonlin(rocWeighted,w0,[],[],options);
% w = fsolve(rocWeighted,w0);
z.w = w;
z.auc = 1 - rocWeighted(w);
end


function z = rocMaximizeLineSearchA(x,dx,w0,w00)
dx(~isfinite(sum(x,2)) | dx>1) = NaN;
i = isfinite(dx);
dx = dx(i); x = x(i,:); xSD=std(x); x=x./xSD;  % trim out invaid data from x & dx, and self-normalize x
N2 = size(x,2); 
if N2 == 1, x = [x, 0*x]; N2=2; end
rocWeighted = @(w) (1-roc1(x*[1,w]',dx));
if nargin < 4, w00 = 0.1; end
if nargin < 3, w0 = w00*ones(1,N2-1); end
if nargin >= 3, w0 = [w0, w00*ones(1,N2-1-length(w0))]; end
if length(w0) > N2-1, w0 = w0(1:N2-1); end

% dw = [0:.02:1];
Ndw = 25; dwd = 1.26; dw=dwd.^[0:Ndw-1]/dwd^((Ndw-1)/2); K = dw(Ndw-1)/dw(1); %span of coeffs to be tested
[dw1, w1] = deal(ones(size(w0)));
repeat = 0;
% disp([N2, size(x), size(w0), size(w1), size(dw1)])
for k2=1:9,
    for k=1:length(w0),
        done=0;
        while ~done,
            done = 1;
            for k1=1:length(dw), dw1(k)=dw(k1); r(k1,k)=rocWeighted(w0.*w1.*dw1); end
            [z.maxk(k2,k),i]=max(r(:,k)); dw1(k)=dw(i);
            if i==1 & w0(k)*w1(k)>0.001,   w1(k)=w1(k)/K; done=0; end
            if i==Ndw, w1(k)=w1(k)*K; done=0; end
        end
    end
    if k2 >= 2, if z.maxk(k2,k) == z.maxk(k2-1,k), break; end; end;
end
z.nRounds = k2;
z.w = w0.*w1.*dw1;
z.w1 = [1,w0.*w1.*dw1]./xSD;
z.xPre = x;
z.x1 = x*z.w1';
z.aucPre = rocWeighted(0*w0);
z.auc = rocWeighted(w0.*w1.*dw1);
end


function z = rocMaximizeLineSearch(x,dx,w0)
dx(~isfinite(sum(x,2)) | dx>1) = NaN;
i = isfinite(dx);
dx = dx(i); x = x(i,:);
N2 = size(x,2); 
if N2 == 1, x = [x, 0*x]; N2=2; end
rocWeighted = @(w) (1-roc1(x*[1,w]',dx));
if nargin < 3, w0 = 0.01*randn(1,N2-1); end
if nargin >= 3, w0 = [w0,zeros(1,N2-1-length(w0))]; end

dw = [0:.02:1];
dw1=0*w0;
for k2=1:5,
    for k=1:length(w0),
        for k1=1:length(dw), dw1(k)=dw(k1); r(k1,k)=rocWeighted(w0+dw1); end
        [z.maxk(k2,k),i]=max(r(:,k)); dw1(k)=dw(i);
    end
    if k2 >= 2, if z.maxk(k2,k) == z.maxk(k2-1,k), break; end; end;
end
z.nRounds = k2;
z.w = w0+dw1;
z.aucPre = rocWeighted(w0);
z.auc = rocWeighted(w0+dw1);

end


function z = rocMaximizePoly(x,dx,w0,polyOrder)
N2 = size(x,2); 
if nargin < 3, w0 = 0*randn(1,N2-1); end
if nargin < 4, polyOrder = 8; end

dx(~isfinite(sum(x,2)) | dx>1) = NaN;
i = isfinite(dx);
dx = dx(i); x = x(i,:);

% (1) get polyfits: 
[pN, pP] = deal(zeros(N2,polyOrder+1));
interval = deal(zeros(N2,2));
for k=1:N2, 
    z{k} = rocPolyFit(x(:,k),dx,[],polyOrder,1);
    pN(k,:) = z{k}.polyN;
    pP(k,:) = z{k}.polyP;
    interval(k,:) = z{k}.interval;
end
rocWeighted = @(w) (rocPoly([1,w]*pN,[1,w]*pP,[1,w]*interval));
disp(['AUROC (pre):  ',num2str(1-[rocWeighted(w0), roc1(x*[1,w0]',dx)],3)]);
options = optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3);
w = lsqnonlin(rocWeighted,w0,[],[],options);
disp(['AUROC (post):  ',num2str(1-[rocWeighted(w), roc1(x*[1,w]',dx)],3)]);
keyboard
% (2) write anon function to call rocPoly with mixed polys
%

end



function z = rocPolyFit(x,y, ux, fitOrder,plotFlag)
y(~isfinite(x) | y>1 ) = NaN;
i = isfinite(y);
y = y(i); x = x(i);

if nargin < 3, ux = prctile(x,1) : prctile(x,99); end
if isempty(ux), ux = prctile(x,1) : prctile(x,99); end
if nargin < 4, fitOrder=8; end
if nargin < 5, plotFlag=0; end

y(~isfinite(x) | x < ux(1) | x > ux(end) | y>1 ) = NaN;
i = isfinite(y);
y = y(i); x = x(i);

qy0 = hist(x(y==0),ux); z.polyN = polyfit(ux,qy0,fitOrder);
qy1 = hist(x(y==1),ux); z.polyP = polyfit(ux,qy1,fitOrder);
% z.interval = [min(ux),max(ux)];
z.interval = ux([1,end]);
z.auc = rocPoly(z.polyP,z.polyN,z.interval);
z.auc1 = roc1(-x,y);
z.auc2 = roc1(-x,y,fliplr(-ux));

if plotFlag,
    figure;
    ux1 = min(ux):0.1:max(ux);
    subplot(211);
    plot(ux,qy0,'o',ux,qy1,'o'); hold on;
    plot(ux1,polyval(z.polyN,ux1),'-b',ux1,polyval(z.polyP,ux1),'-r');
    title (['AUC = ', num2str([z.auc, z.auc1, z.auc2],3)])
    subplot(212);
    totP = diff(polyval(polyint(z.polyP),[min(ux),max(ux)])); pP = z.polyP / totP; % normalized to form a PDF
    totN = diff(polyval(polyint(z.polyN),[min(ux),max(ux)])); pN = z.polyN / totN;
    plot(ux,qy0/totN,'o',ux,qy1/totP,'o'); hold on;
    plot(ux1,polyval(pN,ux1),'-b',ux1,polyval(pP,ux1),'-r');
end
end

function auc = rocPoly(polyP,polyN,interval)
% computes the ROC AUC (AUROC) when given polynomial coeffs (polyP & polyN) 
% describing the test-score distributions of positive and negative data
% points on the interval "interval"
% Use it something like this: 
% auc = rocPolyFit(x,y, ux, fitOrder) which calls: 
% auc = rocPoly(polyP,polyN,[min(x),max(x)]);

pP = polyP / diff(polyval(polyint(polyP),interval)); % normalized to form a PDF
pN = polyN / diff(polyval(polyint(polyN),interval));
p = polyint(conv(polyint(pP),pN));
auc_analytic = diff(polyval(p,interval));

i1 = linspace(min(interval),max(interval),1000);
tpr = polyval(polyint(pP),i1); tpr=tpr-tpr(1);
fpr = polyval(polyint(pN),i1);fpr=fpr-fpr(1);
auc = 0.5 * sum( (tpr(1:end-1)+tpr(2:end)) .* diff(fpr)); 

% figure; 
% subplot(121); plot(fpr,tpr); title(num2str([auc,auc_analytic]))
% subplot(122); plot(i1,polyval(polyint(conv(polyint(pP),pN)),i1))

end

% 
% fitting algorithm:  
% (1) fit pos & neg (dx=0 vs 1) polynominals (RBFs?)to each regressor to obtain PCs
% (2) use lsqnonlin to optimize RCs
% (2a) call function with RCs to weight the mixing of PCs, 
% (2b) and send the mixed PCs AND the x-value range to the ROC fn
% (3) in the ROC fn use polyint() to integrate the mixed PCs from A to xx  = polyval(polyint(PCp or PCn,xx))
% (4) but note (!!!) that AUROC = integ(ROC) and integ(y(x)dx) = integ(y(t)x'(t)dt) = polyint(conv(polyint(PCp),PCn))


function plotItemAUC(y)
% calculates the ROC AUC for each continuous variable in y.dat0
% & also for binary versions of each variable (with a range of different thresholds used for making it binary)
N2 = size(y.dat,2);
sList = [10:10:90]'; % sensitivity in percent: used for thresholds in thresholded binary versions
p = prctile(y.dat0,sList);
for k=1:N2,
    qq = y.dat0(:,k);
    ii = (isfinite(qq) & (y.dx0==0 | y.dx0==1));
    n(k)=sum(ii);
    [~,~,r(k)] = roc1(qq(ii),y.dx0(ii),[0:1:max(qq(ii))],0);
    for k1=1:length(sList),
        [~,~,r1(k1,k)] = roc1(qq(ii)<=p(k1,k),y.dx0(ii),[0:1:max(qq(ii))],0);  %
    end
end

figure; 
colororder glow12;
subplot(211); hold on; h = plot(0.5+abs(r1'-0.5),'LineWidth',2); 
legend(cellstr(num2str(sList)));
set(h([1,2]),'linestyle','--'); 
set(h([end-1,end]),'linestyle',':');
subplot(211); plot(0.5+abs(r-0.5),'.-.b','MarkerSize',30);  ylabel('ROC AUC')
subplot(212); plot(n,'.-.b','MarkerSize',30); ylabel('N')

for k=1:2, subplot(2,1,k); xticks(1:N2); end 
if isfield(y,'vars'), for k=1:2, subplot(2,1,k); xticklabels(y.vars); end; end
end

function plotThresholdedAUCs(y, sList, yBaseData, testName)
% calculates the ROC AUC for each continuous variable in y.dat0
% & also for binary versions of each variable (with a range of different thresholds used for making each binary)
if nargin < 4, testName='NONE'; end 
if nargin >= 3, compareToRef = 1; else compareToRef = 0; end 
if nargin < 2, sList = [10:10:90]'; end % sensitivity in percent: used for thresholds in thresholded binary versions
dx0 = y.dx0a;
N2 = size(y.dat0,2);
if ~compareToRef, 
    xMCI = y.dat0(dx0==1,:); 
else 
    yBase = yBaseData.datSum0;
    yBase1 = yBaseData.datSum0 - yBaseData.dat0(:,1);
    yBase17 = yBaseData.datSum0 - yBaseData.dat0(:,17);
    xMCI = y.dat0(dx0==1 & isfinite(yBase),:); 
end

p = prctile(xMCI,sList);
% keyboard
sActual1 = 100*squeeze(sum(xMCI <= shiftdim(p',-1)) ./ sum(xMCI<=Inf))';
sActual2 = 100*squeeze(sum(xMCI <  shiftdim(p',-1)) ./ sum(xMCI<=Inf))';
pChoice = 1 + ( abs(sActual1-sList(:)) < abs(sActual2-sList(:)) );

for k=1:N2,
    q = y.dat0(:,k);
    if compareToRef, qq = [yBase, q]; yB=yBase; else, qq=q; end
    if (compareToRef & (k==3|k==4) & strcmp(testName,'MOCA')), qq = [yBase1, q]; yB=yBase1; end  % remove Trails B from MoCA when testing trails A or B
    if (compareToRef & (k==1|k==2|k==5|k==6) & strcmp(testName,'MOCA')), qq = [yBase17, q]; yB=yBase17; end % remove LF from MoCA when testing LF or CF
    ii = (isfinite(sum(qq,2)) & (dx0==0 | dx0==1)); % CN & MCI
    iMCI = (isfinite(sum(qq,2)) & (dx0==1)); % MCI only
    iCN = (isfinite(sum(qq,2)) & (dx0==0)); % MCI only
    n(k)=sum(ii);
    if ~compareToRef, r(k) = roc1(qq(ii),dx0(ii)); % non-binary ROC
    else z = rocMaximizeLineSearchA(qq(ii,:),dx0(ii)); r(k) = z.auc-z.aucPre; % disp([0,z.aucPre, z.auc, z.auc - z.aucPre])
    end
    for k1=1:length(sList),
        if pChoice(k1,k) == 1,
            sActual(k1,k) = nanmean(q(iMCI)<p(k1,k)); % sensitivity
            spActual(k1,k) = 1-nanmean(q(iCN)<p(k1,k)); % specificity
            if ~compareToRef, r1(k1,k) = roc1(q(ii)<p(k1,k),dx0(ii));  % binarized ROC % [0:1:max(qq(ii))],0
            else z = rocMaximizeLineSearchA([yB(ii)*0.1, q(ii)>=p(k1,k)],dx0(ii)); r1(k1,k) = z.auc - z.aucPre; % disp([1,k,k1,z.aucPre, z.auc - z.aucPre])
            end
        else
            sActual(k1,k) = nanmean(q(iMCI)<=p(k1,k)); % sensitivity
            spActual(k1,k) = 1-mean(q(iCN)<=p(k1,k)); % specificity
            if ~compareToRef, r1(k1,k) = roc1(qq(ii)<=p(k1,k),dx0(ii));  % binarized ROC
            else z = rocMaximizeLineSearchA([yB(ii)*0.1, q(ii)>p(k1,k)],dx0(ii)); r1(k1,k) = z.auc - z.aucPre; % disp([2,k,k1,z.aucPre, z.auc - z.aucPre])
            end
        end
    end
end
if compareToRef, baseAUC = 0; else baseAUC = 0.5; end
if ~compareToRef, r = 0.5+abs(r-0.5); r1 = 0.5+abs(r1-0.5); end

[~,iSort] = sort(r,'descend'); 
r1 = r1(:,iSort); r = r(iSort); y.names = y.names(iSort); s1 = sActual(:,iSort);

colours = colororder('parula'); colours = colours(y.colorOrder(iSort),:);
colororder (colours);
% h1 = plot(s1,r1,'.','MarkerSize',30); hold on; % plot AUCs for binarized variables vs sensitivity
for k=1:N2, h1(k) = plot(s1(:,k),r1(:,k),'.','MarkerSize',30,'Color',colours(k,:)); hold on; end % plot AUCs for binarized variables vs sensitivity
for k=1:N2, h2(k) = plot(s1(:,k),r1(:,k),'-','LineWidth',1,'Color',colours(k,:)); end % plot AUCs for binarized variables vs sensitivity
% plot AUCs for multi-valued variables 
for k=1:N2, [~,i]=max(r1(:,k)); plot([s1(i,k); s1(end,k)+0.2 + 0.1*(rand(1)-0.5)],[r1(i,k); r(k)],'.:','LineWidth',1,'MarkerSize',30,'Color',colours(k,:)); end
% for k=1:N2, plot(s1(end,k)+[0; 0.2],[r1(end,k); r(k)],'.:','LineWidth',1,'MarkerSize',30,'Color',colours(k,:)); end
xlim([0 sList(end)+20+10]/100);
axPos = get(gca,'Position');
hL = legend(h1, y.names,'AutoUpdate','off'); hL.Position = [sum(axPos([1,3]))+0.01, axPos(2), 0.10, axPos(4)];
break_axis('axis', 'x','position', (sList(end)+10)/100,'length', 0.1);

xticks([sList, sList(end)+20]/100); 
for k=1:length(sList), xticklab{k}=num2str(sList(k)/100); end; xticklab{k+1}=sprintf('granular\\newline   data'); xticklabels(xticklab);
h = xlabel({'Sensitivity for thresholded data                        ','with different cutoffs                        '}); h.Position = h.Position.*[1,0.6,1];
if compareToRef, ylabel('AUC Improvement'); else ylabel('AUC'); end

LW=1.5;
ylim1=get(gca,'ylim'); plot(0.05*[1 1],ylim1,':k','LineWidth',LW); plot(0.125*[1 1],ylim1,':k','LineWidth',LW); 
h=text(0.05,ylim1(end)*0.95,'median sensitivity of MMSE items','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
h=text(0.125,ylim1(end)*0.95,'median sensitivity of MOCA items','rotation',90,'VerticalAlignment','top','HorizontalAlignment','right','FontName','Arial','FontSize',11);
% plot(n,'.-.b','MarkerSize',30); ylabel('N'); xticks(1:N2); xticklabels(y.vars); 
end


function plotItemCombinationAUC(y)
% calculates the ROC AUC fior each continuous variable in y.dat0
% & also for binary versions of each variable (with a range of different thresholds used for making it binary)
N2 = size(y.dat0,2);
d = y.dat0 ./ std(y.dat0,'omitnan');
meandiff = nanmean(d(y.dx0==1,:))-nanmean(d(y.dx0==0,:));
sign1 = sign(meandiff); d1 = d.*sign1; d1=d1-min(d1); 
dprime = meandiff.*sign1

for k=1:N2,
    qq = d1(:,k);
    ii = (isfinite(qq) & (y.dx0==0 | y.dx0==1));
    n(k)=sum(ii);
    [~,~,r(k)] = roc1(qq(ii),y.dx0(ii),[0:0.01:1]*max(qq(ii)),0);
    for k1=1:N2,
        qqq = d1(:,k1).*dprime(k1)/dprime(k);
        iii = (isfinite(qqq) & isfinite(qq) & (y.dx0==0 | y.dx0==1));
        n1(k1,k) = sum(iii);
        [~,~,r1(k1,k)] = roc1(qq(iii)+qqq(iii),y.dx0(iii),[0:0.01:1]*max(qq(iii)+qqq(iii)),0);  %
    end
end
n1
sign1
y.vars
figure; 
colororder parula;
%subplot(211); hold on; h = plot(0.5+abs(r1'-0.5),'LineWidth',2); 
subplot(211); hold on; h = plot(r1','LineWidth',2); 
%legend(cellstr(num2str(sList)));
%set(h([1,2]),'linestyle','--'); 
%set(h([end-1,end]),'linestyle',':');
%subplot(211); plot(0.5+abs(r-0.5),'.-.b','MarkerSize',30);  ylabel('ROC AUC')
subplot(211); plot(r,'.-.b','MarkerSize',30);  ylabel('ROC AUC')
subplot(212); plot(n,'.-.b','MarkerSize',30); ylabel('N')

for k=1:2, subplot(2,1,k); xticks(1:N2); end 
if isfield(y,'vars'), for k=1:2, subplot(2,1,k); xticklabels(y.vars); end; end
end


function plotItemCombinationPairsAUC(dat0a,y,paramList)
% calculates the ROC AUC for each continuous variable in y.dat0
% & also for binary versions of each variable (with a range of different thresholds used for making it binary)
%if nargin ==3, dat0 = y.dat0(:,paramList); end

N2 = size(y.dat0,2);
dat0 = [dat0a,y.dat0];
d =  dat0./ std(dat0,'omitnan');
meandiff = nanmean(d(y.dx0==1,:))-nanmean(d(y.dx0==0,:));
sign1 = sign(meandiff); d1 = d.*sign1; d1=d1-min(d1);
dprime = meandiff.*sign1;
qq = d1(:,1);
c = [0:.05:1];
cc = corr(d1(y.dx0==0,1),d1(y.dx0==0,2:end),'rows','pairwise');
for k=1:N2,
    qqq = d1(:,k+1); 
    iii = (isfinite(qqq) & isfinite(qq) & (y.dx0==0 | y.dx0==1));
    n(k) = sum(iii);
    if n(k) >100
        for k1=1:length(c),
            qtot = qq(iii)*(1-c(k1)) + qqq(iii)*c(k1);
            % keyboard
            [~,~,r1(k1,k)] = roc1(qtot,y.dx0(iii),[0:0.01:1]*max(qtot),0);  %
        end
        qtot = dprime(1)*qq(iii) + dprime(k+1)*qqq(iii);
        [~,~,r(k)] = roc1(qtot,y.dx0(iii),[0:0.01:1]*max(qtot),0);
        qtot1 = (dprime(1)-sqrt(cc(1)))*qq(iii) + (dprime(k+1)-sqrt(cc(k)))*qqq(iii);
        [~,~,rr(k)] = roc1(qtot,y.dx0(iii),[0:0.01:1]*max(qtot),0);
    end
end

figure; 
colororder parula;
subplot(111); hold on; h = plot(c,r1','LineWidth',2); title('MyPlot')

%keyboard
for k=1:N2, qvars{k}=[y.vars{k},'(',num2str(round(cc(k)*100)),'-',num2str(n(k)),')']; end
legend(qvars);
%set(h([1,2]),'linestyle','--'); 
%set(h([end-1,end]),'linestyle',':');
plot(dprime(2:end)./(dprime(1)+dprime(2:end)),r,'.k','MarkerSize',30);  ylabel('ROC AUC')
% plot((dprime(2:end)-sqrt(cc))./(dprime(1)+dprime(2:end)-2*sqrt(cc)),r,'.k','MarkerSize',30);  ylabel('ROC AUC')
disp([dprime(2:end)',cc',dprime(2:end)'-sqrt(cc)'])
end




function plotItemCorr(y)
% calculates the cross correlation between items for NC and for MCI for each continuous variable in y.dat0
% & also for binary versions of each variable (with a range of different thresholds used for making it binary)
N2=size(y.dat0,2);
[c1,p(:,:,1)] = corr(y.dat0(y.dx0==0,:),'rows','pairwise'); c1(eye(N2)>0)=0; c(:,:,1)=c1;
[c1,p(:,:,2)] = corr(y.dat0(y.dx0==1,:),'rows','pairwise'); c1(eye(N2)>0)=0; c(:,:,2)=c1;
[c1,p(:,:,3)] = corr(y.dat0((y.dx0==0 | y.dx0==1),:),'rows','pairwise'); c1(eye(N2)>0)=0; c(:,:,3)=c1;
figure;
for k=1:3, subplot(2,3,k); imagesc(c(:,:,k),[-1,1]); colorbar; end
for k=1:3, subplot(2,3,k+3); imagesc(-log10(p(:,:,k)),[-10 20]); colorbar; end
if isfield(y,'vars'), 
    for k=1:6, subplot(2,3,k); xticks(1:size(c,1)); yticks(1:size(c,1)); end 
    for k=1:6, subplot(2,3,k); xticklabels(y.vars); yticklabels(y.vars); end
end
end

function practiceEffectAnalysis(q)
testRange = 0:5;
xhist = 0:30;
dx0 = q.dx0;
for k=1:length(testRange),
    for k1=1:2
        i = (q.nPrev0==testRange(k) & dx0==k1-1);
        y = q.datSum0(i);
        n(k,k1) = sum(i);
        m(k,k1) = mean(y);
        sd(k,k1) = std(y);
        yhist(:,k,k1) = hist(y,xhist);
        if k1==1, i0=i; end
        if k1==2, i1=i; end
        % ii{k,k1}=i;
    end;
    rocq(k) = 1-roc1(q.datSum0(i0|i1),dx0(i0|i1));
end
se=sd./sqrt(n);
figure;
subplot(311);
plot(xhist,yhist(:,:,1)./n(:,1)','linewidth',0.5); hold on
plot(xhist,yhist(:,:,2)./n(:,2)','linewidth',2); 
ylabel('Score frequency'); xlabel('Score')

subplot(312); hold on;
delta = 0.02; xval = testRange'+delta*[1,-1];
colour = 'br';
for k1=1:2,
    errorbar(xval(:,k1),m(:,k1),1.96*se(:,k1),['.',colour(k1)],'linewidth',2,'markersize',12); 
    errorbar(xval(:,k1),m(:,k1),sd(:,k1),['.-',colour(k1)],'linewidth',0.5);
end
ylabel('Score'); xlabel('Previous tests administered'); xlim([min(testRange)-0.5, max(testRange)+0.5])

subplot(313);
plot(testRange',rocq,'.-','markersize',24); ylim([0.5 1])
ylabel('ROC AUC'); xlabel('Previous tests administered'); xlim([min(testRange)-0.5, max(testRange)+0.5])
for k=1:3, subplot(3,1,k); grid on; end
end

function y = readItemLevelImaging(T)
% note that these variables are all 0 vs 1 (neg vs pos), with unavailable = -4 or +8
% narrow at some future time to remove other diagnises: PARK PSP PSPIF MSA MSAIF etc etc (and how many participants have each?)
y.vars = {'AMYLPET','AMYLCSF','FDGAD','HIPPATR','TAUPETAD','CSFTAU'};  
% {'FDGFTLD','TPETFTLD','MRFTLD','DATSCAN','',''} %(1-3) for FTD, (4) for lewy body disease
% {'OTHBIOM','OTHBIOMX'} %finding (1) and name(2) of other biomerker (for AD?)
for k=1:length(y.vars),
    qq=T.(y.vars{k}); qq(qq<0 | qq>=8) = NaN;
    y.dat0(:,k) = qq;
end
end

function y = readItemLevelCSF(T)
% note that these variables are all 0 vs 1 (neg vs pos), with unavailable = -4 or +8
% narrow at some future time to remove other diagnises: PARK PSP PSPIF MSA MSAIF etc etc (and how many participants have each?
y.vars = {'CSFABETA','CSFPTAU','CSFTTAU'};

y.adc = T.NACCADC;
y.id = str2double(cellstr(erase(T.NACCID,'NACC')));
y.idList = unique(y.id);  % y.i = i;
y.regDate = datetime(T.CSFLPYR,T.CSFLPMO,T.CSFLPDY);
for k=1:length(y.vars),
    qq=T.(y.vars{k}); % qq(qq<0 | qq>=8) = NaN;
    y.dat0(:,k) = qq;
end
end


function y = readItemLevelMOCA(T)

i = find(string(T.Properties.VariableNames) == "MOCATRAI"); % MOCATRAI is one of the MOCA items
dat0 = T{:,i+[0:21]}; dat0(dat0<0 | dat0>=88) = NaN;
iPrimary = setdiff(1:22,[7,15,16]);
i = isfinite(sum(dat0(:,iPrimary),2)); 

y = loadBasicInfoNACC(T,i);
dat = dat0(i,:);
max(dat)
y.dat(:,1:5) = dat(:,1:5);         % Trails B, copy design, draw clock (1EF, 4VS]
y.dat(:,6:8) = dat(:,6)>=[1:3];    % Naming (3 - Language)
y.dat(:,9:10) = dat(:,8)>=[1:2];   % Digit Span Fwd / Bkwd (2 - Attention)
y.dat(:,11) = dat(:,9);             % Vigilance (1 - Attention)
y.dat(:,12:14) = dat(:,10)>=[1:3];  % Serial 7's (3 - Attention)
y.dat(:,15:16) = dat(:,11)>=[1:2];  % Repeat Sentences (2 - Language)
y.dat(:,17) = dat(:,12);            % Letter Fluency (1 - Executive Function)
y.dat(:,18:19) = dat(:,13)>=[1:2];  % Asbraction (2 - Executive Function)
y.dat(:,20:24) = dat(:,14)>=[1:5];  % Delayed Recall (5 - Memory)
y.dat(:,25:30) = dat(:,17:22);      % Orientation (6 - Orientaion)
y.datSum = sum(y.dat,2);

y.datSumMem     = sum(y.dat(:,20:24),2);       % 5 pts
y.datSumOr      = sum(y.dat(:,25:30),2);       % 6 pts
y.datSumVS      = sum(y.dat(:,2:5),2);         % 4 pts
y.datSumEF      = sum(y.dat(:,[1,17:19]),2);   % 4 pts
y.datSumLang    = sum(y.dat(:,[6:8,15:16]),2); % 5 pts
y.datSumAttn    = sum(y.dat(:,[9:14]),2);      % 6 pts
y.datSumNoMem   = sum(y.dat(:,setdiff(1:30,20:24)),2);

y.domainDat  = [y.datSumMem, y.datSumOr, y.datSumVS, y.datSumEF, y.datSumLang, y.datSumAttn];
y.domainNames  = {'5-word recall','Orientation','Visuospatial','Executive Fn','Language','Attention'};

y.totRaw = T.MOCATOTS(i);
y.datMIS = sum(dat(:,14:16),2);
y.datMIS3 = dat(:,14);  % MIS: spontaneous --> 3 pts
y.datMIS2 = dat(:,15); y.datMIS2(~isfinite(y.datMIS2)) = 0; % MIS: category cue --> 2 pts
y.datMIS1 = dat(:,16); y.datMIS1(~isfinite(y.datMIS1)) = 0; % MIS: multiple choice --> 1 pt

ii=find(i);
for k=1:length(ii),
    y.nPrev(k) = sum(y.id==y.id(k) & y.regDate<y.regDate(k));
    y.nPrev2(k) = sum(y.id==y.id(k) & y.regDate<y.regDate(k));
    y.nPrevCraft(k) = sum(y.id==y.id(k) & y.regDate<y.regDate(k));
    
end
L = length(i);
% y.nPrev0 = make0(y.nPrev);
y.nPrev0=NaN(L,1);                          y.nPrev0(i,:)=y.nPrev;   
y.domainDat0=NaN(L,size(y.domainDat,2));    y.domainDat0(i,:)=y.domainDat;     
y.dat0=NaN(L,size(y.dat,2));                y.dat0(i,:)=y.dat;     
y.datSum0=NaN(L,1);                         y.datSum0(i,:)=y.datSum;  

qMIS = [y.datMIS3,y.datMIS2,y.datMIS1].*[3,2,1];
qq=NaN(L,3); qq(i,:)=qMIS; 
y.datMIS0all = qq;
y.datMIS0 = sum(qq,2);

edu = T.EDUC(i);  edu(edu<0 | edu>=99) = NaN; 
y.education0 = edu;
y.totCorrected0 = T.NACCMOCA; % corrected for education-level if poor score


if 1,
    j = sum(y.datSum==T.NACCMOCA(i))
    j = sum(y.datSum==T.MOCATOTS(i))
end
if 1, disp('MOCA:'); dispTestInfo1(y); end
end


function y = readItemLevelMMSE(T)

i = find(string(T.Properties.VariableNames) == "MMSEORDA");
dat0 = T{:,i+[0:2]};
dat0(dat0<0 | dat0>=88) = NaN;
i = isfinite(sum(dat0,2)); 
dat = dat0(i,:);
y = loadBasicInfoNACC(T,i);

y.totRaw = T.NACCMMSE(i); 
y.datSum = T.NACCMMSE(i); 
y.dat(:,1:3) = dat(:,1:3);
y.dat(:,4) = y.totRaw - sum(y.dat,2);

ii=find(i);
for k=1:length(ii),
    y.nPrev(k) = sum(y.id==y.id(k) & y.regDate<y.regDate(k));
end
qq=NaN*i; qq(i)=y.nPrev;               y.nPrev0 = qq;
qq=NaN*i; qq(i,1:size(y.dat,2))=y.dat; y.dat0 = qq;
qq=NaN*i; qq(i)=y.totRaw;              y.totRaw0 = qq;
qq=NaN*i; qq(i)=y.datSum;              y.datSum0 = qq;

edu = T.EDUC(i);  edu(edu<0 | edu>=99) = NaN; 
y.education = edu;

if 1,
    j = sum(isfinite(y.dat(:,4)))
    [sum(y.dx==0),sum(y.dx==1),sum(y.dx==2),sum(y.dx==0.5)]
end
if 1, disp('MMSE:'); dispTestInfo1(y); end
end



function y = readItemLevelFAQ(T)

i = find(string(T.Properties.VariableNames) == "BILLS");
dat0 = T{:,i+[0:9]};
dat0(dat0<0 | dat0>=8) = NaN;
i = isfinite(sum(dat0,2));
y = loadBasicInfoNACC(T,i);
y.dat0 = dat0;
y.dat(:,1:10) = dat0(i,:);
y.totRaw = sum(y.dat,2);
if 1, disp('FAQ:'); dispTestInfo1(y); end
end



function y = readItemLevelCogBat2a(T)

vars = {'CRAFTVRS','CRAFTURS','DIGFORCT','DIGFORSL','DIGBACCT','DIGBACLS',...
    'CRAFTDVR','CRAFTDRE','CRAFTDTI','CRAFTCUE',...
    'MINTTOTS','MINTTOTW','MINTSCNG','MINTSCNC','MINTPCNG','MINTPCNC'};  % all n=16
dat0 = zeros(length(T.(vars{1})), length(vars));
for k=1:length(vars),
    d = T.(vars{k});
    if max(d)>100, 
        d(d<0 | d>=888) = NaN;
    else
        d(d<0 | d>=88) = NaN;
    end
    dat0(:,k)=d;
end
i = isfinite(sum(dat0,2)); 

y = loadBasicInfoNACC(T,i);
y.dat0 = dat0;
y.dat = dat0(i,:);
y.totRaw = sum(y.dat,2); 
y.vars = vars;
if 1, disp('Cogbat2a:'); dispTestInfo1(y); end
end




function y = readItemLevelCogBat2b(T)

vars = {'CRAFTVRS','CRAFTURS','DIGFORSL','DIGBACLS',...
    'CRAFTDVR','CRAFTDRE',...
    'MINTTOTS','MINTTOTW'};  % trimmed to n=8
dat0 = zeros(length(T.(vars{1})), length(vars));
for k=1:length(vars),
    d = T.(vars{k});
    if max(d)>100, 
        d(d<0 | d>=888) = NaN;
    else
        d(d<0 | d>=88) = NaN;
    end
    dat0(:,k)=d;
end
i = isfinite(sum(dat0,2)); 

y = loadBasicInfoNACC(T,i);
y.dat0 = dat0;
y.dat = dat0(i,:);
y.totRaw = sum(y.dat,2); 
y.vars = vars;
if 1, disp('Cogbat2b:'); dispTestInfo1(y); end
end


function y = readItemLevelCogBat1a(T)

vars = {'LOGIMEM','MEMUNITS','MEMTIME','UDSBENTC','UDSBENTD','UDSBENRS','DIGIF','DIGIFLEN','DIGIB','DIGIBLEN',...
    'ANIMALS','VEG','TRAILA','TRAILARR','TRAILALI','TRAILB','TRAILBRR','TRAILBLI',...
    'WAIS','BOSTON','UDSVERFC','UDSVERFN','UDSVERNF','UDSVERLC','UDSVERLR','UDSVERLN','UDSVERTN','UDSVERTE','UDSVERTI'}; % all: n=29
[dat0, dat00] = deal(zeros(length(T.(vars{1})), length(vars)));
for k=1:length(vars),
    d = T.(vars{k}); d00=d;
    if max(d)>100, 
        d(d<0 | d>=888) = NaN;
    else
        d(d<0 | d>=88) = NaN;
    end
    dat0(:,k)=d; dat00(:,k)=d00;
end
i = isfinite(sum(dat0,2)); 

y = loadBasicInfoNACC(T,i);
y.dat0 = dat0;
y.dat = dat0(i,:);
y.totRaw = sum(y.dat,2); 
y.vars = vars;
if 1, disp('Cogbat1a:'); dispTestInfo1(y); end
end


function y = readItemLevelCogBat1b(T)

vars = {'LOGIMEM','MEMUNITS','UDSBENTD','DIGIF','DIGIB',...
    'ANIMALS','VEG','TRAILA','TRAILB',...
    'WAIS','BOSTON','UDSVERFC','UDSVERLC','UDSVERTN'}; % trimmed to n=19
dat0 = zeros(length(T.(vars{1})), length(vars));
for k=1:length(vars),
    d = T.(vars{k});
    if max(d)>100, 
        d(d<0 | d>=888) = NaN;
    else
        d(d<0 | d>=88) = NaN;
    end
    dat0(:,k)=d;
end
i = isfinite(sum(dat0,2)); 

y = loadBasicInfoNACC(T,i);
y.dat0 = dat0;
y.dat = dat0(i,:);
y.totRaw = sum(y.dat,2); 
y.vars = vars;
if 1, disp('Cogbat1b:'); dispTestInfo1(y); end
end

function cb = organizeCogBatData(x2)

cb.MOCA1.vars = {'ANIMALS','VEG','TRAILA','TRAILB','UDSVERFC','UDSVERLC','CRAFTURS','CRAFTDRE','DIGFORSL','DIGBACLS','MINTTOTS','UDSBENTC','UDSBENTD'};
cb.MOCA1.names = {'Cat Fluency (Ani)','Cat Fluency (Veg)','Trails A','Trails B','Letter Fluency (F)','Letter Fluency (L)','Craft story immed','Craft story recall','Digit Span (F)','Digit Span (B)','MINT','Benson fig copy','Benson fig recall',};
% cb.MMSE1.vars = {'ANIMALS','VEG','TRAILA','TRAILB','UDSVERFC','UDSVERLC','LOGIMEM','MEMUNITS','DIGIFLEN','DIGIBLEN','WAIS','BOSTON','UDSBENTC','UDSBENTD'};   
% cb.MMSE1.names = {'Cat Fluency (Ani)','Cat Fluency (Veg)','Trails A','Trails B','Letter Fluency (F)','Letter Fluency (L)','Logimem immed','Logimem recall','Digit Span (F)','Digit Span (B)','WAIS','BNT','Benson fig copy','Benson fig recall',};
cb.MMSE1.vars = {'ANIMALS','VEG','TRAILA','TRAILB','UDSVERFC','UDSVERLC','LOGIMEM','MEMUNITS','DIGIFLEN','DIGIBLEN','WAIS','BOSTON'};   
cb.MMSE1.names = {'Cat Fluency (Ani)','Cat Fluency (Veg)','Trails A','Trails B','Logimem immed','Logimem recall','Digit Span (F)','Digit Span (B)','WAIS','BNT'};
cb.MOCA1 = loadVarData(cb.MOCA1, x2.cb1a, x2.cb2a); 
cb.MMSE1 = loadVarData(cb.MMSE1, x2.cb1a, x2.cb2a);
cb.MOCA1.colorOrder = [11, 5, 12, 3, 7, 8, 1, 9, 6, 6, 10, 2, 4];
cb.MMSE1.colorOrder = [11, 5, 12, 3, 7, 8, 1, 9, 6, 6, 7, 10, 2, 4];

cb.MOCA2.vars = {{'ANIMALS','VEG'},'UDSBENTD','CRAFTURS','CRAFTDRE',{'UDSVERFC','UDSVERLC'},'UDSBENTC',{'TRAILA','TRAILB'},{'DIGFORSL','DIGBACLS'},'MINTTOTS'};
cb.MOCA2.names = {'Cat Fluency (V+A)','Benson recall','Craft story immed','Craft story recall','Letter Fluency (F+L)','Benson copy','Trails (A+B)','Digit Span (F+B)','MINT'};
cb.MOCA2x.vars = {{'ANIMALS','VEG'},'UDSBENTD',{'CRAFTURS','CRAFTDRE'},{'UDSVERFC','UDSVERLC'},'UDSBENTC',{'TRAILA','TRAILB'},{'DIGFORSL','DIGBACLS'},'MINTTOTS'};
cb.MOCA2x.names = {'Cat Fluency (Veg+Ani)','Benson recall','Craft story 21','Letter Fluency (F+L)','Benson copy','Trails (A+B)','Digit Span (Fwd+Bkwd)','MINT'};
cb.MMSE2.vars = {{'ANIMALS','VEG'},'UDSBENTD','LOGIMEM','MEMUNITS','UDSBENTC',{'TRAILA','TRAILB'},{'DIGIFLEN','DIGIBLEN'},'WAIS','BOSTON'};
cb.MMSE2.names = {'Cat Fluency (V+A)','Benson recall','Logimem immed','Logimem recall','Benson copy','Trails (A+B)','Digit Span (F+B)','WAIS','BNT'};
% cb.MMSE2.vars = {{'ANIMALS','VEG'},'UDSBENTD','LOGIMEM','MEMUNITS',{'UDSVERFC','UDSVERLC'},'UDSBENTC',{'TRAILA','TRAILB'},{'DIGIFLEN','DIGIBLEN'},'WAIS','BOSTON'};
% cb.MMSE2.names = {'Cat Fluency (V+A)','Benson recall','Logimem immed','Logimem recall','Letter Fluency (F+L)','Benson copy','Trails (A+B)','Digit Span (F+B)','WAIS','BNT'};
cb.MMSE2a.vars = {{'ANIMALS','VEG'},'LOGIMEM','MEMUNITS',{'TRAILA','TRAILB'},{'DIGIFLEN','DIGIBLEN'},'WAIS','BOSTON'};
cb.MMSE2a.names = {'Cat Fluency (V+A)','Logimem immed','Logimem recall','Trails (A+B)','Digit Span (F+B)','WAIS','BNT'};
% cb.MMSE2a.vars = {{'ANIMALS','VEG'},'LOGIMEM','MEMUNITS',{'UDSVERFC','UDSVERLC'},{'TRAILA','TRAILB'},{'DIGIFLEN','DIGIBLEN'},'WAIS','BOSTON'};
% cb.MMSE2a.names = {'Cat Fluency (V+A)','Logimem immed','Logimem recall','Letter Fluency (F+L)','Trails (A+B)','Digit Span (F+B)','WAIS','BNT'};
cb.MOCA2 = loadVarData(cb.MOCA2, x2.cb1a, x2.cb2a);    
cb.MOCA2x = loadVarData(cb.MOCA2x, x2.cb1a, x2.cb2a);    
cb.MMSE2 = loadVarData(cb.MMSE2, x2.cb1a, x2.cb2a);    
cb.MMSE2a = loadVarData(cb.MMSE2a, x2.cb1a, x2.cb2a); 

% Note that 2s & 2as are shorter versions (with only the better vars)
% Note that 2a & 2as are versions with

cb.MOCA2s1.vars = {'VEG','CRAFTURS','UDSBENTD'}; cb.MOCA2s1.c0 = {1,4,1};
cb.MOCA2s1.names = {'Cat Fluency (Veg)','Craft story immed','Benson recall'};
cb.MOCA2s.vars = {{'ANIMALS','VEG'},'CRAFTURS','CRAFTDRE','UDSBENTD'}; cb.MOCA2s.c0 = {[1, 2],4,1,1};
cb.MOCA2s.names = {'Cat Fluency (V+A)','Craft story immed','Craft story recall','Benson recall'};
cb.MMSE2s.vars = {{'ANIMALS','VEG'},'LOGIMEM','MEMUNITS','UDSBENTD'}; cb.MMSE2s.c0 = {1,1,1,1};
cb.MMSE2s.names = {'Cat Fluency (V+A)','Logimem immed','Logimem recall','Benson recall'};
cb.MMSE2as.vars = {{'ANIMALS','VEG'},'LOGIMEM','MEMUNITS'}; cb.MMSE2as.c0 = cellmat(1,length(cb.MMSE2as.vars),1,1);
cb.MMSE2as.names = {'Cat Fluency (V+A)','Logimem immed','Logimem recall'};
cb.MOCA2s1 = loadVarData(cb.MOCA2s1, x2.cb1a, x2.cb2a);    
cb.MOCA2s = loadVarData(cb.MOCA2s, x2.cb1a, x2.cb2a);    
cb.MMSE2s = loadVarData(cb.MMSE2s, x2.cb1a, x2.cb2a);    
cb.MMSE2as = loadVarData(cb.MMSE2as, x2.cb1a, x2.cb2a);  

% % q = cb.MOCA2; q1 = cb.MOCA1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MOCA2=q;
% % q = cb.MOCA2a; q1 = cb.MOCA1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MOCA2a=q;
% % q = cb.MOCA2s; q1 = cb.MOCA1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MOCA2s=q;
% q = cb.MOCA2s; q1 = cb.MOCA1; for k=1:length(q.vars), q2=q.vars{k}; if iscell(q2), q2=q2{2}; end; ii=strcmp(q1.vars,q2); q.colorOrder(k)=q1.colorOrder(ii); end; cb.MOCA2s=q;
% % q = cb.MMSE2; q1 = cb.MMSE1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MMSE2=q;
% % q = cb.MMSE2a; q1 = cb.MMSE1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MMSE2=q;
% q = cb.MMSE2s; q1 = cb.MMSE1; for k=1:length(q.vars), q2=q.vars{k}; if iscell(q2), q2=q2{2}; end; ii=strcmp(q1.vars,q2); q.colorOrder(k)=q1.colorOrder(ii); end; cb.MMSE2s=q;
% % q = cb.MMSE2as; q1 = cb.MMSE1; for k=1:length(q.vars), ii=strcmp(q1.vars,q.vars{k}); q.colorOrder(k)=q1.colorOder(ii); end; cb.MMSE2=q;
% cb.MOCA2s1 = addColorOrderInfo (cb.MOCA2s1, cb.MOCA1);
% cb.MOCA2s = addColorOrderInfo (cb.MOCA2s, cb.MOCA1);
% cb.MOCA2 = addColorOrderInfo (cb.MOCA2, cb.MOCA1);
% cb.MOCA2x = addColorOrderInfo (cb.MOCA2x, cb.MOCA1);
% cb.MMSE2 = addColorOrderInfo (cb.MMSE2, cb.MMSE1);
% cb.MMSE2a = addColorOrderInfo (cb.MMSE2a, cb.MMSE1);
% cb.MMSE2s = addColorOrderInfo (cb.MMSE2s, cb.MMSE1);
% cb.MMSE2as = addColorOrderInfo (cb.MMSE2as, cb.MMSE1);

end







function roc = educationSplitAnalysis(score,educ,Dx1,DxList,edSplit)
i = ismember(Dx1, DxList);
x = score(i);
ed = educ(i);
Dx1 = Dx1(i); 
Dx = Dx1;
Dx(Dx1==DxList(1))=0;
Dx(Dx1==DxList(2))=1;

iLow = ed < edSplit;
iHigh = ed >= edSplit;
N = max(x);
[~,~,roc.tot]    = roc1(N-x,Dx,[-1:1:N+1],0);
[~,~,roc.lowEd]  = roc1(N-x(iLow),Dx(iLow),[-1:1:N+1],0);
[~,~,roc.highEd] = roc1(N-x(iHigh),Dx(iHigh),[-1:1:N+1],0);
end


function roc = educationBinAnalysis(score,educ,Dx1,DxList,edBins)
i = ismember(Dx1, DxList);
x = score(i);
ed = educ(i);
Dx1 = Dx1(i); 
Dx = Dx1;
Dx(Dx1==DxList(1))=0;
Dx(Dx1==DxList(2))=1;

N = max(x);
[~,~,roc.tot]    = roc1(N-x,Dx,[-1:1:N+1],0);
for k=1:length(edBins)-1,
    iB = (ed>=edBins(k)) & (ed<edBins(k+1));
    N = max(x(iB));
    [~,~,roc.bins(k)]  = roc1(N-x(iB),Dx(iB),[-1:1:N+1],0);
end
end


function ax_=buildUpPlotOpt(x,Dx1,DxList,testName, whatToPlot,textPos, rocList, yAxisLimits)

if nargin<4, testName='Test'; end
if nargin<5, whatToPlot = 'ROC'; end
if nargin<6, textPos = [15 15]; end
if nargin<7, rocList = []; end
if nargin<8, yAxisLimits = []; end
N = size(x,2);

i = ismember(Dx1, DxList);
x = x(i,:);
Dx1 = Dx1(i); 
Dx = Dx1;
Dx(Dx1==DxList(1))=0;
Dx(Dx1==DxList(2))=1;

% x=x<2.5;
[iOptimal,Z] = OptimizeItemOrder(x,Dx);
Z = [0.5, Z]'; rocAUC=Z;
% for k2=1:N+1, x1 = nansum(x(:,iOptimal(1:k2-1))==1,2);  [~,~,rocAUC(k2,1),ssAVG(k2,1)] = roc1(N-x1,Dx,[0:1:N],0); end;

% Make plot
% coloredCategoryPlot(testName, Z, iOptimal, symbols)
coloredCategoryPlot(testName, Z, iOptimal)
groupNameList = {'CN','MCI','AD'};
groupNames = groupNameList(DxList+1);
if strncmp(whatToPlot,'ROC',3), ylab='Area under the ROC curve'; else Z=ssAVG; ylab='Classification Accuracy'; end
ylabel(ylab); xlabel('Items included')
if iscell(testName), testName = strjoin(testName,'+'); end
title(['Optimal item sequence for an abreviated ', testName,' for ',groupNames{1},' (n=',num2str(sum(Dx1==DxList(1))),') vs ',groupNames{2},' (n=',num2str(sum(Dx1==DxList(2))),')'])
shg
plot([0 length(iOptimal)],0.95*Z(end,1)*[1 1],'-k'); h=text(textPos(1),0.95*Z(end,1),' 95% of full-test value','HorizontalAlignment','left','VerticalAlignment','middle','BackgroundColor',[1 1 1]);
plot([0 length(iOptimal)],0.98*Z(end,1)*[1 1],'-k'); text(textPos(2),0.98*Z(end,1),' 98% of full-test value','HorizontalAlignment','left','VerticalAlignment','middle','BackgroundColor',[1 1 1]);
if length(yAxisLimits)>0, ylim(yAxisLimits); end
ax_= gca;
if length(rocList)>0,  % make inset plot showing the ROC curves
    c = loadColors();
    colorListROC = {[0 0 0], c.GREEN};
    widthListROC = [1.5,2];
    P = get(gca,'position');
    P1 = [P(1)+0.25*P(3), P(2)+0.20*P(4), 0.30*P(3), 0.40*P(4)];
    ax = axes('position',P1);
    axis(ax); hold on;
    for k = (1:length(rocList)), 
        x1 = nansum(x(:,iOptimal(1:rocList(k)))==1,2);  [FPR,TPR,AUC,~] = roc1(N-x1,Dx,[0:1:N],0);
        plot(FPR, TPR,'linewidth',widthListROC(k),'color',colorListROC{k}); 
        legendText{k} = [testName,'-',num2str(rocList(k)),' (AUC=',num2str(AUC,2),')'];
        if k==1, FPR1=FPR; TPR1=TPR; end
    end
    i = find((TPR1-FPR1)==max(TPR1-FPR1),1);
    if strcmp(testName,'eCog'),
        ii = max([3,i-6]):3:i;
    else
        ii = max([3,i-2]):i;
    end
    plot(FPR1(ii), TPR1(ii),'.','markersize',45,'color',[0 0 0]); 
    for kk=1:length(ii), text(FPR1(ii(kk)), TPR1(ii(kk)),num2str(N-(ii(kk)-2)), 'Color',[1 1 1],'FontName','Arial Narrow','FontWeight','Bold','HorizontalAlignment','center','VerticalAlignment','middle'); end; 


    xlabel('1-Specificity'); 
    ylabel('Sensitivity');
    title('ROC Curves')
    legend(legendText,'Location','southeast');
    ax_=[ax_,ax];
end % end if
end % end function



function buildUpPlotOptSimple(x,Dx1,DxList,lineColor,lineWidth)
if nargin < 4, lineColor = 0.9*[1 1 1]; end
if nargin < 5, lineWidth = 2; end

N = size(x,2);
i = ismember(Dx1, DxList);
x = x(i,:);
Dx1 = Dx1(i); 
Dx = Dx1;
Dx(Dx1==DxList(1))=0;
Dx(Dx1==DxList(2))=1;
[iOptimal,Z] = OptimizeItemOrder(x,Dx);
Z = [0.5, Z]'; rocAUC=Z;
plot(0:N,Z,'Color',lineColor,'LineWidth',lineWidth)
end



function coloredCategoryPlot(testName, Z, iOptimal, symbols)
if nargin < 4, symbols={'.','pentagram','^'}; end
for k=1:length(symbols), if strcmp('.',symbols{k}), sizes{k}=18; else sizes{k}=12; end; end

N = length(Z)-1;
itemInfo = createItemInfo();
GREY = 0.9*[1 1 1];
if ~iscell(testName), T = {testName}; else T = testName; end;
[colorList, colors, symbolList, symbolc, sizeList, sizec, labelList, typeList] = deal(cell(0));

for k=1:length(T), 
    colors = {colors{:}, itemInfo.(T{k}).colors{:}};
    for k1=1:length(itemInfo.(T{k}).colors),
        symbolc = {symbolc{:}, symbols{k}};
        sizec = {sizec{:}, sizes{k}};
    end
    colorList = {colorList{:}, itemInfo.(T{k}).colorList{:}};
    typeList = {typeList{:}, itemInfo.(T{k}).typeList{:}};
    labelList = {labelList{:}, itemInfo.(T{k}).labels{:}};
    for k1=1:length(itemInfo.(T{k}).colorList),
        symbolList = {symbolList{:}, symbols{k}};
        sizeList = {sizeList{:}, sizes{k}};
    end
end
for k = 1:length(colors), plot(0,0.5,'.','Marker',symbolc{k},'MarkerSize',sizec{k},'Color',colors{k}); hold on; end %FOR LEGEND
plot(0:N,Z,'.-','linewidth',2,'Color',GREY,'MarkerSize',18); hold on;
for k= 1:N, i=iOptimal(k); plot(k,Z(k+1),'Marker',symbolList{i},'Markersize',sizeList{i},'Color',colorList{i},'MarkerFaceColor',colorList{i}); end

set(gca, 'xtick', 1:N, 'xticklabel', labelList(iOptimal));
grid; legend(typeList,'location','southeast','AutoUpdate','off'); %drawnow;
end



function coloredCategoryPlotOld(testName, DxList, Z) %non-cell testname
if strcmp(testName,'eCogMOCA'),
    for k = 1:length(itemInfo.(testName).colors),
        h = plot(0,0.5,'.','color',itemInfo.(testName).colors{k},'MarkerSize',18,'Marker',itemInfo.(testName).shapes{k}); hold on; 
        if h.Marker ~= '.', set(h,'Marker','pentagram','MarkerFaceColor',itemInfo.(testName).colors{k},'MarkerSize',18); end
    end
    h = plot(0:N,Z,'.-','linewidth',2,'Color',GREY,'markersize',18); hold on;
    for k=1:N,
        h = plot(k,Z(k+1),'.','markersize',18,'Color',colorList{iOptimal(k)},'marker',itemInfo.(testName).shapeList{iOptimal(k)});
        if h.Marker ~= '.', set(h,'Marker','pentagram','MarkerFaceColor',colorList{iOptimal(k)},'MarkerSize',18); end
    end
else
    for k = 1:length(itemInfo.(testName).colors), plot(0,0.5,'.','color',itemInfo.(testName).colors{k},'markersize',18); hold on; end
    plot(0:N,Z,'.-','linewidth',2,'Color',GREY,'markersize',18); hold on;
    for k=1:N,plot(k,Z(k+1),'.','markersize',18,'Color',colorList{iOptimal(k)}); end
end
end



function [OptimalSequence, auc] = OptimizeItemOrder(x,Dx,NN) 
currentList = [];
N = size(x,2);
if nargin<3, NN=N; end

for k2=1:NN,
    remainingX = setdiff(1:N,currentList);
    roc = [];
    for k3=1:length(remainingX),
        x1 = nansum(x(:,[currentList,remainingX(k3)])==1,2);
        [~,~,roc(k3)] = roc1(N-x1,Dx,[0:1:N],0);
    end
    [auc(k2),ii] = max(roc);
    currentList = [currentList, remainingX(ii)];
end
OptimalSequence = currentList;
end % end function



function itemInfo = createItemInfo ()
itemInfo.MMSE.numbers = 1:30;
itemInfo.MMSE.descriptions = {'Current Date','Current Year','Current Month','Current Day','Current Season',...
    'Hospital Name','Hospital Floor','Current City','Current County','Current State',...
    'Registering Ball','Registering ','Registering ',...
    'Reverse Spelling 1','Reverse Spelling 2','Reverse Spelling 3','Reverse Spelling 4','Reverse Spelling 5',...
    'Recalling Ball','Recalling Ball','Recalling Ball',...
    'Naming Watch','Naming Pencil',...
    'Repeating Phrase','Following Written Instructions','Write a Sentence','Copy a Design',...
    'Following Verbal Instructions 1','Following Verbal Instructions 2','Following Verbal Instructions 3'...
    };
itemInfo.MMSE.labels = {'Date','Year','Month','Day','Season',...
    'Hospital','Floor','City','County','State',...
    'Register1','Register2','Register3',...
    'Spell1','Spell2','Spell3','Spell4','Spell5',...
    'Recall1','Recall2','Recall3',...
    'WatchName','PencilName',...
    'Repeat Phrase','Written Instructions','Write Sentence','Copy Design',...
    'Verbal Instructions 1','Verbal Instructions 2','Verbal Instructions 3'...
    };
itemInfo.MMSE.type = {'Time','Time','Time','Time','Time',...
    'Place','Place','Place','Place','Place',...
    'Registration','Registration','Registration',...
    'Executive fn','Executive fn','Executive fn','Executive fn','Executive fn',...
    'Recall','Recall','Recall',...
    'Language','Language',...
    'Language','Language','Language','Executive fn'...
    'Language','Language','Language'...
    };
itemInfo.MOCA.labels = {'Trail Making','Copy Design','Draw clock 1','Draw clock 2','Draw clock 3',...
    'Camel Name','Lion Name','Rhino Name'...
    'Digit Span Fwd','Digit Span Bkwd','Vigilance','Serial 7''s 1','Serial 7''s 2-3','Serial 7''s 4-5',...
    'Repeat Phrase 1','Repeat Phrase 2','Verbal Fluency','Abstraction 1','Abstraction 2',...
    'Recall 1','Recall 2','Recall 3','Recall 4','Recall 5',...
    'Date','Month','Year','Day','Hospital','City',...    
    };
itemInfo.MOCA.type = {'Executive fn','Executive fn','Executive fn','Executive fn','Executive fn',...
    'Language','Language','Language',...
    'Executive fn','Executive fn','Executive fn','Executive fn','Executive fn','Executive fn',...
    'Language','Language','Language','Language','Language',...
    'Recall','Recall','Recall','Recall','Recall',...
    'Time','Time','Time','Time','Place','Place',...    
    };
itemInfo.eCog.labels = {'Mem 1','Mem 2','Mem 3','Mem 4','Mem 5','Mem 6','Mem 7','Mem 8',....
    'Lang 1','Lang 2','Lang 3','Lang 4','Lang 5','Lang 6','Lang 7','Lang 8','Lang 9',....
    'VS 1','VS 2','VS 3','VS 4','VS 5','VS 6',....
    'Plan 1','Plan 2','Plan 3','Plan 4','Plan 5',....
    'Org 1','Org 2','Org 3','Org 4','Org 5','Org 6',....
    'DivAtt 1','DivAtt 2','DivAtt 3','DivAtt 4',...    
    };
itemInfo.eCog.type = {'Memory','Memory','Memory','Memory','Memory','Memory','Memory','Memory',...
    'eLanguage','eLanguage','eLanguage','eLanguage','eLanguage','eLanguage','eLanguage','eLanguage','eLanguage',...
    'Visuospatial','Visuospatial','Visuospatial','Visuospatial','Visuospatial','Visuospatial',...
    'Planning','Planning','Planning','Planning','Planning',...
    'Organization','Organization','Organization','Organization','Organization','Organization',...
    'Divided Attn','Divided Attn','Divided Attn','Divided Attn',...
    };
itemInfo.FAQ.labels = {'Bills','Taxes','Shopping','Games','Stove','Mealprep','CurEvents','PayAttn','RemDates','Travel'};
itemInfo.FAQ.type = {'FAQ','FAQ','FAQ','FAQ','FAQ','FAQ','FAQ','FAQ','FAQ','FAQ'};


cogbat1VarsAll = {'LOGIMEM','MEMUNITS','MEMTIME','UDSBENTC','UDSBENTD','UDSBENRS','DIGIF','DIGIFLEN','DIGIB','DIGIBLEN',...
    'ANIMALS','VEG','TRAILA','TRAILARR','TRAILALI','TRAILB','TRAILBRR','TRAILBLI',...
    'WAIS','BOSTON','UDSVERFC','UDSVERFN','UDSVERNF','UDSVERLC','UDSVERLR','UDSVERLN','UDSVERTN','UDSVERTE','UDSVERTI'}; % all: n=29
cogbat1Vars1 = {'LOGIMEM','MEMUNITS','UDSBENTC','UDSBENTD','DIGIFLEN','DIGIBLEN',...
    'ANIMALS','VEG','TRAILA','TRAILARR','TRAILALI','TRAILB','TRAILBRR','TRAILBLI',...
    'WAIS','BOSTON','UDSVERFC','UDSVERLC','UDSVERTN'}; % trimmed to n=19

itemInfo.CB1a.labels = cogbat1VarsAll; for k=1:length(itemInfo.CB1a.labels), itemInfo.CB1a.type{k}='CB1a'; end
itemInfo.CB1b.labels = cogbat1Vars1;   for k=1:length(itemInfo.CB1b.labels), itemInfo.CB1b.type{k}='CB1b'; end

cogbat2VarsAll = {'CRAFTVRS','CRAFTURS','DIGFORCT','DIGFORSL','DIGBACCT','DIGBACLS',...
    'CRAFTDVR','CRAFTDRE','CRAFTDTI','CRAFTCUE',...
    'MINTTOTS','MINTTOTW','MINTSCNG','MINTSCNC','MINTPCNG','MINTPCNC'};  % all n=16
cogbat2Vars1 = {'CRAFTVRS','CRAFTURS','DIGFORSL','DIGBACLS',...
    'CRAFTDVR','CRAFTDRE',...
    'MINTTOTS','MINTTOTW'};  % trimmed to  n=8
itemInfo.CB2a.labels = cogbat1VarsAll; for k=1:length(itemInfo.CB2a.labels), itemInfo.CB2a.type{k}='CB2a'; end
itemInfo.CB2b.labels = cogbat1Vars1;   for k=1:length(itemInfo.CB2b.labels), itemInfo.CB2b.type{k}='CB2b'; end


for k=1:30, itemInfo.BNT.labels{k} = [num2str(k),'.']; end
for k=1:30, itemInfo.BNT.type{k} = '1'; end

itemInfo.MMSE.typeList = {'Recall','Time','Place','Registration','Executive fn','Language'};
itemInfo.MOCA.typeList = {'Recall','Time','Place','Executive fn','Language'};
itemInfo.eCog.typeList = {'Memory','eLanguage','Visuospatial','Planning','Organization','Divided Attn'};
itemInfo.BNT.typeList = {'1'};
itemInfo.FAQ.typeList  = {'FAQ'};
itemInfo.CB1a.typeList  = {'CB1a'};
itemInfo.CB1b.typeList  = {'CB1b'};
itemInfo.CB2a.typeList  = {'CB2a'};
itemInfo.CB2b.typeList  = {'CB2b'};


c = loadColors();
itemInfo.MMSE.colors  = {c.BLUE,     c.RED,    c.ORANGE,    c.LTBLUE,     c.GREEN,    c.PURPLE};
itemInfo.MOCA.colors  = {c.BLUE,     c.RED,    c.ORANGE,    c.GREEN,    c.PURPLE};
itemInfo.eCog.colors  = {c.BLUE,     c.LTPURPLE,  c.ORANGE,   c.RED,      c.GREEN,    c.LTBLUE};
itemInfo.FAQ.colors  = {c.GREEN};
itemInfo.CB1a.colors  = {c.RED};
itemInfo.CB1b.colors  = {c.PURPLE};
itemInfo.CB2a.colors  = {c.ORANGE};
itemInfo.CB2b.colors  = {c.LTPURPLE};

% itemInfo.MOCA.colors  = {c.BLUE,     c.RED,    c.ORANGE,    c.LTPURPLE,   c.GREEN,    c.PURPLE};
itemInfo.BNT.colors  = {c.RED};
% keyboard
itemInfo.eCogMOCA.colors = [itemInfo.eCog.colors, itemInfo.MOCA.colors];
itemInfo.eCogMOCA.labels = [itemInfo.eCog.labels, itemInfo.MOCA.labels];
itemInfo.eCogMOCA.type = [itemInfo.eCog.type, itemInfo.MOCA.type];
itemInfo.eCogMOCA.typeList = [itemInfo.eCog.typeList, itemInfo.MOCA.typeList];
for k=1:length(itemInfo.eCog.colors), itemInfo.eCogMOCA.shapes{k} = '^'; end
for k2=1:length(itemInfo.MOCA.colors), itemInfo.eCogMOCA.shapes{k+k2} = '.'; end

f = fieldnames(itemInfo);
for k = 1:length(f),   % ["MMSE", "MOCA","BNT","eCog","eCogMOCA"],
    testName = f{k};
    %keyboard
    typeList = itemInfo.(testName).typeList;
    for k=1:length(typeList),
        i = find(startsWith(itemInfo.(testName).type,typeList{k}));
        for k1=1:length(i), itemInfo.(testName).colorList{i(k1)} = itemInfo.(testName).colors{k}; end
    end
end
for k=1:length(itemInfo.eCog.colorList), itemInfo.eCogMOCA.shapeList{k} = '^'; end
for k2=1:length(itemInfo.MOCA.colorList), itemInfo.eCogMOCA.shapeList{k+k2} = '.'; end
end % end function


function c = loadColors()
c.RED = [1 0 0];
c.ORANGE = [1 .7 0];
c.LTBLUE = [.5 .7 1];
c.BLUE = [0 0 1];
c.GREEN = [0 0.6 .2];
c.PURPLE = [.6 .1 1];
c.LTPURPLE = [.8 .7 1];
c.YELLOW = [1 0.9 0.05];
c.GREY = 0.8*[1 1 1];
c.BLUE2 = [.07 0.62, 1];
c.RED2 = [.85 .33 .10];
end


function y = testMatchingDifferentDate (qq1,qq2,i01,i02,dateCompStr, timeLimit)
% new function for ADNI paper 1 (1/19/2024)
if nargin < 3, i01=ones(size(qq.dx)); end
if nargin < 4, i02=ones(size(qq.dx)); end
if  nargin < 6, timeLimit = 1.5; end
if  nargin < 5, 
    compType = 0;
else
    if strcmp(lower(dateCompStr),'greater than'), compType=2; end
    if strcmp(lower(dateCompStr),'less than'), compType=1; end
    if strcmp(lower(dateCompStr),'equal'),     compType=0; end
    if strcmp(lower(dateCompStr),'unequal'),   compType=-1; end
    if strcmp(lower(dateCompStr),'any'),       compType=-2; end
end
if isempty(qq2), qq2=qq1; compType = 1; end
if isempty(i01), i01=ones(size(qq1.dx)); end
if isempty(i02), i02=ones(size(qq2.dx)); end

N=size(qq1.dat,1); 
[idt, i1, i2] = deal(zeros(N,1));
c = 0;
for k=1:N,
    if compType==2,              % 'greater than'
        i = find(i01(k) & i02 & (qq1.id(k)==qq2.id) & (qq1.dx(k)==qq2.dx) & (qq1.regDate(k)>qq2.regDate));
    elseif compType==1,          % 'less than'
        i = find(i01(k) & i02 & (qq1.id(k)==qq2.id) & (qq1.dx(k)==qq2.dx) & (qq1.regDate(k)<qq2.regDate));
    elseif compType==0,          % 'equal'
        i = find(i01(k) & i02 & (qq1.id(k)==qq2.id) & (qq1.dx(k)==qq2.dx) & (qq1.regDate(k)==qq2.regDate));
    elseif compType==-1,         % 'unequal'
        i = find(i01(k) & i02 & (qq1.id(k)==qq2.id) & (qq1.dx(k)==qq2.dx) & (qq1.regDate(k)~=qq2.regDate));
    elseif compType==-2,         % 'any'
        i = find(i01(k) & i02 & (qq1.id(k)==qq2.id) & (qq1.dx(k)==qq2.dx));
    end
    ni(k) = length(i);
    if ni(k) > 0,
        [idt(k), iClosest] = min(abs(years(qq1.regDate(k)-qq2.regDate(i))));
        if idt(k) <= timeLimit,
            % i1(k) = 1;
            % i2(i(iClosest)) = 1;
            c = c+1;
            i1(c) = k;
            i2(c) = i(iClosest);
            %if k>10, keyboard
        end
    end
end
[y.i1, y.i2, y.ni, y.idt] = deal(i1(1:c), i2(1:c), ni, idt);
%figure; subplot(211); hist(idt,300); title('Distribution'); subplot(212); plot(sort(idt),1:N); title('Cumulative distribution')
% x1 = sum(qq.dat-1,2); 
% x2 = sum(qq2.dat-1,2); 
% i = (idt<1.2)&ismember(qq.dx,[0 1]); [~,~,rocAUC]  = roc1(x1(i),qq.dx(i),[0:0.5:30],0); disp('MMSE+ecogSP:'); disp([sum(i)/1000, rocAUC])
% i1 = ii(i); keyboard
% [~,~,rocAUC]  = roc1(x2(i1),qq2.dx(i1),[0:.5:120],0); disp([sum(i)/1000, rocAUC])
end
