function scratch(tr,dat,time,x,y,start)
%
% eye trace inspector
%
eval(unpackHeader(dat));

dind = find(tr==dat.c(:,TR))

figure;
subplot(2,1,1)
title(['trial number ' num2str(tr)])
hold on;
plot(dat.t,dat.lfp(dind,:,10),'k');
plot(dat.t,dat.lfp(dind,:,11),'r');
yline(start{dind},'k');
%set(gca,'XLim',dat.c(dind,[T_START T_STOP]));
%set(gca,'XLim',[dat.c(dind,[T_STIMON T_STIMON])+[-100 1000]]);
%set(gca,'XLim',[dat.c(dind,[T_SACRESET T_SACRESET])+[-600 100]]);
set(gca,'XLim',[dat.c(dind,[T_SACWAIT T_SACWAIT])+[-600 100]]);

subplot(2,1,2);
hold on;
plot(time{dind},x{dind},'k');
plot(time{dind},y{dind},'r');
yline(start{dind},'k');
yline(dat.c(dind,T_STIMON),'b');
yline(dat.c(dind,T_SACWAIT),'g');
yline(dat.c(dind,T_SACRESET),'m');
yline(dat.c(dind,T_FIXBREAK),'Color','r','LineWidth',2);
yline(dat.c(dind,T_REWARD),'y');
%set(gca,'XLim',[dat.c(dind,[T_SACRESET T_SACRESET])+[-600 100]]);
set(gca,'XLim',[dat.c(dind,[T_SACWAIT T_SACWAIT])+[-600 100]]);

foo=[]