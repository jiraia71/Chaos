function fig_qzs_kam(dpi,conjunto)
% QZS-ADV: analise KAM no limite conservativo (zeta1 = zeta2 = 0).
% Uso: fig_qzs_kam                 % todas as figuras, 1200 dpi
%      fig_qzs_kam(150)            % previa rapida
%      fig_qzs_kam(1200,'K2')      % so a figura K2
% Compativel com MATLAB R2019a ou posterior; nenhuma toolbox adicional.
% Le dados_kam.mat, exportado do calculo verificado em Python (integrador
% simpletico de Yoshida de 4a ordem; energia conservada em ~1e-6 em 2000
% periodos; omega(E) por dois metodos independentes).
%
% Hamiltoniano: H = p^2/2 + p2^2/(2 mu) + U(X) + mu beta^2 (X-X2)^2/2 - f X cos t,
% U(X) = 1.55 X^2 - 3 sqrt((1.5 eta)^2 + X^2). Modelo reduzido (1.5 GL): absorvedor
% congelado, H1 = p^2/2 + U(X) - f X cos t. Omega = 1.
%
% K1 esqueleto integravel (potencial, ressonancias, torcao domega/dE)
% K2 quebra dos toros: secoes de Poincare estroboscopicas
% K3 mapas FLI sobre as condicoes iniciais
% K4 fracao regular do espaco de fase em funcao de f
% K5 efeito do acoplamento com o absorvedor (1.5 GL x 2.5 GL)
% K6 destino dos toros com amortecimento (KAM nao se aplica a ele)
% PNG: fundo transparente. PDF: vetorial (mapas e nuvens embutidos como imagem).

if nargin<1, dpi=1200; end
if nargin<2, conjunto='todas'; end
validateattributes(dpi,{'numeric'},{'scalar','integer','positive'});
conjunto=validatestring(conjunto,{'todas','K1','K2','K3','K4','K5','K6'});
root=fileparts(mfilename('fullpath'));
S=load(fullfile(root,'dados_kam.mat'));D=S.D;
names={'monostable','shallow_wells','deep_wells'};
labels={'Monostable','Shallow wells','Deep wells'};
etatex={'\it\eta\rm = \it\eta\rm_{QZS}','\it\eta\rm = 0.60','\it\eta\rm = 0.30'};
for k=1:3, check_case(D.(names{k})); end
out=fullfile(root,sprintf('figuras_matlab_%ddpi',dpi));
if ~exist(out,'dir'), mkdir(out); end
want=@(t) any(strcmp(conjunto,{'todas',t}));
if want('K1'), fig_K1(D,names,labels,etatex,out,dpi); end
if want('K2'), fig_K2(D,names,labels,etatex,out,dpi); end
if want('K3'), fig_K3(D,names,labels,etatex,out,dpi); end
if want('K4'), fig_K4(D,names,labels,etatex,out,dpi); end
if want('K5'), fig_K5(D,out,dpi); end
if want('K6'), fig_K6(D,names,labels,etatex,out,dpi); end
fprintf('Figuras salvas em: %s\n',out);
end

function check_case(C)
assert(numel(C.dE)==numel(C.omega_quad),'omega(E) incompleto.');
assert(size(C.sec_pts,1)==numel(C.sec_f) && size(C.sec_pts,3)==2,'Secoes incompletas.');
assert(size(C.fli_maps,1)==numel(C.fli_f),'Mapas FLI incompletos.');
assert(abs(C.Ubar-C.Umin-C.depth)<1e-12,'Profundidade inconsistente.');
end

% ------------------------------------------------------------------ K1
function fig_K1(D,names,labels,etatex,out,dpi)
W=180;H=142;[fig,L]=new_figure(W,H);
for c=1:3
    C=D.(names{c});xw=C.xw;
    ax=mm_axes(fig,W,H,18+(c-1)*57,82,48,38);hold(ax,'on');
    x=linspace(-xw,xw,1200);u=potential(x,C.eta)-C.Umin;
    plot(ax,x,u,'Color',[.13 .13 .13],'LineWidth',1.1);
    if C.depth>0, plot(ax,[-xw xw],[C.depth C.depth],'--','Color',[.4 .4 .4],'LineWidth',.6); end
    for i=1:numel(C.res_dE)
        e=C.res_dE(i);xx=x;xx(u>e)=NaN;
        plot(ax,xx,e*ones(size(xx)),'Color',qcolor(C.res_q(i)),'LineWidth',.9);
    end
    xlim(ax,[-xw xw]);ylim(ax,[-.03 .95]);style_axes(ax,9);
    xlabel(ax,'\itX');
    if c==1, ylabel(ax,'\itU\rm(\itX\rm) - \itU\rm_{min}'); end
    put_text(L,18+(c-1)*57+24,124,sprintf('%s, %s',labels{c},etatex{c}),11);
    ax=mm_axes(fig,W,H,18+(c-1)*57,26,48,38);hold(ax,'on');
    for over=0:1
        m=C.over==over;
        if sum(m)<3, continue; end
        [dd,o]=sort(C.dE(m));w=C.omega_quad(m);w=w(o);wo=C.omega_ode(m);wo=wo(o);
        col=[.70 .15 .12];if over==1, col=[.11 .25 .48]; end
        plot(ax,dd,w,'Color',col,'LineWidth',1.2);
        plot(ax,dd(1:5:end),wo(1:5:end),'o','MarkerSize',2,'MarkerEdgeColor',[.35 .35 .35],'LineWidth',.4);
    end
    for i=1:numel(C.res_dE)
        plot(ax,[1e-5 1.2],C.res_rho(i)*[1 1],':','Color',qcolor(C.res_q(i)),'LineWidth',.4);
        plot(ax,C.res_dE(i),C.res_rho(i),'o','MarkerSize',3.5,'MarkerFaceColor',qcolor(C.res_q(i)),'MarkerEdgeColor','w','LineWidth',.4);
    end
    if C.depth>0, plot(ax,C.depth*[1 1],[0 1.75],'--','Color',[.4 .4 .4],'LineWidth',.6); end
    set(ax,'XScale','log');xlim(ax,[1e-5 1.2]);ylim(ax,[0 1.75]);style_axes(ax,9);
    xlabel(ax,'\itE\rm - \itU\rm_{min}');
    if c==1, ylabel(ax,'\it\omega\rm(\itE\rm)/\it\Omega'); end
    set(ax,'XTick',[1e-4 1e-2 1],'XTickLabel',{'10^{-4}','10^{-2}','10^{0}'});
    if C.depth>0, t={'d\it\omega\rm/d\itE\rm < 0 inside','d\it\omega\rm/d\itE\rm > 0 above'};ty=.30;
    else, t={'d\it\omega\rm/d\itE\rm > 0 everywhere'};ty=.95; end
    text(ax,.04,ty,t,'Units','normalized','VerticalAlignment','top','FontName','Times New Roman','FontSize',8);
end
leg=mm_axes(fig,W,H,18,4,150,12);set(leg,'XLim',[0 150],'YLim',[0 12],'Visible','off');hold(leg,'on');
items={[.70 .15 .12],'inside the wells';[.11 .25 .48],'over the barrier / single well';qcolor(1),'resonance \it\omega\rm/\it\Omega\rm = \itp\rm/1';
       qcolor(2),'resonance \it\omega\rm/\it\Omega\rm = \itp\rm/2';qcolor(3),'resonance \it\omega\rm/\it\Omega\rm = \itp\rm/3'};
xs=[0 0 0 76 76];ys=[9 5 1 5 1];
for i=1:size(items,1)
    plot(leg,xs(i)+[0 6],ys(i)*[1 1],'Color',items{i,1},'LineWidth',1.2);
    text(leg,xs(i)+8,ys(i),items{i,2},'FontName','Times New Roman','FontSize',8,'VerticalAlignment','middle');
end
plot(leg,76+[3 3],9*[1 1],'o','MarkerSize',2,'MarkerEdgeColor',[.35 .35 .35],'LineWidth',.4);
text(leg,76+8,9,'period by ODE event (independent check)','FontName','Times New Roman','FontSize',8,'VerticalAlignment','middle');
put_text(L,W/2,136,'Integrable skeleton (\itf\rm = 0): potential, resonances and twist d\it\omega\rm/d\itE',12);
export_figure(fig,fullfile(out,'K1_esqueleto_integravel'),dpi);close(fig);
end

% ------------------------------------------------------------------ K2
function fig_K2(D,names,labels,etatex,out,dpi)
W=180;H=168;[fig,L]=new_figure(W,H);sz=35;x0=20;gx=4.5;gy=11;y0=18;
cmap=energy_colormap(512);
for r=1:3
    C=D.(names{r});xw=C.xw;vw=C.vw;
    [XG,VG]=meshgrid(linspace(-xw,xw,400),linspace(-vw,vw,400));
    HG=.5*VG.^2+potential(XG,C.eta)-C.Umin;
    for c=1:4
        ax=mm_axes(fig,W,H,x0+(c-1)*(sz+gx),y0+(3-r)*(sz+gy),sz,sz);hold(ax,'on');
        pts=squeeze(C.sec_pts(c,:,:,:));      % 1500 x 2 x 64
        chaotic=C.sec_fli(c,:)>S_chaos();
        xg=squeeze(pts(:,1,chaotic));yg=squeeze(pts(:,2,chaotic));
        if ~isempty(xg), scatter(ax,double(xg(:)),double(yg(:)),.2,[.55 .55 .55],'filled','MarkerEdgeColor','none'); end
        keep=find(~chaotic);
        if ~isempty(keep)
            xr=squeeze(pts(:,1,keep));yr=squeeze(pts(:,2,keep));
            cc=repmat(min(max(C.sec_E0(keep).'/.8,0),1),size(xr,1),1);
            scatter(ax,double(xr(:)),double(yr(:)),.25,double(cc(:)),'filled','MarkerEdgeColor','none');
            colormap(ax,cmap);caxis(ax,[0 1]);
        end
        if C.depth>0, contour(ax,XG,VG,HG,C.depth*[1 1],'LineColor',[.13 .13 .13],'LineWidth',.45,'LineStyle','--'); end
        if c==1
            for i=1:numel(C.res_dE)
                contour(ax,XG,VG,HG,C.res_dE(i)*[1 1],'LineColor',qcolor(C.res_q(i)),'LineWidth',.35);
            end
        end
        xlim(ax,[-xw xw]);ylim(ax,[-vw vw]);style_axes(ax,8);pbaspect(ax,[1 1 1]);
        text(ax,.03,.97,sprintf('chaotic %.0f%%',100*mean(chaotic)),'Units','normalized','VerticalAlignment','top', ...
            'FontName','Times New Roman','FontSize',7,'BackgroundColor','w','Margin',.5);
        if r==1, put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+3*sz+2*gy+3.5,sprintf('\\itf\\rm = %g',C.sec_f(c)),11); end
        if r==3, xlabel(ax,'\itX'); else, set(ax,'XTickLabel',[]); end
        if c==1, ylabel(ax,'d\itX\rm/d\itt'); else, set(ax,'YTickLabel',[]); end
    end
    put_rot(L,5,y0+(3-r)*(sz+gy)+sz/2,sprintf('%s, %s',labels{r},etatex{r}),10);
end
cb=mm_axes(fig,W,H,x0,8,52,2);manual_colorbar(cb,cmap,[0 .4 .8],'%.1f');
put_left(L,x0+54,9,'initial energy \itE\rm_0 - \itU\rm_{min} (regular orbits)',8);
put_left(L,x0,3.5,'gray: chaotic orbits (FLI > 10);  black dashed: separatrix;  colored curves (\itf\rm = 0.002): resonant tori \itp\rm/\itq',8);
put_text(L,W/2,163,['Stroboscopic Poincar' char(233) ' sections (\itt\rm = 2\pi\itn\rm), conservative limit \it\zeta\rm_1 = \it\zeta\rm_2 = 0, absorber frozen'],12);
export_figure(fig,fullfile(out,'K2_quebra_dos_toros'),dpi);close(fig);
end

% ------------------------------------------------------------------ K3
function fig_K3(D,names,labels,etatex,out,dpi)
W=180;H=168;[fig,L]=new_figure(W,H);sz=35;x0=20;gx=4.5;gy=11;y0=18;cmap=fli_colormap(512);
for r=1:3
    C=D.(names{r});
    for c=1:4
        ax=mm_axes(fig,W,H,x0+(c-1)*(sz+gx),y0+(3-r)*(sz+gy),sz,sz);
        reg=fli_panel(ax,C,squeeze(C.fli_maps(c,:,:)),C.fli_xs,C.fli_vs,cmap);
        text(ax,.03,.97,sprintf('regular %.0f%%',100*reg),'Units','normalized','VerticalAlignment','top', ...
            'FontName','Times New Roman','FontSize',7,'Color','w');
        if r==1, put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+3*sz+2*gy+3.5,sprintf('\\itf\\rm = %g',C.fli_f(c)),11); end
        if r==3, xlabel(ax,'\itX\rm_0'); else, set(ax,'XTickLabel',[]); end
        if c==1, ylabel(ax,'d\itX\rm_0/d\itt'); else, set(ax,'YTickLabel',[]); end
    end
    put_rot(L,5,y0+(3-r)*(sz+gy)+sz/2,sprintf('%s, %s',labels{r},etatex{r}),10);
end
cb=mm_axes(fig,W,H,x0,8,52,2);manual_colorbar(cb,cmap,[0 5 10 15 20]/20,'%g',[0 5 10 15 20]);
put_left(L,x0+54,9,'FLI after 400 periods (dark: KAM tori, bright: chaos)',8);
put_left(L,x0,3.5,'white dashed: separatrix;  white dotted: \itH\rm_0 - \itU\rm_{min} = 0.5 (region used for the regular fraction)',8);
put_text(L,W/2,163,'Fast Lyapunov indicator over initial conditions at \itt\rm = 0 (conservative limit)',12);
export_figure(fig,fullfile(out,'K3_mapas_FLI'),dpi);close(fig);
end

% ------------------------------------------------------------------ K4
function fig_K4(D,names,labels,etatex,out,dpi)
W=120;H=86;[fig,L]=new_figure(W,H);ax=mm_axes(fig,W,H,16,14,98,60);hold(ax,'on');
cols=[.11 .49 .84;.91 .35 .05;.18 .62 .27];
patch(ax,[.01 .30 .30 .01],[0 0 1.02 1.02],[.93 .93 .93],'EdgeColor','none');
h=gobjects(1,4);
for k=1:3
    C=D.(names{k});f=C.frac_f;
    patch(ax,[f fliplr(f)],[C.frac_fli6 fliplr(C.frac_fli10)],cols(k,:),'FaceAlpha',.18,'EdgeColor','none');
    h(k)=plot(ax,f,C.frac_regular,'o-','Color',cols(k,:),'MarkerFaceColor',cols(k,:),'MarkerSize',3.2,'LineWidth',1.2);
end
C=D.shallow_wells;
h(4)=plot(ax,C.abs_f,C.abs_regular,'s--','Color',[.48 .17 .75],'MarkerSize',5,'MarkerFaceColor','w','LineWidth',1);
set(ax,'XScale','log');xlim(ax,[8e-4 .35]);ylim(ax,[0 1.02]);style_axes(ax,11);
xlabel(ax,'forcing amplitude \itf');ylabel(ax,'regular fraction (KAM tori)');
text(ax,.0105,.03,'range of the dissipative maps','Rotation',90,'FontName','Times New Roman','FontSize',7,'Color',[.45 .45 .45]);
lg=legend(h,{sprintf('%s, %s',labels{1},etatex{1}),sprintf('%s, %s',labels{2},etatex{2}), ...
    sprintf('%s, %s',labels{3},etatex{3}),'shallow wells + absorber (2.5 DOF)'},'Location','southwest','Box','off','FontSize',7.5);
set(lg,'ItemTokenSize',[12 6]);
put_text(L,W/2,80,'Initial conditions with \itH\rm_0 - \itU\rm_{min} \leq 0.5; band: FLI threshold from 6 to 10',9);
export_figure(fig,fullfile(out,'K4_fracao_regular'),dpi);close(fig);
end

% ------------------------------------------------------------------ K5
function fig_K5(D,out,dpi)
W=180;H=120;[fig,L]=new_figure(W,H);sz=44;x0=22;gx=7;gy=10;y0=18;cmap=fli_colormap(512);
C=D.shallow_wells;rows={'1.5 DOF (absorber frozen)','2.5 DOF (with absorber)'};
for c=1:3
    f=C.abs_f(c);k=find(abs(C.fli_f-f)<1e-12,1);
    for r=1:2
        ax=mm_axes(fig,W,H,x0+(c-1)*(sz+gx),y0+(2-r)*(sz+gy),sz,sz);
        if r==1, M=squeeze(C.fli_maps(k,:,:));xs=C.fli_xs;vs=C.fli_vs;
        else, M=squeeze(C.abs_maps(c,:,:));xs=C.abs_xs;vs=C.abs_vs; end
        reg=fli_panel(ax,C,M,xs,vs,cmap);
        text(ax,.03,.97,sprintf('regular %.0f%%',100*reg),'Units','normalized','VerticalAlignment','top', ...
            'FontName','Times New Roman','FontSize',8,'Color','w');
        if r==1, put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+2*sz+gy+3.5,sprintf('\\itf\\rm = %g',f),11); end
        if r==2, xlabel(ax,'\itX\rm_0'); else, set(ax,'XTickLabel',[]); end
        if c==1, ylabel(ax,'d\itX\rm_0/d\itt');put_rot(L,5,y0+(2-r)*(sz+gy)+sz/2,rows{r},9);
        else, set(ax,'YTickLabel',[]); end
    end
end
cb=mm_axes(fig,W,H,x0,8,52,2);manual_colorbar(cb,cmap,[0 5 10 15 20]/20,'%g',[0 5 10 15 20]);
put_left(L,x0+54,9,'FLI after 400 periods',8);
put_left(L,x0,3.5,'shallow wells, \it\mu\rm = 0.1, \it\beta\rm = 0.35, absorber at rest relative to the primary mass (\itZ\rm = \itW\rm = 0) at \itt\rm = 0',8);
put_text(L,W/2,115,'Coupling to the absorber shrinks the KAM region',12);
export_figure(fig,fullfile(out,'K5_absorvedor'),dpi);close(fig);
end

% ------------------------------------------------------------------ K6
function fig_K6(D,names,labels,etatex,out,dpi)
W=180;H=118;[fig,L]=new_figure(W,H);sz=35;x0=20;gx=4.5;gy=11;y0=18;cmap=energy_colormap(512);
rows=[2 3];
for r=1:2
    C=D.(names{rows(r)});xw=C.xw;vw=C.vw;
    [XG,VG]=meshgrid(linspace(-xw,xw,400),linspace(-vw,vw,400));
    HG=.5*VG.^2+potential(XG,C.eta)-C.Umin;
    for c=1:4
        ax=mm_axes(fig,W,H,x0+(c-1)*(sz+gx),y0+(2-r)*(sz+gy),sz,sz);hold(ax,'on');
        e=squeeze(C.damp_early(c,:,:,:));l=squeeze(C.damp_late(c,:,:,:));
        xe=squeeze(e(:,1,:));ye=squeeze(e(:,2,:));
        scatter(ax,double(xe(:)),double(ye(:)),.2,[.81 .81 .81],'filled','MarkerEdgeColor','none');
        xl=squeeze(l(:,1,:));yl=squeeze(l(:,2,:));
        cc=repmat(min(max(C.damp_E0.'/.8,0),1),size(xl,1),1);
        s=.25;if C.damp_zeta1(c)>0, s=2; end
        scatter(ax,double(xl(:)),double(yl(:)),s,double(cc(:)),'filled','MarkerEdgeColor','none');
        colormap(ax,cmap);caxis(ax,[0 1]);
        contour(ax,XG,VG,HG,C.depth*[1 1],'LineColor',[.13 .13 .13],'LineWidth',.45,'LineStyle','--');
        xlim(ax,[-xw xw]);ylim(ax,[-vw vw]);style_axes(ax,8);pbaspect(ax,[1 1 1]);
        if r==1
            z=C.damp_zeta1(c);
            if z==0, t='\it\zeta\rm_1 = 0'; else, t=sprintf('\\it\\zeta\\rm_1 = 10^{%d}',round(log10(z))); end
            put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+2*sz+gy+3.5,t,11);
        end
        if r==2, xlabel(ax,'\itX'); else, set(ax,'XTickLabel',[]); end
        if c==1, ylabel(ax,'d\itX\rm/d\itt'); else, set(ax,'YTickLabel',[]); end
    end
    put_rot(L,5,y0+(2-r)*(sz+gy)+sz/2,sprintf('%s, %s',labels{rows(r)},etatex{rows(r)}),10);
end
put_left(L,x0,8,'gray: first 150 periods;  color: periods 2701-3000, by initial energy;  dashed: separatrix',8);
put_left(L,x0,3.5,'\itf\rm = 0.05, absorber frozen; \it\zeta\rm_1 = 0 is the conservative limit',8);
put_text(L,W/2,113,'Dissipation destroys the KAM tori: invariant curves collapse onto attractors',12);
export_figure(fig,fullfile(out,'K6_amortecimento'),dpi);close(fig);
end

% ------------------------------------------------------------------ helpers
function v=S_chaos(), v=10; end
function u=potential(x,eta), u=1.55*x.^2-3*sqrt((1.5*eta)^2+x.^2); end

function reg=fli_panel(ax,C,M,xs,vs,cmap)
imagesc(ax,xs,vs,min(max(M,0),20));set(ax,'YDir','normal');colormap(ax,cmap);caxis(ax,[0 20]);hold(ax,'on');
[XG,VG]=meshgrid(xs,vs);HG=.5*VG.^2+potential(XG,C.eta)-C.Umin;
if C.depth>0, contour(ax,XG,VG,HG,C.depth*[1 1],'LineColor','w','LineWidth',.4,'LineStyle','--'); end
contour(ax,XG,VG,HG,.5*[1 1],'LineColor','w','LineWidth',.35,'LineStyle',':');
xlim(ax,[xs(1) xs(end)]);ylim(ax,[vs(1) vs(end)]);style_axes(ax,8);pbaspect(ax,[1 1 1]);
reg=mean(M(HG<=.5)<=8);
end

function cmap=energy_colormap(n)
cmap=interp_colors({'2b1a6f','3b5bdb','1c9fd6','20b2aa','6cc24a','e0c300','f08c00','d7263d'},n);
end
function cmap=fli_colormap(n)
cmap=interp_colors({'0b1d3a','1f4e8c','3a8fb7','9ad0c2','f4e285','f4a259','bc4b51','5b1a18'},n);
end
function cmap=interp_colors(hex,n)
rgb=zeros(numel(hex),3);
for i=1:numel(hex), rgb(i,:)=sscanf(hex{i},'%2x%2x%2x',[1 3])/255; end
cmap=interp1(linspace(0,1,numel(hex)),rgb,linspace(0,1,n));
end
function c=qcolor(q)
switch q
    case 1, c=[.84 .15 .24];
    case 2, c=[.94 .55 0];
    case 3, c=[.18 .62 .27];
    otherwise, c=[.11 .49 .84];
end
end

function manual_colorbar(ax,cmap,ticks,fmt,labels)
n=size(cmap,1);v=((1:n)-.5)/n;
imagesc(ax,[v(1) v(end)],[.25 .75],[v;v]);colormap(ax,cmap);caxis(ax,[0 1]);
xlim(ax,[0 1]);ylim(ax,[0 1]);style_axes(ax,8);
set(ax,'YTick',[],'XTick',ticks,'TickLength',[.004 .004]);
if nargin<5, xtickformat(ax,fmt); else, set(ax,'XTickLabel',compose(fmt,labels(:).')); end
end

function [fig,L]=new_figure(W,H)
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 W/10 H/10], ...
    'PaperUnits','centimeters','PaperPosition',[0 0 W/10 H/10],'PaperSize',[W/10 H/10], ...
    'PaperPositionMode','manual','InvertHardcopy','off','Renderer','painters');
L=axes('Parent',fig,'Units','normalized','Position',[0 0 1 1],'XLim',[0 W],'YLim',[0 H],'Visible','off');
end
function ax=mm_axes(fig,W,H,x,y,w,h)
ax=axes('Parent',fig,'Units','normalized','Position',[x/W y/H w/W h/H]);
end
function style_axes(ax,fs)
set(ax,'FontName','Times New Roman','FontSize',fs,'LineWidth',.6,'Color','none','Box','on','Layer','top', ...
    'TickDir','in','TickLength',[.012 .012],'TickLabelInterpreter','tex');
end
function put_text(L,x,y,str,fs)
text(L,x,y,str,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Times New Roman', ...
    'FontSize',fs,'Interpreter','tex','Clipping','off');
end
function put_left(L,x,y,str,fs)
text(L,x,y,str,'HorizontalAlignment','left','VerticalAlignment','middle','FontName','Times New Roman', ...
    'FontSize',fs,'Interpreter','tex','Clipping','off');
end
function put_rot(L,x,y,str,fs)
text(L,x,y,str,'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90, ...
    'FontName','Times New Roman','FontSize',fs,'Interpreter','tex','Clipping','off');
end

function export_figure(fig,base,dpi)
% R2019a nao exporta PNG transparente diretamente: a figura e renderizada sobre
% branco e sobre preto e o canal alfa e recuperado (W = a*C + 1 - a, B = a*C).
whitefile=[tempname '.png'];blackfile=[tempname '.png'];
clean=onCleanup(@() remove_temp(whitefile,blackfile));
set(fig,'Color','w');drawnow;
print(fig,whitefile,'-dpng','-painters',sprintf('-r%d',dpi));
set(fig,'Color','k');drawnow;
print(fig,blackfile,'-dpng','-painters',sprintf('-r%d',dpi));
Wi=imread(whitefile);Bi=imread(blackfile);
assert(isequal(size(Wi),size(Bi)) && size(Wi,3)==3,'Renderizacoes incompativeis.');
alpha=zeros(size(Wi,1),size(Wi,2),'uint8');
for a=1:256:size(Wi,1)
    rows=a:min(a+255,size(Wi,1));
    wc=single(Wi(rows,:,:))/255;bc=single(Bi(rows,:,:))/255;
    ac=max(0,min(1,1-mean(wc-bc,3)));
    rgb=bsxfun(@rdivide,bc,max(ac,1/255));
    Bi(rows,:,:)=uint8(255*max(0,min(1,rgb)));
    alpha(rows,:)=uint8(255*ac);
end
imwrite(Bi,sprintf('%s_%ddpi.png',base,dpi),'Alpha',alpha,'ResolutionUnit','meter', ...
    'XResolution',round(dpi/0.0254),'YResolution',round(dpi/0.0254));
clear Wi Bi alpha;
set(fig,'Color','w');
print(fig,[base '.pdf'],'-dpdf','-painters',sprintf('-r%d',min(dpi,600)));
print(fig,[base '_preview.png'],'-dpng','-painters','-r170');
fprintf('Exportado: %s (%d dpi)\n',base,dpi);
end
function remove_temp(a,b)
if exist(a,'file'), delete(a); end
if exist(b,'file'), delete(b); end
end
