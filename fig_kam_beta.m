function fig_kam_beta(dpi)
% QZS-ADV: mapas KAM (FLI 2.5 GL conservativo) versus a sintonia beta do absorvedor.
% Uso: fig_kam_beta            % 1900 dpi
%      fig_kam_beta(150)       % previa rapida
% Le dados_kam_beta.mat (exportado de work/compute_kam_beta_mat.py), com os mapas FLI
% do sistema completo 2.5 GL no limite conservativo (zeta=0), f=0.05, para tres
% sintonias: beta=0.35 (baseline), 0.70 (vale de caos) e 1.00 (=Omega, antirressonancia).
% Linhas = casos (monostable, shallow, deep); colunas = beta.
% Compativel com MATLAB R2019a+; nenhuma toolbox adicional.
% PNG: fundo transparente. PDF: vetorial (mapas embutidos como imagem).

if nargin<1, dpi=1900; end
validateattributes(dpi,{'numeric'},{'scalar','integer','positive'});
root=fileparts(mfilename('fullpath'));
S=load(fullfile(root,'dados_kam_beta.mat'));B=S.B;
names={'monostable','shallow_wells','deep_wells'};
labels={'Monostable, \it\eta\rm = \it\eta\rm_{QZS}','Shallow wells, \it\eta\rm = 0.60','Deep wells, \it\eta\rm = 0.30'};
betas=B.betas(:).';
sub={'baseline','vale de caos','\it\beta\rm = \it\Omega\rm (antirressonancia)'};
out=fullfile(root,sprintf('figuras_matlab_%ddpi',dpi));
if ~exist(out,'dir'), mkdir(out); end

W=180;H=176;[fig,L]=new_figure(W,H);sz=42;x0=22;gx=6;gy=12;y0=16;cmap=fli_colormap(512);
for r=1:3
    C=B.(names{r});xs=C.xs(:).';vs=C.vs(:).';eta=C.eta;Um=C.Umin;depth=C.depth;
    [XG,VG]=meshgrid(xs,vs);HG=.5*VG.^2+potential(XG,eta)-Um;
    for c=1:numel(betas)
        ax=mm_axes(fig,W,H,x0+(c-1)*(sz+gx),y0+(3-r)*(sz+gy),sz,sz);hold(ax,'on');
        M=squeeze(C.maps(c,:,:));
        imagesc(ax,xs,vs,min(max(M,0),20));set(ax,'YDir','normal');colormap(ax,cmap);caxis(ax,[0 20]);
        if depth>0, contour(ax,XG,VG,HG,depth*[1 1],'LineColor','w','LineWidth',.4,'LineStyle','--'); end
        contour(ax,XG,VG,HG,.5*[1 1],'LineColor','w','LineWidth',.35,'LineStyle',':');
        xlim(ax,[xs(1) xs(end)]);ylim(ax,[vs(1) vs(end)]);style_axes(ax,8);pbaspect(ax,[1 1 1]);
        text(ax,.03,.97,sprintf('regular %.0f%%',100*C.reg(c)),'Units','normalized','VerticalAlignment','top', ...
            'FontName','Times New Roman','FontSize',8,'Color','w');
        if r==1
            put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+3*sz+2*gy+7,sprintf('\\it\\beta\\rm = %.2f',betas(c)),11);
            put_text(L,x0+(c-1)*(sz+gx)+sz/2,y0+3*sz+2*gy+2.5,sub{c},8.5);
        end
        if r==3, xlabel(ax,'\itX\rm_0'); else, set(ax,'XTickLabel',[]); end
        if c==1, ylabel(ax,'d\itX\rm_0/d\itt'); else, set(ax,'YTickLabel',[]); end
    end
    put_rot(L,5,y0+(3-r)*(sz+gy)+sz/2,labels{r},9.5);
end
cb=mm_axes(fig,W,H,x0,7,52,2);manual_colorbar(cb,cmap,[0 5 10 15 20]/20,'%g',[0 5 10 15 20]);
put_left(L,x0+54,8,'FLI apos 400 periodos (escuro: regularidade; claro: caos)',8);
put_left(L,x0,3,['sistema completo 2.5 GL (com absorvedor), limite conservativo \it\zeta\rm = 0, \itf\rm = 0.05, \it\mu\rm = 0.1;' ...
    '  branco pontilhado: \itH\rm_0 - \itU\rm_{min} = 0.5'],8);
put_text(L,W/2,171,'Regularidade (indicador FLI) do QZS-ADV versus a sintonia \it\beta\rm',12.5);
export_figure(fig,fullfile(out,'KAM_sintonia_beta'),dpi);close(fig);
fprintf('Figura salva em: %s\n',out);
end

% ------------------------------------------------------------------ helpers
function u=potential(x,eta), u=1.55*x.^2-3*sqrt((1.5*eta)^2+x.^2); end

function cmap=fli_colormap(n)
cmap=interp_colors({'0b1d3a','1f4e8c','3a8fb7','9ad0c2','f4e285','f4a259','bc4b51','5b1a18'},n);
end
function cmap=interp_colors(hex,n)
rgb=zeros(numel(hex),3);
for i=1:numel(hex), rgb(i,:)=sscanf(hex{i},'%2x%2x%2x',[1 3])/255; end
cmap=interp1(linspace(0,1,numel(hex)),rgb,linspace(0,1,n));
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
% R2019a nao exporta PNG transparente diretamente: renderiza sobre branco e sobre preto
% e recupera o canal alfa (W = a*C + 1 - a, B = a*C).
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
