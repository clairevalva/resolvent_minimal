% last edited by CV, October 2023

%%%%
% file to plot Rossler system approximate eigenfunctions
%
% set eigplt to an array of indices of the eigenfunctions to plot. Recommended to use no more than 4(?) plots due to an aliasing issue when plotting on the attractor if there are too many scatterplots
%%%%

eigplt = pltTiles; %[4, 28, 54], index of eigenfunctions to plot
timeplot_short = 40; % length of time to plot if freq > 0.6
timeplot_long = 160; % length of time to plot if freq < 0.6, otherwise might not get several full cycles

% get data / time coords
toplot = zeta; 
re_toplot = real(toplot);
lenplot = ceil(timeplot_short/delt);
x  =  getSrcData( model );
plottimes = (1:size(toplot,1))*delt;
ncols = size(eigplt,2);

lenattractor = 62000;
if size(x,2) < lenattractor
    lenattractor = 6200;
end

% name figure
phistr = string(eigplt(1));
for i  = 2:ncols
    phistr = phistr + "_" + string(eigplt(i));
end
savename = "figs/" + savestart + "/"  + namecon + "_allplt_" + phistr  + "_tiled.png";

% get ok layout
setfontsize = 16;
setlw = 1.5;
markerSize = 7; 

h_top = 4;
h_bot = 5;
totrows = 2*h_top + h_bot;
fscale = 0.8

clf 
figure('position', [1, 1, 500*fscale*ncols + 60, 1040*fscale])
t = tiledlayout(totrows, ncols);
t.TileSpacing = 'compact';

t_txt = title(t,'Rossler system')
t_txt.FontSize = 20
for iRow = 1:3
    for iCol = 1:ncols
        k = eigplt(iCol);

        if abs(frequencies(k)) < 0.6
            lenplotuse = ceil(timeplot_long/delt);
        else
            lenplotuse = lenplot
        end

        if iRow == 1
            nexttile([h_top 1])
            plot((plottimes(1:lenplotuse)), real(toplot(1:lenplotuse,k)), 'LineWidth',setlw);
            ylim([-2,2])
            test2 = sprintf( '  \\psi_{%i}', k ) + ", \omega: " + string(abs(frequencies(k)))
            title( test2 )
            xlabel("time")
            if iCol == 1
                ylabel(sprintf('re(\\psi)'));
            end
    
        elseif iRow == 2
            nexttile([h_top 1])
            plot(real(toplot(1:lenplotuse,k)), imag(toplot(1:lenplotuse,k)), 'LineWidth',setlw);
            xlabel(sprintf('re(\\psi)'))
            ylabel(sprintf('im(\\psi)'))
            axis equal

        elseif iRow == 3
            nexttile([h_bot 1])
            scatter3( x( 1, 1:(lenattractor) ), x( 2, 1:(lenattractor) ), x( 3, 1:(lenattractor) ), markerSize, re_toplot( 1:lenattractor, k ), ...
          'filled'  )
            set( gca, 'cLim', [ -1 1 ] *1.3 )
            axis off

            if iCol == 1
                cb = colorbar
                cb.Layout.Tile = 'east';
                cb.Label.String = sprintf('re(\\psi)');
            end
        end

        ax = gca; 
        ax.FontSize = setfontsize; 
    end
end

print('-dpng','-r400',savename)
