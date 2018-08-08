load Chn.mat;
Num_Cell = 7;
B2Bdist = 0.8;
BSLoc = 0.5*B2Bdist* [0; sqrt(3)+j; 2*j; -sqrt(3)+j; -sqrt(3)-j; -2*j; sqrt(3)-j;
    sqrt(12); sqrt(12)+2*j; sqrt(3)+3*j; 4*j; -sqrt(3)+3*j; -sqrt(12)+2*j; -sqrt(12);
    -sqrt(12)-2*j; -sqrt(3)-3*j; -4*j; sqrt(3)-3*j; sqrt(12)-2*j];

PicoLoc = 0.5*B2Bdist*[1/sqrt(3); (cos(2*pi/3)+j*sin(2*pi/3))/sqrt(3); (cos(2*pi/3)+j*sin(-2*pi/3))/sqrt(3)];

AllBSLoc = zeros(Num_Cell,Num_BS);

for l = 1:Num_Cell
    for n = 1:Num_BS
        if n == 1
            AllBSLoc(l,n) = BSLoc(l);
        else
            AllBSLoc(l,n) = BSLoc(l) + PicoLoc(n-1);
        end
        
    end
end

figure;
MacroBS_Loc = AllBSLoc(:,1);
PicoBS_Loc = AllBSLoc(:,[2:end]);

scatter(real(MacroBS_Loc(:)),imag(MacroBS_Loc(:)),'rs');

hold on;

scatter(real(PicoBS_Loc(:)),imag(PicoBS_Loc(:)),'r*');

scatter(real(MULoc(:)),imag(MULoc(:)),'bo');

cell_edge = B2Bdist/sqrt(3)*exp(j*[0:pi/3:2*pi]);

for l = 1:Num_Cell
    cell_edge_point = cell_edge + BSLoc(l);
    plot(cell_edge_point,'k');
end


grid on;
legend('Macro BS', 'Pico BS','Mobile User');
saveas(gcf,'HetNetTopology.fig','fig');