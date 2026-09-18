function custom_colormap = osclormap1new()
%OSCOLORMAP1NEW Return colors for the four oxygen-flux components.
rx=[0.78 0.66 0.63]; % Zonal lateral induction
ry=[0.76 0.79 0.46]; % Meridional lateral induction
g=[1.00 0.83 0.45];  % Vertical velocity
y=[0.47 0.80 0.88];  % Eddy-induced subduction
custom_colormap=[rx;ry;g;y];
end

