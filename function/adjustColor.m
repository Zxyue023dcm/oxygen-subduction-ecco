function adjusted_rgb = adjustColor(rgb, s_scale, v_scale)
%ADJUSTCOLOR Scale the saturation and value of an RGB color.
%   RGB is a three-element vector in [0, 1].

hsv = rgb2hsv(rgb);
hsv(2) = min(1, max(0, hsv(2) * s_scale));
hsv(3) = min(1, max(0, hsv(3) * v_scale));
adjusted_rgb = hsv2rgb(hsv);
end

