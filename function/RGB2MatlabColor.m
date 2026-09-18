function c = RGB2MatlabColor(rgb)
%RGB2MATLABCOLOR : RGB值转化为Matlab颜色向量
%   此处显示详细说明
c = round(100 * rgb/255) / 100;
end

