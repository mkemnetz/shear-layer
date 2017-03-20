%%
newF14_data = F14_data;
newF14_data(:, 3) = (F14_data(:, 3)-1.9)./2;

figure(); 
plot(newF14_data(163624:164500, :));
%%
newF14_data2 = zeros(1, 1);
k = 1;
for i = 163628:1125000
    if (newF14_data(i, 3) > 0 && newF14_data(i-1, 3) <0)
        newF14_data2(k) = newF14_data(i, 1);
        k = k+1;
    end
end