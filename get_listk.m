function[listk] = get_listk(l)

listk = zeros(l+1,1);
counter = 1;

for i=-l:1:l
    if mod(i+l,2) == 0
        listk(counter) = i;
        counter = counter + 1;
    end
end

end