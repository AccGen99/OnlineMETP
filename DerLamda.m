function val = DerLamda(lamda, tf)
    t_arr = linspace(0, tf, length(lamda));
    t_step = t_arr(2) - t_arr(1);
    range = 1:length(t_arr);
    val = 0*lamda;
    
    for j=range
        if j == length(range)
            val(j,:) = lamda(j,:);%val(j-1);
        else
            val(j,:) = (lamda(j,:)-lamda(j+1,:))/t_step;
        end
    end
end