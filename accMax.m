function result = accMax()
    u_max = [0, -1, 1]';
    metvParams;
    acc = k2*inv(R')*inv(A)*(B')*u_max;
    result = acc(1);
end