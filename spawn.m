function [y1,y2] = spawn(n,x1,x2)

% Convert integer representation of parents into n-bit binary, choose
% crossover point, and split, rearrange to get binary children

for i = 1:size(x1,2);

    x1_bin = bitget(x1(i),n:-1:1);
    x2_bin = bitget(x2(i),n:-1:1);

    a = randperm(n);
    split_pt = a(1);

    y_01 = x1_bin(1:split_pt-1);
    y_11 = x2_bin(split_pt:n);

    y_02 = x2_bin(1:split_pt-1);
    y_12 = x1_bin(split_pt:n);

    y1_bin = [y_01 y_11];
    y2_bin = [y_02 y_12];
%--------------------

% Form test array, compare bit to bit with each child, decide whether
% to mutate each bit, perform xor to mutate

    m_index1 = rand(1,n);
    m_index2 = rand(1,n);
    for s = 1:n;
        if ((m_index1(s)<0.65)&&(m_index1(s)>0.5))
        m_index1(s)= 1;
        else
        m_index1(s)= 0;
        end
    end

    for s = 1:n;
        if ((m_index2(s)<0.65)&(m_index2(s)>0.50))
        m_index2(s)= 1;
        else
        m_index2(s)= 0;
        end
    end
    y1_bin_1 = xor(y1_bin,m_index1);
    y2_bin_1 = xor(y2_bin,m_index2);
%--------------------

% Convert binary to integer and return integer children to calling program.
yy1 = 0;
yy2 = 0;

    for s = 1:n;
        yy1 = yy1 + (2^(n-s))*y1_bin_1(s);
        yy2 = yy2 + (2^(n-s))*y2_bin_1(s);
    end
    y1(i) = yy1;
    y2(i) = yy2;
end

%---------------------
