function f = v_ind(vector,r)
a=[];
    for i=1:r
            a=[a sum(vector==i)]; 
    end
f=a;
end